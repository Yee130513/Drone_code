#include <cstdlib>
#include <fstream>
#include <plan_manage/planner_manager.h>
#include <string>
#include <thread>
#include "visualization_msgs/Marker.h" // zx-todo

namespace ego_planner
{

  // SECTION interfaces for setup and query

  EGOPlannerManager::EGOPlannerManager() {}

  EGOPlannerManager::~EGOPlannerManager() {}

  void EGOPlannerManager::initPlanModules(ros::NodeHandle &nh, PlanningVisualization::Ptr vis)
  {
    /* read algorithm parameters */

    nh.param("manager/max_vel", pp_.max_vel_, -1.0);
    nh.param("manager/max_acc", pp_.max_acc_, -1.0);
    nh.param("manager/max_jerk", pp_.max_jerk_, -1.0);
    nh.param("manager/feasibility_tolerance", pp_.feasibility_tolerance_, 0.0);
    nh.param("manager/control_points_distance", pp_.ctrl_pt_dist, -1.0);
    nh.param("manager/planning_horizon", pp_.planning_horizen_, 5.0);
    nh.param("manager/use_distinctive_trajs", pp_.use_distinctive_trajs, false);
    nh.param("manager/drone_id", pp_.drone_id, -1);
    nh.param("manager/time_scaling_method", time_scaling_method_, std::string("Global"));
    /* read algorithm parameters锛屽弬鏁拌鍙?*/
    nh.param("manager/max_vel", pp_.max_vel_, -1.0);
    // ... 鍘熸湁鍙傛暟璇诲彇 ...
    nh.param("manager/drone_id", pp_.drone_id, -1);

    // [鏂板] 璇诲彇娑堣瀺瀹為獙鐨勭畻娉曞悕绉帮紝榛樿涓?"A_star"
    nh.param("manager/planner_type", planner_type_, std::string("A_star"));

    local_data_.traj_id_ = 0;
    grid_map_.reset(new GridMap);
    grid_map_->initMap(nh);

    // obj_predictor_.reset(new fast_planner::ObjPredictor(nh));
    // obj_predictor_->init();
    // obj_pub_ = nh.advertise<visualization_msgs::Marker>("/dynamic/obj_prdi", 10); // zx-todo

    bspline_optimizer_.reset(new BsplineOptimizer);
    bspline_optimizer_->setParam(nh);
    bspline_optimizer_->setEnvironment(grid_map_, obj_predictor_);
    //bspline_optimizer_->a_star_.reset(new ThetaStar);
    //bspline_optimizer_->a_star_.reset(new JPSThetaPlanner);
    bspline_optimizer_->a_star_.reset(new AStar);
    bspline_optimizer_->a_star_->setPlannerType(planner_type_);
    bspline_optimizer_->a_star_->initGridMap(grid_map_, Eigen::Vector3i(100, 100, 100));

    visualization_ = vis;
  }

  // !SECTION

  // SECTION rebond replanning

  bool EGOPlannerManager::reboundReplan(Eigen::Vector3d start_pt, Eigen::Vector3d start_vel,
                                        Eigen::Vector3d start_acc, Eigen::Vector3d local_target_pt,
                                        Eigen::Vector3d local_target_vel, bool flag_polyInit, bool flag_randomPolyTraj)
  {
    static int count = 0;
    printf("\033[47;30m\n[drone %d replan %d]==============================================\033[0m\n", pp_.drone_id, count++);
    // cout.precision(3);
    // cout << "start: " << start_pt.transpose() << ", " << start_vel.transpose() << "\ngoal:" << local_target_pt.transpose() << ", " << local_target_vel.transpose()
    //      << endl;

    if ((start_pt - local_target_pt).norm() < 0.2)
    {
      cout << "Close to goal" << endl;
      continous_failures_count_++;
      return false;
    }

    bspline_optimizer_->setLocalTargetPt(local_target_pt);

    ros::Time t_start = ros::Time::now();
    ros::Duration t_init, t_opt, t_refine;

    /*** STEP 1: INIT ***/
    double ts = (start_pt - local_target_pt).norm() > 0.1 ? pp_.ctrl_pt_dist / pp_.max_vel_ * 1.5 : pp_.ctrl_pt_dist / pp_.max_vel_ * 5; // pp_.ctrl_pt_dist / pp_.max_vel_ is too tense, and will surely exceed the acc/vel limits
    vector<Eigen::Vector3d> point_set, start_end_derivatives;
    static bool flag_first_call = true, flag_force_polynomial = false;
    bool flag_regenerate = false;
    do
    {
      point_set.clear();
      start_end_derivatives.clear();
      flag_regenerate = false;

      if (flag_first_call || flag_polyInit || flag_force_polynomial /*|| ( start_pt - local_target_pt ).norm() < 1.0*/) // Initial path generated from a min-snap traj by order.
      {
        flag_first_call = false;
        flag_force_polynomial = false;

        PolynomialTraj gl_traj;

        double dist = (start_pt - local_target_pt).norm();
        double time = pow(pp_.max_vel_, 2) / pp_.max_acc_ > dist ? sqrt(dist / pp_.max_acc_) : (dist - pow(pp_.max_vel_, 2) / pp_.max_acc_) / pp_.max_vel_ + 2 * pp_.max_vel_ / pp_.max_acc_;

        if (!flag_randomPolyTraj)
        {
          gl_traj = PolynomialTraj::one_segment_traj_gen(start_pt, start_vel, start_acc, local_target_pt, local_target_vel, Eigen::Vector3d::Zero(), time);
        }
        else
        {
          Eigen::Vector3d horizen_dir = ((start_pt - local_target_pt).cross(Eigen::Vector3d(0, 0, 1))).normalized();
          Eigen::Vector3d vertical_dir = ((start_pt - local_target_pt).cross(horizen_dir)).normalized();
          Eigen::Vector3d random_inserted_pt = (start_pt + local_target_pt) / 2 +
                                               (((double)rand()) / RAND_MAX - 0.5) * (start_pt - local_target_pt).norm() * horizen_dir * 0.8 * (-0.978 / (continous_failures_count_ + 0.989) + 0.989) +
                                               (((double)rand()) / RAND_MAX - 0.5) * (start_pt - local_target_pt).norm() * vertical_dir * 0.4 * (-0.978 / (continous_failures_count_ + 0.989) + 0.989);
          Eigen::MatrixXd pos(3, 3);
          pos.col(0) = start_pt;
          pos.col(1) = random_inserted_pt;
          pos.col(2) = local_target_pt;
          Eigen::VectorXd t(2);
          t(0) = t(1) = time / 2;
          gl_traj = PolynomialTraj::minSnapTraj(pos, start_vel, local_target_vel, start_acc, Eigen::Vector3d::Zero(), t);
        }

        double t;
        bool flag_too_far;
        ts *= 1.5; // ts will be divided by 1.5 in the next
        do
        {
          ts /= 1.5;
          point_set.clear();
          flag_too_far = false;
          Eigen::Vector3d last_pt = gl_traj.evaluate(0);
          for (t = 0; t < time; t += ts)
          {
            Eigen::Vector3d pt = gl_traj.evaluate(t);
            if ((last_pt - pt).norm() > pp_.ctrl_pt_dist * 1.5)
            {
              flag_too_far = true;
              break;
            }
            last_pt = pt;
            point_set.push_back(pt);
          }
        } while (flag_too_far || point_set.size() < 7); // To make sure the initial path has enough points.
        t -= ts;
        start_end_derivatives.push_back(gl_traj.evaluateVel(0));
        start_end_derivatives.push_back(local_target_vel);
        start_end_derivatives.push_back(gl_traj.evaluateAcc(0));
        start_end_derivatives.push_back(gl_traj.evaluateAcc(t));
      }
      else // Initial path generated from previous trajectory.
      {

        double t;
        double t_cur = (ros::Time::now() - local_data_.start_time_).toSec();

        vector<double> pseudo_arc_length;
        vector<Eigen::Vector3d> segment_point;
        pseudo_arc_length.push_back(0.0);
        for (t = t_cur; t < local_data_.duration_ + 1e-3; t += ts)
        {
          segment_point.push_back(local_data_.position_traj_.evaluateDeBoorT(t));
          if (t > t_cur)
          {
            pseudo_arc_length.push_back((segment_point.back() - segment_point[segment_point.size() - 2]).norm() + pseudo_arc_length.back());
          }
        }
        t -= ts;

        double poly_time = (local_data_.position_traj_.evaluateDeBoorT(t) - local_target_pt).norm() / pp_.max_vel_ * 2;
        if (poly_time > ts)
        {
          PolynomialTraj gl_traj = PolynomialTraj::one_segment_traj_gen(local_data_.position_traj_.evaluateDeBoorT(t),
                                                                        local_data_.velocity_traj_.evaluateDeBoorT(t),
                                                                        local_data_.acceleration_traj_.evaluateDeBoorT(t),
                                                                        local_target_pt, local_target_vel, Eigen::Vector3d::Zero(), poly_time);

          for (t = ts; t < poly_time; t += ts)
          {
            if (!pseudo_arc_length.empty())
            {
              segment_point.push_back(gl_traj.evaluate(t));
              pseudo_arc_length.push_back((segment_point.back() - segment_point[segment_point.size() - 2]).norm() + pseudo_arc_length.back());
            }
            else
            {
              ROS_ERROR("pseudo_arc_length is empty, return!");
              continous_failures_count_++;
              return false;
            }
          }
        }

        double sample_length = 0;
        double cps_dist = pp_.ctrl_pt_dist * 1.5; // cps_dist will be divided by 1.5 in the next
        size_t id = 0;
        do
        {
          cps_dist /= 1.5;
          point_set.clear();
          sample_length = 0;
          id = 0;
          while ((id <= pseudo_arc_length.size() - 2) && sample_length <= pseudo_arc_length.back())
          {
            if (sample_length >= pseudo_arc_length[id] && sample_length < pseudo_arc_length[id + 1])
            {
              point_set.push_back((sample_length - pseudo_arc_length[id]) / (pseudo_arc_length[id + 1] - pseudo_arc_length[id]) * segment_point[id + 1] +
                                  (pseudo_arc_length[id + 1] - sample_length) / (pseudo_arc_length[id + 1] - pseudo_arc_length[id]) * segment_point[id]);
              sample_length += cps_dist;
            }
            else
              id++;
          }
          point_set.push_back(local_target_pt);
        } while (point_set.size() < 7); // If the start point is very close to end point, this will help

        start_end_derivatives.push_back(local_data_.velocity_traj_.evaluateDeBoorT(t_cur));
        start_end_derivatives.push_back(local_target_vel);
        start_end_derivatives.push_back(local_data_.acceleration_traj_.evaluateDeBoorT(t_cur));
        start_end_derivatives.push_back(Eigen::Vector3d::Zero());

        if (point_set.size() > pp_.planning_horizen_ / pp_.ctrl_pt_dist * 3) // The initial path is unnormally too long!
        {
          flag_force_polynomial = true;
          flag_regenerate = true;
        }
      }
    } while (flag_regenerate);

    Eigen::MatrixXd ctrl_pts, ctrl_pts_temp;
    UniformBspline::parameterizeToBspline(ts, point_set, start_end_derivatives, ctrl_pts);

    /* ================= ABLATION STUDY: DATA LOGGER ================= */
    auto escape_csv = [](const std::string &s) -> std::string
    {
      std::string out = s;
      for (char &c : out)
      {
        if (c == ',')
        {
          c = ';';
        }
      }
      return out;
    };

    auto log_ablation_data = [&](bool optimization_success, double jerk_cost)
    {
      const char *home_env = std::getenv("HOME");
      if (home_env == nullptr)
      {
        ROS_WARN_THROTTLE(1.0, "[DataLogger] HOME env is null, skip logging.");
        return;
      }

      const std::string file_name = std::string(home_env) + "/.ros/sca_ablation.csv";

      std::ofstream file(file_name, std::ios::out | std::ios::app);
      if (!file.is_open())
      {
        ROS_WARN_THROTTLE(1.0, "[DataLogger] failed to open %s", file_name.c_str());
        return;
      }

      file.seekp(0, std::ios::end);
      if (file.tellp() == 0)
      {
        file << "Timestamp,"
             << "Planner_Type,"
             << "Init_Time_ms,"
             << "Search_Time_ms,"
             << "Optimization_Time_ms,"
             << "Nodes_Expanded,"
             << "Path_Length_m,"
             << "Minimum_Clearance_m,"
             << "Segment_Recheck_Fail_Count,"
             << "Trajectory_Jerk_Cost,"
             << "Optimization_Success,"
             << "Opt_Iterations,"
             << "AStar_Used_Fallback,"
             << "AStar_Fallback_Reason,"
             << "Guide_Path_Count,"
             << "Guide_Invalid_Count,"
             << "Init_Collision_Pts,"
             << "Frontend_Failed,"
             << "Frontend_Fail_Reason,"
             << "Dynamic_Clearance_Enabled,"
             << "Dynamic_Clearance_Count,"
             << "Dynamic_Clearance_Fallback_Count,"
             << "Dynamic_Clearance_Min_m,"
             << "Dynamic_Clearance_Mean_m,"
             << "Dynamic_Clearance_Max_m\n";
      }

      int nodes_expanded = 0;
      double path_length_m = 0.0;
      double search_time_ms = 0.0;
      double minimum_clearance_m = -1.0;
      int segment_recheck_fail_count = 0;
      int astar_used_fallback = 0;
      std::string astar_fallback_reason = "NO_COLLISION";

      int opt_iterations = -1;
      int guide_path_count = -1;
      int guide_invalid_count = -1;
      int init_collision_pts = -1;
      int frontend_failed = -1;
      std::string frontend_fail_reason = "NA";
      int dyn_clearance_enabled = 0;
      int dyn_clearance_count = 0;
      int dyn_clearance_fallback_count = 0;
      double dyn_clearance_min_m = 0.0;
      double dyn_clearance_mean_m = 0.0;
      double dyn_clearance_max_m = 0.0;

      if (bspline_optimizer_ != nullptr)
      {
        opt_iterations = bspline_optimizer_->getIterNum();
        guide_path_count = bspline_optimizer_->getLastGuidePathCount();
        guide_invalid_count = bspline_optimizer_->getLastGuideInvalidCount();
        init_collision_pts = bspline_optimizer_->getInitCollisionPts();
        frontend_failed = bspline_optimizer_->getLastFrontendFailed() ? 1 : 0;
        frontend_fail_reason = bspline_optimizer_->getLastFrontendFailReason();

        dyn_clearance_enabled = bspline_optimizer_->isDynamicClearanceEnabled() ? 1 : 0;
        dyn_clearance_count = bspline_optimizer_->getDynamicClearanceCount();
        dyn_clearance_fallback_count = bspline_optimizer_->getDynamicClearanceFallbackCount();
        dyn_clearance_min_m = bspline_optimizer_->getDynamicClearanceMin();
        dyn_clearance_mean_m = bspline_optimizer_->getDynamicClearanceMean();
        dyn_clearance_max_m = bspline_optimizer_->getDynamicClearanceMax();

        if (bspline_optimizer_->a_star_ != nullptr && (guide_path_count > 0 || frontend_failed == 1))
        {
          nodes_expanded = bspline_optimizer_->a_star_->getExpandedNodes();
          path_length_m = bspline_optimizer_->a_star_->getPathLength();
          search_time_ms = bspline_optimizer_->a_star_->getLastSearchTimeMs();
          minimum_clearance_m = bspline_optimizer_->a_star_->getLastMinClearance();
          segment_recheck_fail_count = bspline_optimizer_->a_star_->getSegmentCheckFailCount();
          astar_used_fallback = bspline_optimizer_->a_star_->usedFallback() ? 1 : 0;
          astar_fallback_reason = bspline_optimizer_->a_star_->getLastFallbackReason();
        }
      }

      file << ros::Time::now().toSec() << ","
           << escape_csv(planner_type_) << ","
           << t_init.toSec() * 1000.0 << ","
           << search_time_ms << ","
           << t_opt.toSec() * 1000.0 << ","
           << nodes_expanded << ","
           << path_length_m << ","
           << minimum_clearance_m << ","
           << segment_recheck_fail_count << ","
           << jerk_cost << ","
           << (optimization_success ? 1 : 0) << ","
           << opt_iterations << ","
           << astar_used_fallback << ","
           << escape_csv(astar_fallback_reason) << ","
           << guide_path_count << ","
           << guide_invalid_count << ","
           << init_collision_pts << ","
           << frontend_failed << ","
           << escape_csv(frontend_fail_reason) << ","
           << dyn_clearance_enabled << ","
           << dyn_clearance_count << ","
           << dyn_clearance_fallback_count << ","
           << dyn_clearance_min_m << ","
           << dyn_clearance_mean_m << ","
           << dyn_clearance_max_m << "\n";

      file.close();
    };
    /* ============================================================== */

    vector<std::pair<int, int>> segments;
    segments = bspline_optimizer_->initControlPoints(ctrl_pts, true);

    if (bspline_optimizer_->getLastFrontendFailed())
    {
      t_init = ros::Time::now() - t_start;
      t_opt = ros::Duration(0.0);

      ROS_ERROR_THROTTLE(0.5,
                         "[PlannerManager] frontend guide failed before optimization. reason=%s",
                         bspline_optimizer_->getLastFrontendFailReason().c_str());
      continous_failures_count_++;
      log_ablation_data(false, -1.0);
      return false;
    }

    t_init = ros::Time::now() - t_start;
    t_start = ros::Time::now();

    /*** STEP 2: OPTIMIZE ***/
    bool flag_step_1_success = false;
    vector<vector<Eigen::Vector3d>> vis_trajs;

    if (pp_.use_distinctive_trajs)
    {
      // cout << "enter" << endl;
      std::vector<ControlPoints> trajs = bspline_optimizer_->distinctiveTrajs(segments);
      cout << "\033[1;33m"
           << "multi-trajs=" << trajs.size() << "\033[1;0m" << endl;

      double final_cost, min_cost = 999999.0;
      for (int i = trajs.size() - 1; i >= 0; i--)
      {
        if (bspline_optimizer_->BsplineOptimizeTrajRebound(ctrl_pts_temp, final_cost, trajs[i], ts))
        {

          cout << "traj " << trajs.size() - i << " success." << endl;

          flag_step_1_success = true;
          if (final_cost < min_cost)
          {
            min_cost = final_cost;
            ctrl_pts = ctrl_pts_temp;
          }

          // visualization
          point_set.clear();
          for (int j = 0; j < ctrl_pts_temp.cols(); j++)
          {
            point_set.push_back(ctrl_pts_temp.col(j));
          }
          vis_trajs.push_back(point_set);
        }
        else
        {
          cout << "traj " << trajs.size() - i << " failed." << endl;
        }
      }

      t_opt = ros::Time::now() - t_start;

      visualization_->displayMultiInitPathList(vis_trajs, 0.2); // This visuallization will take up several milliseconds.
    }
    else
    {
      flag_step_1_success = bspline_optimizer_->BsplineOptimizeTrajRebound(ctrl_pts, ts);
      t_opt = ros::Time::now() - t_start;
      //static int vis_id = 0;
      visualization_->displayInitPathList(point_set, 0.2, 0);
    }

    cout << "plan_success=" << flag_step_1_success << endl;
    if (!flag_step_1_success)
    {
      ROS_ERROR_THROTTLE(0.5, "[PlannerManager] all candidate trajectories failed. plan_success=0");
      visualization_->displayOptimalList(ctrl_pts, 0);
      continous_failures_count_++;
      log_ablation_data(false, -1.0); // 鍩嬬偣锛氫紭鍖栧け璐ワ紝Jerk 鏃犳晥璁颁负 -1.0
      return false;
    }

    t_start = ros::Time::now();

    UniformBspline pos = UniformBspline(ctrl_pts, 3, ts);
    pos.setPhysicalLimits(pp_.max_vel_, pp_.max_acc_, pp_.feasibility_tolerance_);

    /*** STEP 3: REFINE(RE-ALLOCATE TIME) IF NECESSARY ***/
    // Note: Only adjust time in single drone mode. But we still allow drone_0 to adjust its time profile.
    bool flag_step_2_success = true;
    if (pp_.drone_id <= 0)
    {

      double ratio;
      

      if (!pos.checkFeasibility(ratio, false))
      {
        cout << "Need to reallocate time." << endl;

        Eigen::MatrixXd optimal_control_points;
        flag_step_2_success = refineTrajAlgo(pos, start_end_derivatives, ratio, ts, optimal_control_points);
        if (flag_step_2_success)
          pos = UniformBspline(optimal_control_points, 3, ts);
      }

      if (!flag_step_2_success)
      {
        printf("\033[34mThis refined trajectory hits obstacles. It doesn't matter if appeares occasionally. But if continously appearing, Increase parameter \"lambda_fitness\".\n\033[0m");
        continous_failures_count_++;
        log_ablation_data(false, -1.0); // 鍩嬬偣锛歊efine 闃舵鎾炲澶辫触
        return false;
      }
    }
    else
    {
      static bool print_once = true;
      if (print_once)
      {
        print_once = false;
        ROS_ERROR("IN SWARM MODE, REFINE DISABLED!");
      }
    }

    t_refine = ros::Time::now() - t_start;

    // save planned results
    updateTrajInfo(pos, ros::Time::now());

    static double sum_time = 0;
    static int count_success = 0;
    sum_time += (t_init + t_opt + t_refine).toSec();


    // =====================================================================
    // 馃殌 ABLATION STUDY LOGGER: 鍔ㄥ姏瀛︿笌鏃堕棿閲嶅垎閰嶅疄楠屽煁鐐?    // =====================================================================
    if (flag_step_1_success &&  (flag_step_2_success)) {
        UniformBspline vel_traj = pos.getDerivative();
        UniformBspline acc_traj = vel_traj.getDerivative();
        
        double duration = pos.getTimeSum();
        double max_v = 0.0, max_a = 0.0, sum_v = 0.0;
        int sample_count = 0;

        // 1. 閲囨牱鎻愬彇寰鍔ㄥ姏瀛︽洸绾挎暟鎹?(姣?0.02 绉掗噰涓€涓偣)
        std::string home_dir = getenv("HOME");
        std::ofstream profile_file(home_dir + "/.ros/ego_vel_profile.csv", std::ios::out | std::ios::app);
        bool write_profile = profile_file.is_open();
        
        if (write_profile && profile_file.tellp() == 0) {
            profile_file << "Traj_ID,Time_t,Velocity,Acceleration,Method\n";
        }

        for (double t = 0; t <= duration; t += 0.02) {
            double v_norm = vel_traj.evaluateDeBoorT(t).norm();
            double a_norm = acc_traj.evaluateDeBoorT(t).norm();
            
            max_v = std::max(max_v, v_norm);
            max_a = std::max(max_a, a_norm);
            sum_v += v_norm;
            sample_count++;

            // 灏嗗眬閮ㄩ€熷害鏇茬嚎鍐欏叆鏂囦欢
            if (write_profile) {
                profile_file << local_data_.traj_id_ << "," << t << "," 
                             << v_norm << "," << a_norm << "," << time_scaling_method_ << "\n";
            }
        }
        if (write_profile) profile_file.close();

        double avg_v = sum_v / sample_count;

        // 2. 鍐欏叆瀹忚缁熻鏁版嵁
        std::ofstream macro_stats(home_dir + "/.ros/ego_macro_stats.csv", std::ios::out | std::ios::app);
        if (macro_stats.is_open()) {
            if (macro_stats.tellp() == 0) {
                macro_stats << "Timestamp,Method,Duration,Avg_Vel,Max_Vel,Max_Acc,Refine_Time_ms\n";
            }
            macro_stats << ros::Time::now().toSec() << ","
                       << time_scaling_method_ << ","
                       << duration << ","
                       << avg_v << ","
                       << max_v << ","
                       << max_a << ","
                       << t_refine.toSec() * 1000.0 << "\n";
            macro_stats.close();
        }
    }
    // =====================================================================

    count_success++;
    cout << "total time:\033[42m" << (t_init + t_opt + t_refine).toSec() << "\033[0m,optimize:" << (t_init + t_opt).toSec() << ",refine:" << t_refine.toSec() << ",avg_time=" << sum_time / count_success << endl;
    log_ablation_data(true, pos.getJerk()); 
    // success. YoY
    continous_failures_count_ = 0;
    return true;
  }

  bool EGOPlannerManager::EmergencyStop(Eigen::Vector3d stop_pos)
  {
    Eigen::MatrixXd control_points(3, 6);
    for (int i = 0; i < 6; i++)
    {
      control_points.col(i) = stop_pos;
    }

    updateTrajInfo(UniformBspline(control_points, 3, 1.0), ros::Time::now());

    return true;
  }

  bool EGOPlannerManager::checkCollision(int drone_id)
  {
    if (local_data_.start_time_.toSec() < 1e9) // It means my first planning has not started
      return false;

    double my_traj_start_time = local_data_.start_time_.toSec();
    double other_traj_start_time = swarm_trajs_buf_[drone_id].start_time_.toSec();

    double t_start = max(my_traj_start_time, other_traj_start_time);
    double t_end = min(my_traj_start_time + local_data_.duration_ * 2 / 3, other_traj_start_time + swarm_trajs_buf_[drone_id].duration_);

    for (double t = t_start; t < t_end; t += 0.03)
    {
      if ((local_data_.position_traj_.evaluateDeBoorT(t - my_traj_start_time) - swarm_trajs_buf_[drone_id].position_traj_.evaluateDeBoorT(t - other_traj_start_time)).norm() < bspline_optimizer_->getSwarmClearance())
      {
        return true;
      }
    }

    return false;
  }

  bool EGOPlannerManager::planGlobalTrajWaypoints(const Eigen::Vector3d &start_pos, const Eigen::Vector3d &start_vel, const Eigen::Vector3d &start_acc,
                                                  const std::vector<Eigen::Vector3d> &waypoints, const Eigen::Vector3d &end_vel, const Eigen::Vector3d &end_acc)
  {

    // generate global reference trajectory

    vector<Eigen::Vector3d> points;
    points.push_back(start_pos);

    for (size_t wp_i = 0; wp_i < waypoints.size(); wp_i++)
    {
      points.push_back(waypoints[wp_i]);
    }

    double total_len = 0;
    total_len += (start_pos - waypoints[0]).norm();
    for (size_t i = 0; i < waypoints.size() - 1; i++)
    {
      total_len += (waypoints[i + 1] - waypoints[i]).norm();
    }

    // insert intermediate points if too far
    vector<Eigen::Vector3d> inter_points;
    double dist_thresh = max(total_len / 8, 4.0);

    for (size_t i = 0; i < points.size() - 1; ++i)
    {
      inter_points.push_back(points.at(i));
      double dist = (points.at(i + 1) - points.at(i)).norm();

      if (dist > dist_thresh)
      {
        int id_num = floor(dist / dist_thresh) + 1;

        for (int j = 1; j < id_num; ++j)
        {
          Eigen::Vector3d inter_pt =
              points.at(i) * (1.0 - double(j) / id_num) + points.at(i + 1) * double(j) / id_num;
          inter_points.push_back(inter_pt);
        }
      }
    }

    inter_points.push_back(points.back());

    // for ( int i=0; i<inter_points.size(); i++ )
    // {
    //   cout << inter_points[i].transpose() << endl;
    // }

    // write position matrix
    int pt_num = inter_points.size();
    Eigen::MatrixXd pos(3, pt_num);
    for (int i = 0; i < pt_num; ++i)
      pos.col(i) = inter_points[i];

    Eigen::Vector3d zero(0, 0, 0);
    Eigen::VectorXd time(pt_num - 1);
    for (int i = 0; i < pt_num - 1; ++i)
    {
      time(i) = (pos.col(i + 1) - pos.col(i)).norm() / (pp_.max_vel_);
    }

    time(0) *= 2.0;
    time(time.rows() - 1) *= 2.0;

    PolynomialTraj gl_traj;
    if (pos.cols() >= 3)
      gl_traj = PolynomialTraj::minSnapTraj(pos, start_vel, end_vel, start_acc, end_acc, time);
    else if (pos.cols() == 2)
      gl_traj = PolynomialTraj::one_segment_traj_gen(start_pos, start_vel, start_acc, pos.col(1), end_vel, end_acc, time(0));
    else
      return false;

    auto time_now = ros::Time::now();
    global_data_.setGlobalTraj(gl_traj, time_now);

    return true;
  }

  bool EGOPlannerManager::planGlobalTraj(const Eigen::Vector3d &start_pos, const Eigen::Vector3d &start_vel, const Eigen::Vector3d &start_acc,
                                         const Eigen::Vector3d &end_pos, const Eigen::Vector3d &end_vel, const Eigen::Vector3d &end_acc)
  {

    // generate global reference trajectory

    vector<Eigen::Vector3d> points;
    points.push_back(start_pos);
    points.push_back(end_pos);

    // insert intermediate points if too far
    vector<Eigen::Vector3d> inter_points;
    const double dist_thresh = 4.0;

    for (size_t i = 0; i < points.size() - 1; ++i)
    {
      inter_points.push_back(points.at(i));
      double dist = (points.at(i + 1) - points.at(i)).norm();

      if (dist > dist_thresh)
      {
        int id_num = floor(dist / dist_thresh) + 1;

        for (int j = 1; j < id_num; ++j)
        {
          Eigen::Vector3d inter_pt =
              points.at(i) * (1.0 - double(j) / id_num) + points.at(i + 1) * double(j) / id_num;
          inter_points.push_back(inter_pt);
        }
      }
    }

    inter_points.push_back(points.back());

    // write position matrix
    int pt_num = inter_points.size();
    Eigen::MatrixXd pos(3, pt_num);
    for (int i = 0; i < pt_num; ++i)
      pos.col(i) = inter_points[i];

    Eigen::Vector3d zero(0, 0, 0);
    Eigen::VectorXd time(pt_num - 1);
    for (int i = 0; i < pt_num - 1; ++i)
    {
      time(i) = (pos.col(i + 1) - pos.col(i)).norm() / (pp_.max_vel_);
    }

    time(0) *= 2.0;
    time(time.rows() - 1) *= 2.0;

    PolynomialTraj gl_traj;
    if (pos.cols() >= 3)
      gl_traj = PolynomialTraj::minSnapTraj(pos, start_vel, end_vel, start_acc, end_acc, time);
    else if (pos.cols() == 2)
      gl_traj = PolynomialTraj::one_segment_traj_gen(start_pos, start_vel, start_acc, end_pos, end_vel, end_acc, time(0));
    else
      return false;

    auto time_now = ros::Time::now();
    global_data_.setGlobalTraj(gl_traj, time_now);

    return true;
  }


    bool EGOPlannerManager::refineTrajAlgo(UniformBspline &traj, vector<Eigen::Vector3d> &start_end_derivative, double ratio, double &ts, Eigen::MatrixXd &optimal_control_points)
  {
    double t_inc;

    Eigen::MatrixXd ctrl_pts; // = traj.getControlPoint()

    // Planner layer always generates a conservative uniform fallback trajectory.
    // Execution-layer speed changes are handled by traj_server EVTR.
    reparamBsplineGlobal(traj, start_end_derivative, ratio, ctrl_pts, ts, t_inc);

    traj = UniformBspline(ctrl_pts, 3, ts);

    double t_step = traj.getTimeSum() / (ctrl_pts.cols() - 3);
    bspline_optimizer_->ref_pts_.clear();
    for (double t = 0; t < traj.getTimeSum() + 1e-4; t += t_step)
      bspline_optimizer_->ref_pts_.push_back(traj.evaluateDeBoorT(t));

    bool success = bspline_optimizer_->BsplineOptimizeTrajRefine(ctrl_pts, ts, optimal_control_points);

    return success;
  }

  void EGOPlannerManager::updateTrajInfo(const UniformBspline &position_traj, const ros::Time time_now)
  {
    local_data_.start_time_ = time_now;
    local_data_.position_traj_ = position_traj;
    local_data_.velocity_traj_ = local_data_.position_traj_.getDerivative();
    local_data_.acceleration_traj_ = local_data_.velocity_traj_.getDerivative();
    local_data_.start_pos_ = local_data_.position_traj_.evaluateDeBoorT(0.0);
    local_data_.duration_ = local_data_.position_traj_.getTimeSum();
    local_data_.traj_id_ += 1;
  }

  void EGOPlannerManager::reparamBsplineGlobal(UniformBspline &bspline, vector<Eigen::Vector3d> &start_end_derivative, double ratio,
                                         Eigen::MatrixXd &ctrl_pts, double &dt, double &time_inc)
  {
    double time_origin = bspline.getTimeSum();
    int seg_num = bspline.getControlPoint().cols() - 3;
    // double length = bspline.getLength(0.1);
    // int seg_num = ceil(length / pp_.ctrl_pt_dist);

    bspline.lengthenTime(ratio);
    double duration = bspline.getTimeSum();
    dt = duration / double(seg_num);
    time_inc = duration - time_origin;

    vector<Eigen::Vector3d> point_set;
    for (double time = 0.0; time <= duration + 1e-4; time += dt)
    {
      point_set.push_back(bspline.evaluateDeBoorT(time));
    }
    UniformBspline::parameterizeToBspline(dt, point_set, start_end_derivative, ctrl_pts);
  }

// namespace ego_planner


void EGOPlannerManager::reparamBsplineElastic(UniformBspline &bspline, vector<Eigen::Vector3d> &start_end_derivative, double ratio, Eigen::MatrixXd &ctrl_pts, double &dt, double &time_inc)
  {
    // Deprecated:
    // Earlier Elastic versions tried to construct non-uniform knots in the
    // planner layer. The current system keeps planner output uniform and lets
    // traj_server perform execution-layer virtual time re-mapping.
    reparamBsplineGlobal(bspline, start_end_derivative, ratio, ctrl_pts, dt, time_inc);
  }
}
