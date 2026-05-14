// #include <fstream>
#include <plan_manage/planner_manager.h>
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
    nh.param("manager/time_scaling_method", time_scaling_method_, std::string("Elastic"));
    /* read algorithm parameters，参数读取 */
    nh.param("manager/max_vel", pp_.max_vel_, -1.0);
    // ... 原有参数读取 ...
    nh.param("manager/drone_id", pp_.drone_id, -1);

    // [新增] 读取消融实验的算法名称，默认为 "A_star"
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

    vector<std::pair<int, int>> segments;
    segments = bspline_optimizer_->initControlPoints(ctrl_pts, true);

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

    /* ================= ABLATION STUDY: DATA LOGGER ================= */
   auto log_ablation_data = [&](bool is_success, double jerk_cost) {
        std::string home_dir = getenv("HOME");
        std::string file_name = home_dir + "/.ros/sca.csv";
        
        std::ofstream file(file_name, std::ios::out | std::ios::app);
        if (file.is_open()) {
            file.seekp(0, std::ios::end);
            // 🚀 更新表头：加入了前端与耦合指标
            if (file.tellp() == 0) {
                file << "Timestamp,Planner_Type,Search_Time_ms,Opt_Time_ms,Opt_Iterations,Nodes_Expanded,Path_Length_m,Init_Collision_Pts,Is_Opt_Success,Traj_Jerk_Cost\n";
            }
            
            // 安全获取前端数据 (防止 a_star_ 指针为空)
            int nodes = -1;
            double path_len = -1.0;
            int init_cols = -1;
            
            if (bspline_optimizer_ != nullptr) {
                init_cols = bspline_optimizer_->getInitCollisionPts();
                if (bspline_optimizer_->a_star_ != nullptr) {
                    nodes = bspline_optimizer_->a_star_->getExpandedNodes();
                    path_len = bspline_optimizer_->a_star_->getPathLength();
                }
            }
            // 写入所有核心评价指标
            file << ros::Time::now().toSec() << ","
                 << planner_type_ << ","
                 << t_init.toSec() * 1000.0 << ","       
                 << t_opt.toSec() * 1000.0 << ","        
                 << bspline_optimizer_->getIterNum() << "," 
                 << nodes << ","                         // 📈 新增：扩展节点数
                 << path_len << ","                      // 📈 新增：拓扑路径长度
                 << bspline_optimizer_->getInitCollisionPts() << "," // 📈 新增：碰撞控制点数
                 << (is_success ? 1 : 0) << ","          
                 << jerk_cost << "\n";                   
            file.close();
        }
    };
    /* ============================================================== */

    cout << "plan_success=" << flag_step_1_success << endl;
    if (!flag_step_1_success)
    {
      visualization_->displayOptimalList(ctrl_pts, 0);
      continous_failures_count_++;
      log_ablation_data(false, -1.0); // 埋点：优化失败，Jerk 无效记为 -1.0
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
        log_ablation_data(false, -1.0); // 埋点：Refine 阶段撞墙失败
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
    // 🚀 ABLATION STUDY LOGGER: 动力学与时间重分配实验埋点
    // =====================================================================
    if (flag_step_1_success &&  (flag_step_2_success)) {
        UniformBspline vel_traj = pos.getDerivative();
        UniformBspline acc_traj = vel_traj.getDerivative();
        
        double duration = pos.getTimeSum();
        double max_v = 0.0, max_a = 0.0, sum_v = 0.0;
        int sample_count = 0;

        // 1. 采样提取微观动力学曲线数据 (每 0.02 秒采一个点)
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

            // 将局部速度曲线写入文件
            if (write_profile) {
                profile_file << local_data_.traj_id_ << "," << t << "," 
                             << v_norm << "," << a_norm << "," << time_scaling_method_ << "\n";
            }
        }
        if (write_profile) profile_file.close();

        double avg_v = sum_v / sample_count;

        // 2. 写入宏观统计数据
        std::ofstream stats_file(home_dir + "/.ros/ego_macro_stats.csv", std::ios::out | std::ios::app);
        if (stats_file.is_open()) {
            if (stats_file.tellp() == 0) {
                stats_file << "Timestamp,Method,Duration,Avg_Vel,Max_Vel,Max_Acc,Refine_Time_ms\n";
            }
            stats_file << ros::Time::now().toSec() << ","
                       << time_scaling_method_ << ","
                       << duration << ","
                       << avg_v << ","
                       << max_v << ","
                       << max_a << ","
                       << t_refine.toSec() * 1000.0 << "\n";
            stats_file.close();
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

    // std::cout << "ratio: " << ratio << std::endl;
    //reparamBspline(traj, start_end_derivative, ratio, ctrl_pts, ts, t_inc);
    if (time_scaling_method_ == "Global") {
    reparamBsplineGlobal(traj, start_end_derivative, ratio, ctrl_pts, ts, t_inc);
} else {
    reparamBsplineElastic(traj, start_end_derivative, ratio, ctrl_pts, ts, t_inc);
}

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

// ---------------- 修改位置: planner_manager.cpp ----------------
// 用以下代码完全替换原有的 reparamBspline 函数

// ---------------- 修改位置: planner_manager.cpp ----------------
// 彻底重构的 True ETWR 算法：纯时域非均匀时间重分配
void EGOPlannerManager::reparamBsplineElastic(UniformBspline &bspline, vector<Eigen::Vector3d> &start_end_derivative, double ratio, Eigen::MatrixXd &ctrl_pts, double &dt, double &time_inc)
{
    // 获取优化后的原始控制点（绝对不改变其空间几何位置）
    ctrl_pts = bspline.getControlPoint();
    int p = 3;
    int n = ctrl_pts.cols() - 1; 
    int m = n + p + 1; 

    // 1. 基于弦长参数化 (Chord-Length) 的初始局部时间分配
    // 强制目标速度逼近 pp_.max_vel_
    std::vector<double> dt_intervals(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double dist = (ctrl_pts.col(i+1) - ctrl_pts.col(i)).norm();
        dt_intervals[i] = dist / pp_.max_vel_;
        // 设定最低安全时间间隔
        dt_intervals[i] = std::max(dt_intervals[i], 0.02);
    }

    // 2. 局部时间膨胀 (Local Time Dilation) —— 解决加速度超标问题
    // 在急转弯处，速度方向变化剧烈，需要拉长局部时间来降低向心加速度
    for (int pass = 0; pass < 2; ++pass) { // 两次滤波迭代，保证动力学平滑过渡
        for (int i = 0; i < n - 1; ++i) {
            Eigen::Vector3d v1 = (ctrl_pts.col(i+1) - ctrl_pts.col(i)) / dt_intervals[i];
            Eigen::Vector3d v2 = (ctrl_pts.col(i+2) - ctrl_pts.col(i+1)) / dt_intervals[i+1];
            
            // 近似向心加速度计算
            double acc = (v2 - v1).norm() / ((dt_intervals[i] + dt_intervals[i+1]) / 2.0);

            if (acc > pp_.max_acc_) {
                // 计算需要膨胀的时间比率
                double dilation_ratio = sqrt(acc / pp_.max_acc_);
                dt_intervals[i] *= dilation_ratio;
                dt_intervals[i+1] *= dilation_ratio;
            }
        }
    }

    // 3. 构建非均匀节点向量 (Non-Uniform Knot Vector)
    Eigen::VectorXd knots = Eigen::VectorXd::Zero(m + 1);
    
    // 起点约束 (Clamped Start)
    for (int i = 0; i <= p; ++i) knots(i) = 0.0;

    // 内部非均匀时间映射
    for (int i = p + 1; i <= n; ++i) {
        knots(i) = knots(i-1) + dt_intervals[i - p];
    }

    // 终点约束 (Clamped End)
    for (int i = n + 1; i <= m; ++i) {
        knots(i) = knots(i-1) + dt_intervals.back();
    }

    // 4. 将非均匀时间属性赋予轨迹，直接推给底层飞控！
    bspline.setKnot(knots);
    
    // 更新状态机所需的时间增量
    time_inc = knots(m - p) - knots(p);
    
    // 这里的 dt 仅作为格式兼容保留，因为底层 traj_server 实际上是依赖 knots 解析的
    dt = dt_intervals[0]; 
}
}
