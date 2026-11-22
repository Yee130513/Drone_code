#ifndef _PLANNER_MANAGER_H_
#define _PLANNER_MANAGER_H_

#include <fstream>
#include <cstdlib>
#include <string>

#include <stdlib.h>

#include <bspline_opt/bspline_optimizer.h>
#include <bspline_opt/uniform_bspline.h>
#include <traj_utils/DataDisp.h>
#include <plan_env/grid_map.h>
#include <plan_env/obj_predictor.h>
#include <traj_utils/plan_container.hpp>
#include <ros/ros.h>
#include <traj_utils/planning_visualization.h>

namespace ego_planner
{

  // Fast Planner Manager
  // Key algorithms of mapping and planning are called

  class EGOPlannerManager
  {
    // SECTION stable
  public:
    EGOPlannerManager();
    ~EGOPlannerManager();

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /* main planning interface */
    bool reboundReplan(Eigen::Vector3d start_pt, Eigen::Vector3d start_vel, Eigen::Vector3d start_acc,
                       Eigen::Vector3d end_pt, Eigen::Vector3d end_vel, bool flag_polyInit, bool flag_randomPolyTraj);
    bool EmergencyStop(Eigen::Vector3d stop_pos);
    bool planGlobalTraj(const Eigen::Vector3d &start_pos, const Eigen::Vector3d &start_vel, const Eigen::Vector3d &start_acc,
                        const Eigen::Vector3d &end_pos, const Eigen::Vector3d &end_vel, const Eigen::Vector3d &end_acc);
    bool planGlobalTrajWaypoints(const Eigen::Vector3d &start_pos, const Eigen::Vector3d &start_vel, const Eigen::Vector3d &start_acc,
                                 const std::vector<Eigen::Vector3d> &waypoints, const Eigen::Vector3d &end_vel, const Eigen::Vector3d &end_acc);

    void initPlanModules(ros::NodeHandle &nh, PlanningVisualization::Ptr vis = NULL);

    void deliverTrajToOptimizer(void) { bspline_optimizer_->setSwarmTrajs(&swarm_trajs_buf_); };

    void setDroneIdtoOpt(void) { bspline_optimizer_->setDroneId(pp_.drone_id); }

    double getSwarmClearance(void) { return bspline_optimizer_->getSwarmClearance(); }

    bool checkCollision(int drone_id);
    inline int getContinuousFailuresCount() const { return continous_failures_count_; }
    inline void resetContinuousFailuresCount() { continous_failures_count_ = 0; }
    

    PlanParameters pp_;
    LocalTrajData local_data_;
    GlobalTrajData global_data_;
    GridMap::Ptr grid_map_;
    fast_planner::ObjPredictor::Ptr obj_predictor_;    
    SwarmTrajData swarm_trajs_buf_;

  private:
    /* main planning algorithms & modules */
    PlanningVisualization::Ptr visualization_;

    // ros::Publisher obj_pub_; //zx-todo 

    BsplineOptimizer::Ptr bspline_optimizer_;

    int continous_failures_count_{0};

    void updateTrajInfo(const UniformBspline &position_traj, const ros::Time time_now);

    void reparamBspline(UniformBspline &bspline, vector<Eigen::Vector3d> &start_end_derivative, double ratio, Eigen::MatrixXd &ctrl_pts, double &dt,
                        double &time_inc);

    bool refineTrajAlgo(UniformBspline &traj, vector<Eigen::Vector3d> &start_end_derivative, double ratio, double &ts, Eigen::MatrixXd &optimal_control_points);

    std::string planner_type_;
    
    // !SECTION stable

    // SECTION developing
    // Deprecated name. This now only represents conservative planner-layer
    // time allocation for the uniform fallback trajectory.
    std::string time_scaling_method_;

// 声明两种不同的重采样函数
void reparamBsplineGlobal(UniformBspline &bspline, vector<Eigen::Vector3d> &start_end_derivative, double ratio, Eigen::MatrixXd &ctrl_pts, double &dt, double &time_inc);
// Deprecated wrapper: current implementation directly calls reparamBsplineGlobal.
// Non-uniform knots are no longer generated in the planner layer.
void reparamBsplineElastic(UniformBspline &bspline, vector<Eigen::Vector3d> &start_end_derivative, double ratio, Eigen::MatrixXd &ctrl_pts, double &dt, double &time_inc);

  public:
    typedef unique_ptr<EGOPlannerManager> Ptr;

    // !SECTION
  };
} // namespace ego_planner

#endif
