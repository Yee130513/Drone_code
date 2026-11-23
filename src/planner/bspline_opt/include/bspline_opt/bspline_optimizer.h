#ifndef _BSPLINE_OPTIMIZER_H_
#define _BSPLINE_OPTIMIZER_H_

#include <Eigen/Eigen>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
//#include <path_searching/jps_theta_planner.h>
#include <path_searching/dyn_a_star.h>
//#include <path_searching/theta_star.h> 
#include <bspline_opt/uniform_bspline.h>
#include <plan_env/grid_map.h>
#include <plan_env/obj_predictor.h>
#include <ros/ros.h>
#include "bspline_opt/lbfgs.hpp"
#include <traj_utils/plan_container.hpp>

// Gradient and elasitc band optimization

// Input: a signed distance field and a sequence of points
// Output: the optimized sequence of points
// The format of points: N x 3 matrix, each row is a point
namespace ego_planner
{

  class ControlPoints
  {
  public:
    double clearance;
    int size;
    Eigen::MatrixXd points;
    std::vector<std::vector<Eigen::Vector3d>> base_point; // The point at the statrt of the direction vector (collision point)
    std::vector<std::vector<Eigen::Vector3d>> direction;  // Direction vector, must be normalized.
    std::vector<bool> flag_temp;                          // A flag that used in many places. Initialize it everytime before using it.
    // std::vector<bool> occupancy;
    std::vector<std::vector<double>> local_clearance; // per control point soft collision threshold, fallback to clearance when invalid

    void resize(const int size_set)
    {
      size = size_set;

      base_point.clear();
      direction.clear();
      local_clearance.clear(); // 馃殌 鏂板
      flag_temp.clear();
      // occupancy.clear();

      points.resize(3, size_set);
      base_point.resize(size);
      direction.resize(size);
      local_clearance.resize(size); // 馃殌 鏂板
      flag_temp.resize(size);
      // occupancy.resize(size);
    }

    void segment(ControlPoints &buf, const int start, const int end)
    {

      if (start < 0 || end >= size || points.rows() != 3)
      {
        ROS_ERROR("Wrong segment index! start=%d, end=%d", start, end);
        return;
      }

      buf.resize(end - start + 1);
      buf.points = points.block(0, start, 3, end - start + 1);
      buf.clearance = clearance;
      buf.size = end - start + 1;
      for (int i = start; i <= end; i++)
      {
        buf.base_point[i - start] = base_point[i];
        buf.direction[i - start] = direction[i];
         buf.local_clearance[i - start] = local_clearance[i]; // 馃殌 鏂板鏁版嵁鍚屾
        // if ( buf.base_point[i - start].size() > 1 )
        // {
        //   ROS_ERROR("buf.base_point[i - start].size()=%d, base_point[i].size()=%d", buf.base_point[i - start].size(), base_point[i].size());
        // }
      }

      // cout << "RichInfoOneSeg_temp, insede" << endl;
      // for ( int k=0; k<buf.size; k++ )
      //   if ( buf.base_point[k].size() > 0 )
      //   {
      //     cout << "###" << buf.points.col(k).transpose() << endl;
      //     for (int k2 = 0; k2 < buf.base_point[k].size(); k2++)
      //     {
      //       cout << "      " << buf.base_point[k][k2].transpose() << " @ " << buf.direction[k][k2].transpose() << endl;
      //     }
      //   }
    }
  };

  class BsplineOptimizer
  {

  public:
    BsplineOptimizer() {}
    ~BsplineOptimizer() {}

    /* main API */
    void setEnvironment(const GridMap::Ptr &map);
    void setEnvironment(const GridMap::Ptr &map, const fast_planner::ObjPredictor::Ptr mov_obj);
    void setParam(ros::NodeHandle &nh);
    Eigen::MatrixXd BsplineOptimizeTraj(const Eigen::MatrixXd &points, const double &ts,
                                        const int &cost_function, int max_num_id, int max_time_id);

    /* helper function */

    // required inputs
    void setControlPoints(const Eigen::MatrixXd &points);
    void setBsplineInterval(const double &ts);
    void setSwarmTrajs(SwarmTrajData *swarm_trajs_ptr);
    void setDroneId(const int drone_id);

    // optional inputs
    void setGuidePath(const vector<Eigen::Vector3d> &guide_pt);
    void setWaypoints(const vector<Eigen::Vector3d> &waypts,
                      const vector<int> &waypt_idx); // N-2 constraints at most
    void setLocalTargetPt(const Eigen::Vector3d local_target_pt) { local_target_pt_ = local_target_pt; };

    void optimize();

    ControlPoints getControlPoints() { return cps_; };
    //ThetaStar::Ptr a_star_;//=========================================
    AStar::Ptr a_star_;//=======================================
    //JPSThetaPlanner::Ptr a_star_;
    std::vector<Eigen::Vector3d> ref_pts_;

    std::vector<ControlPoints> distinctiveTrajs(vector<std::pair<int, int>> segments);
    std::vector<std::pair<int, int>> initControlPoints(Eigen::MatrixXd &init_points, bool flag_first_init = true);
    bool BsplineOptimizeTrajRebound(Eigen::MatrixXd &optimal_points, double ts); // must be called after initControlPoints()
    bool BsplineOptimizeTrajRebound(Eigen::MatrixXd &optimal_points, double &final_cost, const ControlPoints &control_points, double ts);
    bool BsplineOptimizeTrajRefine(const Eigen::MatrixXd &init_points, const double ts, Eigen::MatrixXd &optimal_points);

    inline int getOrder(void) { return order_; }
    inline double getSwarmClearance(void) { return swarm_clearance_; }
    inline int getIterNum(void) const { return iter_num_; }
    inline int getInitCollisionPts() const { return init_collision_pts_; }
    inline bool getLastFrontendFailed() const { return last_frontend_failed_; }
    inline std::string getLastFrontendFailReason() const { return last_frontend_fail_reason_; }
    inline int getLastGuidePathCount() const { return last_guide_path_count_; }
    inline int getLastGuideInvalidCount() const { return last_guide_invalid_count_; }
    inline bool isDynamicClearanceEnabled() const { return enable_dynamic_clearance_; }
    inline int getDynamicClearanceCount() const { return dynamic_clearance_count_; }
    inline int getDynamicClearanceFallbackCount() const { return dynamic_clearance_fallback_count_; }
    inline double getDynamicClearanceMin() const
    {
      return dynamic_clearance_count_ > 0 ? dynamic_clearance_min_ : cps_.clearance;
    }
    inline double getDynamicClearanceMean() const
    {
      return dynamic_clearance_count_ > 0 ? dynamic_clearance_sum_ / dynamic_clearance_count_ : cps_.clearance;
    }
    inline double getDynamicClearanceMax() const
    {
      return dynamic_clearance_count_ > 0 ? dynamic_clearance_max_ : cps_.clearance;
    }

  private:
    int init_collision_pts_{0};

    bool last_frontend_failed_{false};
    std::string last_frontend_fail_reason_{"NONE"};
    int last_guide_path_count_{0};
    int last_guide_invalid_count_{0};

    bool enable_guide_path_validation_{true};
    bool enable_dynamic_clearance_{false};

    int guide_min_points_{3};
    double guide_endpoint_tolerance_{0.30};
    double guide_validate_resolution_ratio_{0.50};

    // dynamic safe threshold parameters
    // Only affects rebound soft collision cost. It must not change grid_map hard collision checks.
    double dynamic_clearance_gain_{0.25};
    double dynamic_clearance_beta_{0.40};
    double dynamic_clearance_filter_alpha_{0.70};
    double dynamic_clearance_upper_ratio_{1.20};
    double dynamic_clearance_lower_ratio_{1.00};
    double dynamic_clearance_global_max_{0.65};
    double dynamic_clearance_min_guide_dist_{0.02};
    bool dynamic_clearance_debug_{false};

    // runtime state for smoothing and logging
    std::vector<double> dynamic_clearance_prev_;
    std::vector<bool> dynamic_clearance_has_prev_;
    int dynamic_clearance_count_{0};
    int dynamic_clearance_fallback_count_{0};
    double dynamic_clearance_min_{0.0};
    double dynamic_clearance_max_{0.0};
    double dynamic_clearance_sum_{0.0};

    GridMap::Ptr grid_map_;
    fast_planner::ObjPredictor::Ptr moving_objs_;
    SwarmTrajData *swarm_trajs_{NULL}; // Can not use shared_ptr and no need to free
    int drone_id_;

    enum FORCE_STOP_OPTIMIZE_TYPE
    {
      DONT_STOP,
      STOP_FOR_REBOUND,
      STOP_FOR_ERROR
    } force_stop_type_;

    // main input
    // Eigen::MatrixXd control_points_;     // B-spline control points, N x dim
    double bspline_interval_; // B-spline knot span
    Eigen::Vector3d end_pt_;  // end of the trajectory
    // int             dim_;                // dimension of the B-spline
    //
    vector<Eigen::Vector3d> guide_pts_; // geometric guiding path points, N-6
    vector<Eigen::Vector3d> waypoints_; // waypts constraints
    vector<int> waypt_idx_;             // waypts constraints index
                                        //
    int max_num_id_, max_time_id_;      // stopping criteria
    int cost_function_;                 // used to determine objective function
    double start_time_;                 // global time for moving obstacles

    /* optimization parameters */
    int order_;                    // bspline degree
    double lambda1_;               // jerk smoothness weight
    double lambda2_, new_lambda2_; // distance weight
    double lambda3_;               // feasibility weight
    double lambda4_;               // curve fitting

    int a;
    //
    double dist0_, swarm_clearance_; // safe distance
    double max_vel_, max_acc_;       // dynamic limits

    int variable_num_;              // optimization variables
    int iter_num_;                  // iteration of the solver
    Eigen::VectorXd best_variable_; //
    double min_cost_;               //

    Eigen::Vector3d local_target_pt_; 

#define INIT_min_ellip_dist_ 123456789.0123456789
    double min_ellip_dist_;

    ControlPoints cps_;

    /* cost function */
    /* calculate each part of cost function with control points q as input */

    static double costFunction(const std::vector<double> &x, std::vector<double> &grad, void *func_data);
    void combineCost(const std::vector<double> &x, vector<double> &grad, double &cost);

    // q contains all control points
    void calcSmoothnessCost(const Eigen::MatrixXd &q, double &cost, Eigen::MatrixXd &gradient, bool falg_use_jerk = true);
    void calcFeasibilityCost(const Eigen::MatrixXd &q, double &cost, Eigen::MatrixXd &gradient);
    void calcTerminalCost(const Eigen::MatrixXd &q, double &cost, Eigen::MatrixXd &gradient);
    void calcDistanceCostRebound(const Eigen::MatrixXd &q, double &cost, Eigen::MatrixXd &gradient, int iter_num, double smoothness_cost);
    void calcMovingObjCost(const Eigen::MatrixXd &q, double &cost, Eigen::MatrixXd &gradient);
    void calcSwarmCost(const Eigen::MatrixXd &q, double &cost, Eigen::MatrixXd &gradient);
    void calcFitnessCost(const Eigen::MatrixXd &q, double &cost, Eigen::MatrixXd &gradient);
    bool check_collision_and_rebound(void);
    bool buildGuidePathForSegment(const Eigen::Vector3d& in,
                                  const Eigen::Vector3d& out,
                                  std::vector<Eigen::Vector3d>& guide_path,
                                  std::string& fail_reason);
    bool validateGuidePath(const std::vector<Eigen::Vector3d>& guide_path,
                           const Eigen::Vector3d& in,
                           const Eigen::Vector3d& out,
                           std::string& fail_reason) const;
    bool segmentFreeSampledForGuide(const Eigen::Vector3d& p0,
                                    const Eigen::Vector3d& p1,
                                    bool ignore_start_occupancy = false,
                                    bool ignore_end_occupancy = false) const;
    bool isFiniteVec(const Eigen::Vector3d& p) const;
    bool isValidDirectionAndBase(const Eigen::Vector3d& base_pt,
                                 const Eigen::Vector3d& direction) const;
    inline double clampDouble(double x, double lo, double hi) const
    {
      return std::max(lo, std::min(hi, x));
    }
    void resetDynamicClearanceRuntime(int ctrl_pt_num);
    double sanitizeLocalClearance(double sf_local) const;
    void recordDynamicClearanceStat(double sf_local, bool fallback);
    double estimateGuidePointClearance(const Eigen::Vector3d& guide_pt,
                                       const Eigen::Vector3d& push_dir) const;
    double computeLocalClearanceThreshold(double guide_deviation,
                                          double guide_clearance,
                                          int ctrl_pt_id);
    double computeLocalClearanceThreshold(double guide_deviation, int ctrl_pt_id);
    bool matchGuidePointForControlPoint(const std::vector<Eigen::Vector3d>& guide_path,
                                        const Eigen::Vector3d& ctrl_pt,
                                        const Eigen::Vector3d& ctrl_pts_law,
                                        double progress_ratio,
                                        int& last_matched_idx,
                                        Eigen::Vector3d& matched_guide_pt) const;

    static int earlyExit(void *func_data, const double *x, const double *g, const double fx, const double xnorm, const double gnorm, const double step, int n, int k, int ls);
    static double costFunctionRebound(void *func_data, const double *x, double *grad, const int n);
    static double costFunctionRefine(void *func_data, const double *x, double *grad, const int n);

    bool rebound_optimize(double &final_cost);
    bool refine_optimize();
    void combineCostRebound(const double *x, double *grad, double &f_combine, const int n);
    void combineCostRefine(const double *x, double *grad, double &f_combine, const int n);

    /* for benckmark evaluation only */
  public:
    typedef unique_ptr<BsplineOptimizer> Ptr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  };

} // namespace ego_planner
#endif
