#include "bspline_opt/uniform_bspline.h"
#include "nav_msgs/Odometry.h"
#include "traj_utils/Bspline.h"
#include "quadrotor_msgs/PositionCommand.h"
#include "std_msgs/Empty.h"
#include "visualization_msgs/Marker.h"
#include <ros/ros.h>
#include <fstream>
#include <cstdlib>
#include <algorithm>

ros::Publisher pos_cmd_pub;

quadrotor_msgs::PositionCommand cmd;
double pos_gain[3] = {0, 0, 0};
double vel_gain[3] = {0, 0, 0};

using ego_planner::UniformBspline;
using std::min;
using std::max;

bool receive_traj_ = false;
vector<UniformBspline> traj_;
double traj_duration_;
ros::Time start_time_;
int traj_id_;

// yaw control
double last_yaw_, last_yaw_dot_;
double time_forward_;

// ==========================================================
// 🚀 新增：ETWR 虚拟时间追踪与状态记录全局变量
// ==========================================================
std::string time_scaling_method_;
double max_vel_, max_acc_;
double s_virtual_ = 0.0;
int last_traj_id_ = -1;

bool is_flight_started_ = false;
ros::Time start_flight_time_;
// ==========================================================

void bsplineCallback(traj_utils::BsplineConstPtr msg)
{
  // parse pos traj

  Eigen::MatrixXd pos_pts(3, msg->pos_pts.size());

  Eigen::VectorXd knots(msg->knots.size());
  for (size_t i = 0; i < msg->knots.size(); ++i)
  {
    knots(i) = msg->knots[i];
  }

  for (size_t i = 0; i < msg->pos_pts.size(); ++i)
  {
    pos_pts(0, i) = msg->pos_pts[i].x;
    pos_pts(1, i) = msg->pos_pts[i].y;
    pos_pts(2, i) = msg->pos_pts[i].z;
  }

  UniformBspline pos_traj(pos_pts, msg->order, 0.1);
  pos_traj.setKnot(knots);

  // parse yaw traj

  // Eigen::MatrixXd yaw_pts(msg->yaw_pts.size(), 1);
  // for (int i = 0; i < msg->yaw_pts.size(); ++i) {
  //   yaw_pts(i, 0) = msg->yaw_pts[i];
  // }

  //UniformBspline yaw_traj(yaw_pts, msg->order, msg->yaw_dt);

  start_time_ = msg->start_time;
  traj_id_ = msg->traj_id;

  traj_.clear();
  traj_.push_back(pos_traj);
  traj_.push_back(traj_[0].getDerivative());
  traj_.push_back(traj_[1].getDerivative());

  traj_duration_ = traj_[0].getTimeSum();

  receive_traj_ = true;
}

std::pair<double, double> calculate_yaw(double progress_time, Eigen::Vector3d &pos, ros::Time &time_now, ros::Time &time_last)
{
  constexpr double PI = 3.1415926;
  constexpr double YAW_DOT_MAX_PER_SEC = PI;
  // constexpr double YAW_DOT_DOT_MAX_PER_SEC = PI;
  std::pair<double, double> yaw_yawdot(0, 0);
  double yaw = 0;
  double yawdot = 0;

  Eigen::Vector3d dir = progress_time + time_forward_ <= traj_duration_ ? traj_[0].evaluateDeBoorT(progress_time + time_forward_) - pos : traj_[0].evaluateDeBoorT(traj_duration_) - pos;
  double yaw_temp = dir.norm() > 0.1 ? atan2(dir(1), dir(0)) : last_yaw_;

  double max_yaw_change = YAW_DOT_MAX_PER_SEC * (time_now - time_last).toSec();
  if (yaw_temp - last_yaw_ > PI)
  {
    if (yaw_temp - last_yaw_ - 2 * PI < -max_yaw_change)
    {
      yaw = last_yaw_ - max_yaw_change;
      if (yaw < -PI)
        yaw += 2 * PI;

      yawdot = -YAW_DOT_MAX_PER_SEC;
    }
    else
    {
      yaw = yaw_temp;
      if (yaw - last_yaw_ > PI)
        yawdot = -YAW_DOT_MAX_PER_SEC;
      else
        yawdot = (yaw_temp - last_yaw_) / (time_now - time_last).toSec();
    }
  }
  else if (yaw_temp - last_yaw_ < -PI)
  {
    if (yaw_temp - last_yaw_ + 2 * PI > max_yaw_change)
    {
      yaw = last_yaw_ + max_yaw_change;
      if (yaw > PI)
        yaw -= 2 * PI;

      yawdot = YAW_DOT_MAX_PER_SEC;
    }
    else
    {
      yaw = yaw_temp;
      if (yaw - last_yaw_ < -PI)
        yawdot = YAW_DOT_MAX_PER_SEC;
      else
        yawdot = (yaw_temp - last_yaw_) / (time_now - time_last).toSec();
    }
  }
  else
  {
    if (yaw_temp - last_yaw_ < -max_yaw_change)
    {
      yaw = last_yaw_ - max_yaw_change;
      if (yaw < -PI)
        yaw += 2 * PI;

      yawdot = -YAW_DOT_MAX_PER_SEC;
    }
    else if (yaw_temp - last_yaw_ > max_yaw_change)
    {
      yaw = last_yaw_ + max_yaw_change;
      if (yaw > PI)
        yaw -= 2 * PI;

      yawdot = YAW_DOT_MAX_PER_SEC;
    }
    else
    {
      yaw = yaw_temp;
      if (yaw - last_yaw_ > PI)
        yawdot = -YAW_DOT_MAX_PER_SEC;
      else if (yaw - last_yaw_ < -PI)
        yawdot = YAW_DOT_MAX_PER_SEC;
      else
        yawdot = (yaw_temp - last_yaw_) / (time_now - time_last).toSec();
    }
  }

  if (fabs(yaw - last_yaw_) <= max_yaw_change)
    yaw = 0.5 * last_yaw_ + 0.5 * yaw; // nieve LPF
  yawdot = 0.5 * last_yaw_dot_ + 0.5 * yawdot;
  last_yaw_ = yaw;
  last_yaw_dot_ = yawdot;

  yaw_yawdot.first = yaw;
  yaw_yawdot.second = yawdot;

  return yaw_yawdot;
}

void cmdCallback(const ros::TimerEvent &e)
{
  if (!receive_traj_)
    return;

  ros::Time time_now = ros::Time::now();
  double t_cur = (time_now - start_time_).toSec();

  static ros::Time time_last = ros::Time::now();
  double dt = (time_now - time_last).toSec();
  // 🚨 安全断言：防止由于仿真系统卡顿导致的离谱 dt 毁掉积分
  if (dt <= 0.0 || dt > 0.05) dt = 0.01; 
  time_last = time_now;

  Eigen::Vector3d pos(Eigen::Vector3d::Zero()), vel(Eigen::Vector3d::Zero()), acc(Eigen::Vector3d::Zero()), pos_f;
  std::pair<double, double> yaw_yawdot(0, 0);

  double progress_time = t_cur; // 默认为真实时间 (Global 模式)
  double alpha = 1.0;           // 时间拉伸系数

  // =====================================================================
  // 🚀 核心机制：ETWR 虚拟时空映射追踪器 (Virtual Time Tracker)
  // =====================================================================
  if (time_scaling_method_ == "Elastic") 
  {
    // RHC 防截断：如果是新轨迹，虚拟时间重新从 0 开始积分
    if (traj_id_ != last_traj_id_) {
        s_virtual_ = 0.0;
        last_traj_id_ = traj_id_;
    }

    if (s_virtual_ < traj_duration_) 
    {
        // 🚨 防 Segfault 极限界：提取参考物理量时绝不能超界
        double safe_s = max(0.0, min(s_virtual_, traj_duration_ - 1e-4));
        
        Eigen::Vector3d v_spline = traj_[1].evaluateDeBoorT(safe_s);
        Eigen::Vector3d a_spline = traj_[2].evaluateDeBoorT(safe_s);

        // 计算物理极限下的时间扭曲系数
        double v_norm = v_spline.norm();
        double a_norm = a_spline.norm();
        
        double alpha_v = (v_norm > 0.01) ? (max_vel_ / v_norm) : 5.0;
        double alpha_a = (a_norm > 0.01) ? std::sqrt(max_acc_ / a_norm) : 5.0;

        // 取木桶效应的最短板，并预留 5% 的动力学安全余量
        alpha = min(alpha_v, alpha_a) * 0.95;
        
        // 钳制扭曲系数，防止极端奇异跳跃或卡死
        alpha = max(0.2, min(alpha, 3.0));

        // 积分推进虚拟进度条
        s_virtual_ += alpha * dt;
    }
    // 最终安全的虚拟时间
    progress_time = max(0.0, min(s_virtual_, traj_duration_)); 
  }

  // =====================================================================
  // 轨迹评估与下发
  // =====================================================================
  if (progress_time < traj_duration_ && progress_time >= 0.0)
  {
    pos = traj_[0].evaluateDeBoorT(progress_time);

    if (time_scaling_method_ == "Global") {
        vel = traj_[1].evaluateDeBoorT(progress_time);
        acc = traj_[2].evaluateDeBoorT(progress_time);
    } else {
        // 🚀 链式求导：当时间进度按 alpha 缩放时，速度乘以 alpha，加速度乘以 alpha^2
        double safe_p = max(0.0, min(progress_time, traj_duration_ - 1e-4));
        vel = traj_[1].evaluateDeBoorT(safe_p) * alpha;
        acc = traj_[2].evaluateDeBoorT(safe_p) * (alpha * alpha);
    }

    /*** calculate yaw ***/
    yaw_yawdot = calculate_yaw(progress_time, pos, time_now, time_last);

    double tf = min(traj_duration_, progress_time + 2.0);
    pos_f = traj_[0].evaluateDeBoorT(tf);
  }
  else if (progress_time >= traj_duration_)
  {
    /* hover when finish traj_ */
    pos = traj_[0].evaluateDeBoorT(traj_duration_);
    vel.setZero();
    acc.setZero();

    yaw_yawdot.first = last_yaw_;
    yaw_yawdot.second = 0;

    pos_f = pos;
  }
  else
  {
    cout << "[Traj server]: invalid time." << endl;
  }

  // =====================================================================
  // 📊 真机遥测数据埋点记录 (Real-world Telemetry)
  // =====================================================================
  if (!is_flight_started_) {
      start_flight_time_ = ros::Time::now();
      is_flight_started_ = true;
  }

  double t_current = (ros::Time::now() - start_flight_time_).toSec();
  
  // 仅在轨迹未结束，且每 0.1 秒记录一次
  static double last_record_time = 0;
  if (progress_time < traj_duration_ && t_current - last_record_time > 0.1) {
      std::string home_dir = getenv("HOME");
      std::ofstream real_profile(home_dir + "/.ros/ego_real_flight_log.csv", std::ios::out | std::ios::app);
      if (real_profile.is_open()) {
          if (real_profile.tellp() == 0) {
              real_profile << "Real_Time,Velocity,Acceleration,Method\n";
          }
          real_profile << t_current << "," 
                       << vel.norm() << "," 
                       << acc.norm() << "," 
                       << time_scaling_method_ << "\n";
          real_profile.close();
      }
      last_record_time = t_current;
  }

  // ---------------- 发送指令 ----------------
  cmd.header.stamp = time_now;
  cmd.header.frame_id = "world";
  cmd.trajectory_flag = quadrotor_msgs::PositionCommand::TRAJECTORY_STATUS_READY;
  cmd.trajectory_id = traj_id_;

  cmd.position.x = pos(0);
  cmd.position.y = pos(1);
  cmd.position.z = pos(2);

  cmd.velocity.x = vel(0);
  cmd.velocity.y = vel(1);
  cmd.velocity.z = vel(2);

  cmd.acceleration.x = acc(0);
  cmd.acceleration.y = acc(1);
  cmd.acceleration.z = acc(2);

  cmd.yaw = yaw_yawdot.first;
  cmd.yaw_dot = yaw_yawdot.second;

  last_yaw_ = cmd.yaw;

  pos_cmd_pub.publish(cmd);
}

int main(int argc, char **argv)
{
  ros::init(argc, argv, "traj_server");
  // ros::NodeHandle node;
  ros::NodeHandle nh("~");
  nh.param("manager/max_vel", max_vel_, 2.0);
  nh.param("manager/max_acc", max_acc_, 5.0);
  nh.param("manager/time_scaling_method", time_scaling_method_, std::string("Global"));

  ros::Subscriber bspline_sub = nh.subscribe("planning/bspline", 10, bsplineCallback);

  pos_cmd_pub = nh.advertise<quadrotor_msgs::PositionCommand>("/position_cmd", 50);

  ros::Timer cmd_timer = nh.createTimer(ros::Duration(0.01), cmdCallback);
  
  /* control parameter */
  cmd.kx[0] = pos_gain[0];
  cmd.kx[1] = pos_gain[1];
  cmd.kx[2] = pos_gain[2];

  cmd.kv[0] = vel_gain[0];
  cmd.kv[1] = vel_gain[1];
  cmd.kv[2] = vel_gain[2];

  nh.param("traj_server/time_forward", time_forward_, -1.0);
  last_yaw_ = 0.0;
  last_yaw_dot_ = 0.0;

  ros::Duration(1.0).sleep();

  ROS_WARN("[Traj server]: ready.");

  ros::spin();

  return 0;
}