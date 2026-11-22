#include <algorithm>
#include <string>
#include <vector>

#include <geometry_msgs/Point.h>
#include <geometry_msgs/PoseStamped.h>
#include <ros/ros.h>
#include <visualization_msgs/Marker.h>

class StartGoalMarkerNode
{
public:
  StartGoalMarkerNode() : nh_("~")
  {
    nh_.param<std::string>("frame_id", frame_id_, "world");
    nh_.param("drone_id", drone_id_, 0);
    nh_.param("flight_type", flight_type_, 1);
    nh_.param("init_x", init_x_, 0.0);
    nh_.param("init_y", init_y_, 0.0);
    nh_.param("init_z", init_z_, 0.0);
    nh_.param("point_num", point_num_, 0);
    nh_.param("marker_rate", marker_rate_, 2.0);

    start_.x = init_x_;
    start_.y = init_y_;
    start_.z = init_z_;

    readWaypoints();

    const std::string prefix = "/drone_" + std::to_string(drone_id_) + "_debug/";
    start_pub_ = nh_global_.advertise<visualization_msgs::Marker>(prefix + "start_marker", 10, true);
    goal_pub_ = nh_global_.advertise<visualization_msgs::Marker>(prefix + "local_goal_marker", 10, true);
    goal_sub_ = nh_global_.subscribe("/move_base_simple/goal", 10, &StartGoalMarkerNode::goalCallback, this);

    publishStartMarkers();
    clearGoalMarkers();
  }

  void spin()
  {
    ros::Rate rate(std::max(0.1, marker_rate_));
    while (ros::ok())
    {
      publishStartMarkers();

      if (flight_type_ == 1)
      {
        if (has_goal_)
        {
          publishFlightTypeOneGoal();
        }
      }
      else if (flight_type_ == 2)
      {
        clearGoalMarkers();
        publishWaypoints();
      }

      ros::spinOnce();
      rate.sleep();
    }
  }

private:
  visualization_msgs::Marker baseMarker(const std::string& ns, const int id, const int type) const
  {
    visualization_msgs::Marker marker;
    marker.header.frame_id = frame_id_;
    marker.header.stamp = ros::Time::now();
    marker.ns = ns;
    marker.id = id;
    marker.type = type;
    marker.action = visualization_msgs::Marker::ADD;
    marker.pose.orientation.w = 1.0;
    marker.lifetime = ros::Duration(0.0);
    return marker;
  }

  void setColor(visualization_msgs::Marker& marker, const double r, const double g,
                const double b, const double a) const
  {
    marker.color.r = r;
    marker.color.g = g;
    marker.color.b = b;
    marker.color.a = a;
  }

  void goalCallback(const geometry_msgs::PoseStampedConstPtr& msg)
  {
    if (flight_type_ != 1)
    {
      return;
    }

    clearGoalMarkers();
    goal_ = msg->pose.position;
    has_goal_ = true;
    publishFlightTypeOneGoal();
  }

  void readWaypoints()
  {
    waypoints_.clear();
    const int valid_num = std::max(0, std::min(point_num_, 5));

    for (int i = 0; i < valid_num; ++i)
    {
      geometry_msgs::Point point;
      const std::string key = "point" + std::to_string(i);
      nh_.param(key + "_x", point.x, 0.0);
      nh_.param(key + "_y", point.y, 0.0);
      nh_.param(key + "_z", point.z, 0.0);
      waypoints_.push_back(point);
    }
  }

  void clearGoalMarkers()
  {
    visualization_msgs::Marker clear;
    clear.header.frame_id = frame_id_;
    clear.header.stamp = ros::Time::now();
    clear.action = visualization_msgs::Marker::DELETEALL;
    goal_pub_.publish(clear);
  }

  void publishStartMarkers()
  {
    visualization_msgs::Marker sphere = baseMarker("start", 0, visualization_msgs::Marker::SPHERE);
    sphere.pose.position = start_;
    sphere.scale.x = 0.35;
    sphere.scale.y = 0.35;
    sphere.scale.z = 0.35;
    setColor(sphere, 0.05, 0.55, 0.15, 1.0);
    start_pub_.publish(sphere);

    visualization_msgs::Marker text = baseMarker("start_text", 1, visualization_msgs::Marker::TEXT_VIEW_FACING);
    text.pose.position = start_;
    text.pose.position.z += 0.45;
    text.scale.z = 0.35;
    text.text = "START";
    setColor(text, 0.05, 0.05, 0.05, 1.0);
    start_pub_.publish(text);
  }

  void publishFlightTypeOneGoal()
  {
    visualization_msgs::Marker sphere = baseMarker("local_goal", 0, visualization_msgs::Marker::SPHERE);
    sphere.pose.position = goal_;
    sphere.scale.x = 0.45;
    sphere.scale.y = 0.45;
    sphere.scale.z = 0.45;
    setColor(sphere, 1.0, 0.75, 0.05, 1.0);
    goal_pub_.publish(sphere);

    visualization_msgs::Marker text = baseMarker("local_goal_text", 1, visualization_msgs::Marker::TEXT_VIEW_FACING);
    text.pose.position = goal_;
    text.pose.position.z += 0.50;
    text.scale.z = 0.35;
    text.text = "GOAL";
    setColor(text, 0.05, 0.05, 0.05, 1.0);
    goal_pub_.publish(text);

  }

  void publishWaypoints()
  {
    for (std::size_t i = 0; i < waypoints_.size(); ++i)
    {
      visualization_msgs::Marker sphere =
          baseMarker("waypoint", 100 + static_cast<int>(i), visualization_msgs::Marker::SPHERE);
      sphere.pose.position = waypoints_[i];
      sphere.scale.x = 0.45;
      sphere.scale.y = 0.45;
      sphere.scale.z = 0.45;
      if (i == 0)
      {
        setColor(sphere, 1.0, 0.75, 0.05, 1.0);
      }
      else
      {
        setColor(sphere, 0.0, 0.65, 0.9, 1.0);
      }
      goal_pub_.publish(sphere);

      visualization_msgs::Marker text =
          baseMarker("waypoint_text", 200 + static_cast<int>(i), visualization_msgs::Marker::TEXT_VIEW_FACING);
      text.pose.position = waypoints_[i];
      text.pose.position.z += 0.50;
      text.scale.z = 0.32;
      text.text = "WP" + std::to_string(i);
      setColor(text, 0.05, 0.05, 0.05, 1.0);
      goal_pub_.publish(text);
    }

  }

  ros::NodeHandle nh_;
  ros::NodeHandle nh_global_;
  ros::Publisher start_pub_;
  ros::Publisher goal_pub_;
  ros::Subscriber goal_sub_;

  std::string frame_id_;
  int drone_id_ = 0;
  int flight_type_ = 1;
  int point_num_ = 0;
  double init_x_ = 0.0;
  double init_y_ = 0.0;
  double init_z_ = 0.0;
  double marker_rate_ = 2.0;
  bool has_goal_ = false;

  geometry_msgs::Point start_;
  geometry_msgs::Point goal_;
  std::vector<geometry_msgs::Point> waypoints_;
};

int main(int argc, char** argv)
{
  ros::init(argc, argv, "start_goal_marker_node");
  StartGoalMarkerNode node;
  node.spin();
  return 0;
}
