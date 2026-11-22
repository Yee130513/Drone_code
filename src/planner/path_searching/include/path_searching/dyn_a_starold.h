#ifndef _DYN_A_STAR_H_
#define _DYN_A_STAR_H_

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <memory>
#include <ros/ros.h>
#include <ros/console.h>
#include <Eigen/Eigen>
#include <plan_env/grid_map.h>
#include <queue>

constexpr double inf = 1 >> 20;
struct GridNode;
typedef GridNode *GridNodePtr;

struct GridNode
{
	enum enum_state
	{
		OPENSET = 1,
		CLOSEDSET = 2,
		UNDEFINED = 3
	};

	int rounds{0}; // Distinguish every call
	enum enum_state state
	{
		UNDEFINED
	};
	Eigen::Vector3i index;

	double gScore{inf}, fScore{inf};
	GridNodePtr cameFrom{NULL};
};

class NodeComparator
{
public:
	bool operator()(GridNodePtr node1, GridNodePtr node2)
	{
		return node1->fScore > node2->fScore;
	}
};

class AStar
{
private:
	GridMap::Ptr grid_map_;

	// ===== Ablation / logging statistics, compatible with SCA-A* version =====
	std::string last_fallback_reason_{"NONE"};
	double last_search_time_ms_{0.0};
	double last_min_clearance_{0.0};
	int segment_check_fail_count_{0};
	bool used_fallback_{false};

	int expanded_nodes_{0};
	double path_length_{0.0};
	int edge_reject_count_{0};
	int output_edge_fail_count_{0};

	// Metric-only parameters. These must not affect vanilla A* search behavior.
	int max_probe_steps_{3};
	bool reject_unknown_for_metric_{false};

	inline void coord2gridIndexFast(const double x, const double y, const double z, int &id_x, int &id_y, int &id_z);

	double getDiagHeu(GridNodePtr node1, GridNodePtr node2);
	double getManhHeu(GridNodePtr node1, GridNodePtr node2);
	double getEuclHeu(GridNodePtr node1, GridNodePtr node2);
	inline double getHeu(GridNodePtr node1, GridNodePtr node2);

	bool ConvertToIndexAndAdjustStartEndPoints(const Eigen::Vector3d start_pt,
											   const Eigen::Vector3d end_pt,
											   Eigen::Vector3i &start_idx,
											   Eigen::Vector3i &end_idx,
											   Eigen::Vector3d *adjusted_start_pt,
											   Eigen::Vector3d *adjusted_end_pt);

	inline Eigen::Vector3d Index2Coord(const Eigen::Vector3i &index) const;
	inline bool Coord2Index(const Eigen::Vector3d &pt, Eigen::Vector3i &idx) const;

	//bool (*checkOccupancyPtr)( const Eigen::Vector3d &pos );

	inline bool checkOccupancy(const Eigen::Vector3d &pos) { return (bool)grid_map_->getInflateOccupancy(pos); }

	std::vector<GridNodePtr> retrievePath(GridNodePtr current);
	double estimateClearanceMetricOnly(const Eigen::Vector3d &pt);
	bool segmentFreeSampledMetricOnly(const Eigen::Vector3d &p1,
									  const Eigen::Vector3d &p2,
									  bool ignore_start_occupancy = false,
									  bool ignore_end_occupancy = false);
	double computePathLength(const std::vector<Eigen::Vector3d> &path) const;
	void updatePathMetricsFromPath(const std::vector<Eigen::Vector3d> &path);
	void updatePathMetricsFromGridPath();
	std::vector<Eigen::Vector3d> buildFinalPathFromGridPath() const;
	bool validateFinalPathAndUpdateMetrics(const std::vector<Eigen::Vector3d> &path,
										   std::string *reason);

	double step_size_, inv_step_size_;
	Eigen::Vector3d center_;
	Eigen::Vector3i CENTER_IDX_, POOL_SIZE_;
	const double tie_breaker_ = 1.0 + 1.0 / 10000;

	std::vector<GridNodePtr> gridPath_;
	std::vector<Eigen::Vector3d> final_path_;

	GridNodePtr ***GridNodeMap_;
	std::priority_queue<GridNodePtr, std::vector<GridNodePtr>, NodeComparator> openSet_;

	int rounds_{0};

	Eigen::Vector3d last_query_start_pt_{Eigen::Vector3d::Zero()};
	Eigen::Vector3d last_query_end_pt_{Eigen::Vector3d::Zero()};
	Eigen::Vector3d last_adjusted_start_pt_{Eigen::Vector3d::Zero()};
	Eigen::Vector3d last_adjusted_end_pt_{Eigen::Vector3d::Zero()};
	bool has_last_query_{false};
	std::string last_path_invalid_reason_{"NONE"};

public:
	typedef std::shared_ptr<AStar> Ptr;

	AStar(){};
	~AStar();

	void initGridMap(GridMap::Ptr occ_map, const Eigen::Vector3i pool_size);

	bool AstarSearch(const double step_size, Eigen::Vector3d start_pt, Eigen::Vector3d end_pt);

	std::vector<Eigen::Vector3d> getPath();
	inline int getExpandedNodes() const { return expanded_nodes_; }
	inline double getPathLength() const { return path_length_; }
	inline double getLastSearchTimeMs() const { return last_search_time_ms_; }
	inline double getLastMinClearance() const { return last_min_clearance_; }
	inline bool usedFallback() const { return used_fallback_; }
	inline std::string getLastFallbackReason() const { return last_fallback_reason_; }
	inline int getSegmentCheckFailCount() const { return segment_check_fail_count_; }
};

inline double AStar::getHeu(GridNodePtr node1, GridNodePtr node2)
{
	return tie_breaker_ * getDiagHeu(node1, node2);
}

inline Eigen::Vector3d AStar::Index2Coord(const Eigen::Vector3i &index) const
{
	return ((index - CENTER_IDX_).cast<double>() * step_size_) + center_;
};

inline bool AStar::Coord2Index(const Eigen::Vector3d &pt, Eigen::Vector3i &idx) const
{
	idx = ((pt - center_) * inv_step_size_ + Eigen::Vector3d(0.5, 0.5, 0.5)).cast<int>() + CENTER_IDX_;

	if (idx(0) < 0 || idx(0) >= POOL_SIZE_(0) || idx(1) < 0 || idx(1) >= POOL_SIZE_(1) || idx(2) < 0 || idx(2) >= POOL_SIZE_(2))
	{
		ROS_ERROR("Ran out of pool, index=%d %d %d", idx(0), idx(1), idx(2));
		return false;
	}

	return true;
};

#endif
