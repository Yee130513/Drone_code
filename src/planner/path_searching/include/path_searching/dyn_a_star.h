#ifndef _DYN_A_STAR_H_
#define _DYN_A_STAR_H_

#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>
#include <ros/ros.h>
#include <ros/console.h>
#include <Eigen/Eigen>
#include <plan_env/grid_map.h>
#include <queue>

constexpr double inf = 1e18;

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

	int rounds{0};
	enum enum_state state{UNDEFINED};
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

	/* SCA-A* parameters: conservative defaults for front-end guide path. */
	int expand_step_{2};
	int dense_expand_step_{1};
	int max_probe_steps_{3};

	double epsilon_{1.4};
	double weight_clearance_{0.02};//0.03
	double weight_turn_{0.02};

	bool conservative_edge_check_{true};
	bool reject_unknown_{false};

	int max_expand_count_sca_{3000};
	int max_expand_count_fallback_{8000};

	double max_search_time_sca_{0.10};
	double max_search_time_fallback_{0.15};
	double los_prune_max_dist_{1.5};
	double dense_sample_resolution_ratio_{1.0};

	double min_path_clearance_limit_{0.0};

	std::string last_fallback_reason_{"NONE"};
	double last_search_time_ms_{0.0};
	double last_min_clearance_{0.0};
	int segment_check_fail_count_{0};
	bool used_fallback_{false};

	bool searchOnce(const Eigen::Vector3d &start_pt,
					const Eigen::Vector3d &end_pt,
					bool use_sca);

	bool edgeFeasible(const Eigen::Vector3i &idx1,
					  const Eigen::Vector3i &idx2,
					  bool conservative_check);

	bool segmentFreeSampled(const Eigen::Vector3d &p1,
							const Eigen::Vector3d &p2);

	bool lineOfSight(const Eigen::Vector3i &idx1,
					 const Eigen::Vector3i &idx2);

	double estimateClearance(const Eigen::Vector3d &pt);
	double getClearanceCost(const Eigen::Vector3d &pt);

	double getTurnPenalty(GridNodePtr current,
						  const Eigen::Vector3i &neighbor_idx) const;

	bool validateGridPath(const std::vector<GridNodePtr> &path,
						  double *min_clearance);

	double computePathLength(const std::vector<Eigen::Vector3d> &path) const;

	double getDiagHeu(GridNodePtr node1, GridNodePtr node2);
	double getManhHeu(GridNodePtr node1, GridNodePtr node2);
	double getEuclHeu(GridNodePtr node1, GridNodePtr node2);
	inline double getHeu(GridNodePtr node1, GridNodePtr node2);

	bool ConvertToIndexAndAdjustStartEndPoints(const Eigen::Vector3d start_pt,
											   const Eigen::Vector3d end_pt,
											   Eigen::Vector3i &start_idx,
											   Eigen::Vector3i &end_idx);

	inline Eigen::Vector3d Index2Coord(const Eigen::Vector3i &index) const;
	inline bool Coord2Index(const Eigen::Vector3d &pt, Eigen::Vector3i &idx) const;

	inline bool checkOccupancy(const Eigen::Vector3d &pos) { return (bool)grid_map_->getInflateOccupancy(pos); }

	std::vector<GridNodePtr> retrievePath(GridNodePtr current);

	double step_size_, inv_step_size_;
	Eigen::Vector3d center_;
	Eigen::Vector3i CENTER_IDX_, POOL_SIZE_;
	const double tie_breaker_ = 1.0 + 1.0 / 10000;

	std::vector<GridNodePtr> gridPath_;

	GridNodePtr ***GridNodeMap_;
	std::priority_queue<GridNodePtr, std::vector<GridNodePtr>, NodeComparator> openSet_;

	int rounds_{0};
	int expanded_nodes_{0};
	double path_length_{0.0};

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
	return tie_breaker_ * getDiagHeu(node1, node2) * step_size_;
}

inline Eigen::Vector3d AStar::Index2Coord(const Eigen::Vector3i &index) const
{
	return ((index - CENTER_IDX_).cast<double>() * step_size_) + center_;
};

inline bool AStar::Coord2Index(const Eigen::Vector3d &pt, Eigen::Vector3i &idx) const
{
	idx = ((pt - center_) * inv_step_size_ + Eigen::Vector3d(0.5, 0.5, 0.5)).cast<int>() + CENTER_IDX_;

	if (idx(0) < 0 || idx(0) >= POOL_SIZE_(0) ||
		idx(1) < 0 || idx(1) >= POOL_SIZE_(1) ||
		idx(2) < 0 || idx(2) >= POOL_SIZE_(2))
	{
		ROS_ERROR("Ran out of pool, index=%d %d %d", idx(0), idx(1), idx(2));
		return false;
	}

	return true;
};

#endif
