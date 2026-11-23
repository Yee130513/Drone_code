#ifndef _DYN_A_STAR_H_
#define _DYN_A_STAR_H_

#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <cctype>
#include <ros/ros.h>
#include <Eigen/Eigen>
#include <plan_env/grid_map.h>

#include "path_searching/dyn_a_star_sca.h"
#include "path_searching/dyn_a_star_a.h"

class AStar
{
private:
	AStarSCA::Ptr sca_planner_;
	AStarA::Ptr a_planner_;

	std::string planner_type_{"SCA_ASTAR"};
	bool use_sca_{true};

	GridMap::Ptr grid_map_{nullptr};
	Eigen::Vector3i pool_size_{Eigen::Vector3i(100, 100, 100)};
	bool has_init_map_{false};

	static std::string normalizeType(const std::string &raw)
	{
		std::string out;
		out.reserve(raw.size());
		for (char c : raw)
		{
			if (c == ' ' || c == '\t' || c == '-' || c == '_')
				continue;
			out.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
		}
		return out;
	}

	void ensureActivePlannerInitialized()
	{
		if (use_sca_)
		{
			if (!sca_planner_)
			{
				sca_planner_.reset(new AStarSCA());
			}
			if (has_init_map_ && grid_map_)
			{
				sca_planner_->initGridMap(grid_map_, pool_size_);
			}
		}
		else
		{
			if (!a_planner_)
			{
				a_planner_.reset(new AStarA());
			}
			if (has_init_map_ && grid_map_)
			{
				a_planner_->initGridMap(grid_map_, pool_size_);
			}
		}
	}

public:
	typedef std::shared_ptr<AStar> Ptr;

	AStar()
	{
		sca_planner_.reset(new AStarSCA());
		a_planner_.reset(new AStarA());
	}

	~AStar() = default;

	void setPlannerType(const std::string &planner_type)
	{
		planner_type_ = planner_type;
		const std::string norm = normalizeType(planner_type);

		if (norm == "SCAASTAR" || norm == "SCAA*" || norm == "SCA" || norm == "DYNASTARSCA")
		{
			use_sca_ = true;
		}
		else if (norm == "ASTAR" || norm == "A*" || norm == "A" || norm == "BASELINE" ||
				 norm == "ORIGINASTAR" || norm == "DYNASTARA")
		{
			use_sca_ = false;
		}
		else
		{
			ROS_WARN("[AStar] Unknown planner_type '%s', fallback to SCA_ASTAR (dyn_a_star_sca).",
					 planner_type.c_str());
			use_sca_ = true;
		}

		ensureActivePlannerInitialized();
	}

	inline std::string getPlannerType() const { return planner_type_; }
	inline bool isUsingSCA() const { return use_sca_; }

	void initGridMap(GridMap::Ptr occ_map, const Eigen::Vector3i pool_size)
	{
		grid_map_ = occ_map;
		pool_size_ = pool_size;
		has_init_map_ = true;
		ensureActivePlannerInitialized();

		ROS_INFO("[AStar] Front-end planner_type='%s' -> active implementation: %s",
				 planner_type_.c_str(),
				 use_sca_ ? "SCA-A* (dyn_a_star_sca)" : "Baseline A* (dyn_a_star_a)");
	}

	bool AstarSearch(const double step_size, Eigen::Vector3d start_pt, Eigen::Vector3d end_pt)
	{
		ensureActivePlannerInitialized();
		if (use_sca_)
		{
			return sca_planner_->AstarSearch(step_size, start_pt, end_pt);
		}
		return a_planner_->AstarSearch(step_size, start_pt, end_pt);
	}

	std::vector<Eigen::Vector3d> getPath()
	{
		if (use_sca_)
		{
			return sca_planner_->getPath();
		}
		return a_planner_->getPath();
	}

	inline int getExpandedNodes() const
	{
		return use_sca_ ? sca_planner_->getExpandedNodes() : a_planner_->getExpandedNodes();
	}

	inline double getPathLength() const
	{
		return use_sca_ ? sca_planner_->getPathLength() : a_planner_->getPathLength();
	}

	inline double getLastSearchTimeMs() const
	{
		return use_sca_ ? sca_planner_->getLastSearchTimeMs() : a_planner_->getLastSearchTimeMs();
	}

	inline double getLastMinClearance() const
	{
		return use_sca_ ? sca_planner_->getLastMinClearance() : a_planner_->getLastMinClearance();
	}

	inline bool usedFallback() const
	{
		return use_sca_ ? sca_planner_->usedFallback() : a_planner_->usedFallback();
	}

	inline std::string getLastFallbackReason() const
	{
		return use_sca_ ? sca_planner_->getLastFallbackReason() : a_planner_->getLastFallbackReason();
	}

	inline int getSegmentCheckFailCount() const
	{
		return use_sca_ ? sca_planner_->getSegmentCheckFailCount() : a_planner_->getSegmentCheckFailCount();
	}
};

#endif
