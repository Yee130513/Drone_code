#ifndef _JPS_THETA_PLANNER_H_
#define _JPS_THETA_PLANNER_H_

#include <iostream>
#include <ros/ros.h>
#include <Eigen/Eigen>
#include <plan_env/grid_map.h>
#include <queue>
#include <vector>

constexpr double inf = 1e20;
struct GridNode;
typedef GridNode *GridNodePtr;

struct GridNode
{
    int rounds{0}; 
    enum enum_state { UNDEFINED = 0, OPENSET = 1, CLOSEDSET = 2 };
    enum_state state{UNDEFINED};
    
    Eigen::Vector3i index;
    double gScore{inf}, fScore{inf};
    GridNodePtr cameFrom{nullptr};
};

class NodeComparator
{
public:
    bool operator()(GridNodePtr node1, GridNodePtr node2) {
        return node1->fScore > node2->fScore;
    }
};

class JPSThetaPlanner
{
private:
    GridMap::Ptr grid_map_;

    double step_size_, inv_step_size_;
    Eigen::Vector3d center_;
    Eigen::Vector3i CENTER_IDX_, POOL_SIZE_;
    
    std::vector<GridNodePtr> gridPath_;
    GridNodePtr ***GridNodeMap_;
    std::priority_queue<GridNodePtr, std::vector<GridNodePtr>, NodeComparator> openSet_;
    int rounds_{0};

    inline Eigen::Vector3d Index2Coord(const Eigen::Vector3i &index) const;
    inline bool Coord2Index(const Eigen::Vector3d &pt, Eigen::Vector3i &idx) const;
    inline bool checkOccupancy(const Eigen::Vector3d &pos) { return (bool)grid_map_->getInflateOccupancy(pos); }
    
    bool ConvertToIndexAndAdjustStartEndPoints(Eigen::Vector3d &start_pt, Eigen::Vector3d &end_pt, Eigen::Vector3i &start_idx, Eigen::Vector3i &end_idx);
    std::vector<GridNodePtr> retrievePath(GridNodePtr current);

    // Heuristic: Octile Distance
    double getOctileHeu(const Eigen::Vector3i &idx1, const Eigen::Vector3i &idx2);

    // Optimized Fast Line-of-Sight
    bool fastLineOfSight(const Eigen::Vector3d &pt1, const Eigen::Vector3d &pt2);

public:
    typedef std::shared_ptr<JPSThetaPlanner> Ptr;

    JPSThetaPlanner(){};
    ~JPSThetaPlanner();

    void initGridMap(GridMap::Ptr occ_map, const Eigen::Vector3i pool_size);
    bool AstarSearch(const double step_size, Eigen::Vector3d start_pt, Eigen::Vector3d end_pt);
    std::vector<Eigen::Vector3d> getPath();
};

#endif