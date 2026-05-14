#ifndef _THETA_STAR_H_
#define _THETA_STAR_H_

#include <iostream>
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
    enum enum_state state { UNDEFINED };
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

class ThetaStar
{
private:
    GridMap::Ptr grid_map_;

    double getEuclHeu(GridNodePtr node1, GridNodePtr node2);
    inline double getHeu(GridNodePtr node1, GridNodePtr node2);

    bool ConvertToIndexAndAdjustStartEndPoints(const Eigen::Vector3d start_pt, const Eigen::Vector3d end_pt, Eigen::Vector3i &start_idx, Eigen::Vector3i &end_idx);

    inline Eigen::Vector3d Index2Coord(const Eigen::Vector3i &index) const;
    inline bool Coord2Index(const Eigen::Vector3d &pt, Eigen::Vector3i &idx) const;

    // 安全约束：严格调用膨胀图进行碰撞判定
    inline bool checkOccupancy(const Eigen::Vector3d &pos) { return (bool)grid_map_->getInflateOccupancy(pos); }
    
    // Theta* 核心逻辑：3D Line-of-Sight 检测
    bool lineOfSight(const Eigen::Vector3d &pt1, const Eigen::Vector3d &pt2);

    std::vector<GridNodePtr> retrievePath(GridNodePtr current);

    double step_size_, inv_step_size_;
    Eigen::Vector3d center_;
    Eigen::Vector3i CENTER_IDX_, POOL_SIZE_;
    const double tie_breaker_ = 1.0 + 1.0 / 10000;

    std::vector<GridNodePtr> gridPath_;
    GridNodePtr ***GridNodeMap_;
    std::priority_queue<GridNodePtr, std::vector<GridNodePtr>, NodeComparator> openSet_;

    int rounds_{0};

public:
    typedef std::shared_ptr<ThetaStar> Ptr;

    ThetaStar(){};
    ~ThetaStar();

    void initGridMap(GridMap::Ptr occ_map, const Eigen::Vector3i pool_size);
    bool AstarSearch(const double step_size, Eigen::Vector3d start_pt, Eigen::Vector3d end_pt); // 保持原名以便替换
    std::vector<Eigen::Vector3d> getPath();
};

inline double ThetaStar::getHeu(GridNodePtr node1, GridNodePtr node2)
{
    return tie_breaker_ * getEuclHeu(node1, node2);
}

inline Eigen::Vector3d ThetaStar::Index2Coord(const Eigen::Vector3i &index) const
{
    return ((index - CENTER_IDX_).cast<double>() * step_size_) + center_;
};

inline bool ThetaStar::Coord2Index(const Eigen::Vector3d &pt, Eigen::Vector3i &idx) const
{
    idx = ((pt - center_) * inv_step_size_ + Eigen::Vector3d(0.5, 0.5, 0.5)).cast<int>() + CENTER_IDX_;
    if (idx(0) < 0 || idx(0) >= POOL_SIZE_(0) || idx(1) < 0 || idx(1) >= POOL_SIZE_(1) || idx(2) < 0 || idx(2) >= POOL_SIZE_(2))
    {
        return false;
    }
    return true;
};

#endif