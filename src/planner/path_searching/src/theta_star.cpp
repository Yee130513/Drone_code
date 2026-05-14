#include "path_searching/theta_star.h"

using namespace std;
using namespace Eigen;

ThetaStar::~ThetaStar()
{
    for (int i = 0; i < POOL_SIZE_(0); i++)
        for (int j = 0; j < POOL_SIZE_(1); j++)
            for (int k = 0; k < POOL_SIZE_(2); k++)
                delete GridNodeMap_[i][j][k];
}

void ThetaStar::initGridMap(GridMap::Ptr occ_map, const Eigen::Vector3i pool_size)
{
    POOL_SIZE_ = pool_size;
    CENTER_IDX_ = pool_size / 2;

    GridNodeMap_ = new GridNodePtr **[POOL_SIZE_(0)];
    for (int i = 0; i < POOL_SIZE_(0); i++)
    {
        GridNodeMap_[i] = new GridNodePtr *[POOL_SIZE_(1)];
        for (int j = 0; j < POOL_SIZE_(1); j++)
        {
            GridNodeMap_[i][j] = new GridNodePtr[POOL_SIZE_(2)];
            for (int k = 0; k < POOL_SIZE_(2); k++)
            {
                GridNodeMap_[i][j][k] = new GridNode;
            }
        }
    }
    grid_map_ = occ_map;
}

double ThetaStar::getEuclHeu(GridNodePtr node1, GridNodePtr node2)
{
    return (node2->index - node1->index).norm();
}

vector<GridNodePtr> ThetaStar::retrievePath(GridNodePtr current)
{
    vector<GridNodePtr> path;
    path.push_back(current);
    while (current->cameFrom != NULL)
    {
        current = current->cameFrom;
        path.push_back(current);
    }
    return path;
}

// ---------------- 核心 1：3D Line-of-Sight 视线检测 ----------------
bool ThetaStar::lineOfSight(const Eigen::Vector3d &pt1, const Eigen::Vector3d &pt2)
{
    Eigen::Vector3d dir = pt2 - pt1;
    double dist = dir.norm();
    if (dist < 1e-4) return true;
    
    dir /= dist;
    // 采样步长设为体素的一半，保证绝对不会穿模漏检
    double step = step_size_ / 2.0; 
    for (double d = step; d < dist; d += step)
    {
        // 严格遵守规范：仅查询 getInflateOccupancy
        if (checkOccupancy(pt1 + dir * d)) 
            return false;
    }
    return true;
}

bool ThetaStar::ConvertToIndexAndAdjustStartEndPoints(Vector3d start_pt, Vector3d end_pt, Vector3i &start_idx, Vector3i &end_idx)
{
    if (!Coord2Index(start_pt, start_idx) || !Coord2Index(end_pt, end_idx))
        return false;

    if (checkOccupancy(Index2Coord(start_idx))) {
        do {
            start_pt = (start_pt - end_pt).normalized() * step_size_ + start_pt;
            if (!Coord2Index(start_pt, start_idx)) return false;
        } while (checkOccupancy(Index2Coord(start_idx)));
    }

    if (checkOccupancy(Index2Coord(end_idx))) {
        do {
            end_pt = (end_pt - start_pt).normalized() * step_size_ + end_pt;
            if (!Coord2Index(end_pt, end_idx)) return false;
        } while (checkOccupancy(Index2Coord(end_idx)));
    }
    return true;
}

// 注意：这里保留了 AstarSearch 的方法名，以避免你修改上层 bspline_optimizer.cpp 太多代码
bool ThetaStar::AstarSearch(const double step_size, Vector3d start_pt, Vector3d end_pt)
{
    ros::Time time_1 = ros::Time::now();
    ++rounds_;

    step_size_ = step_size;
    inv_step_size_ = 1 / step_size;
    center_ = (start_pt + end_pt) / 2;

    Vector3i start_idx, end_idx;
    if (!ConvertToIndexAndAdjustStartEndPoints(start_pt, end_pt, start_idx, end_idx))
    {
        ROS_ERROR("Unable to handle the initial or end point, force return!");
        return false;
    }

    GridNodePtr startPtr = GridNodeMap_[start_idx(0)][start_idx(1)][start_idx(2)];
    GridNodePtr endPtr = GridNodeMap_[end_idx(0)][end_idx(1)][end_idx(2)];

    std::priority_queue<GridNodePtr, std::vector<GridNodePtr>, NodeComparator> empty;
    openSet_.swap(empty);

    startPtr->index = start_idx;
    startPtr->rounds = rounds_;
    startPtr->gScore = 0;
    startPtr->fScore = getHeu(startPtr, endPtr);
    startPtr->state = GridNode::OPENSET;
    startPtr->cameFrom = NULL;
    openSet_.push(startPtr);

    endPtr->index = end_idx;

    GridNodePtr neighborPtr = NULL;
    GridNodePtr current = NULL;
    int num_iter = 0;

    while (!openSet_.empty())
    {
        num_iter++;
        current = openSet_.top();
        openSet_.pop();

        if (current->index == endPtr->index)
        {
            gridPath_ = retrievePath(current);
            return true;
        }
        current->state = GridNode::CLOSEDSET;

        for (int dx = -1; dx <= 1; dx++)
            for (int dy = -1; dy <= 1; dy++)
                for (int dz = -1; dz <= 1; dz++)
                {
                    if (dx == 0 && dy == 0 && dz == 0) continue;

                    Vector3i neighborIdx = current->index + Vector3i(dx, dy, dz);
                    if (neighborIdx(0) < 1 || neighborIdx(0) >= POOL_SIZE_(0) - 1 || neighborIdx(1) < 1 || neighborIdx(1) >= POOL_SIZE_(1) - 1 || neighborIdx(2) < 1 || neighborIdx(2) >= POOL_SIZE_(2) - 1)
                        continue;

                    neighborPtr = GridNodeMap_[neighborIdx(0)][neighborIdx(1)][neighborIdx(2)];
                    neighborPtr->index = neighborIdx;

                    bool flag_explored = neighborPtr->rounds == rounds_;
                    if (flag_explored && neighborPtr->state == GridNode::CLOSEDSET) continue;

                    neighborPtr->rounds = rounds_;
                    if (checkOccupancy(Index2Coord(neighborPtr->index))) continue;

                    // ---------------- 核心 2：Theta* 的松弛操作 (Relaxation) ----------------
                    double tentative_gScore;
                    GridNodePtr tentative_parent = current;

                    // 尝试将邻居节点直接与 current 的“父节点”连接
                    if (current->cameFrom != NULL) {
                        Vector3d parent_pos = Index2Coord(current->cameFrom->index);
                        Vector3d neighbor_pos = Index2Coord(neighborPtr->index);
                        
                        // 如果父节点到邻居之间无障碍 (LoS)，则直接跨越 current 拉直线！
                        if (lineOfSight(parent_pos, neighbor_pos)) {
                            tentative_parent = current->cameFrom;
                        }
                    }

                    // 计算欧式距离作为 gScore
                    double dist_to_parent = (neighborPtr->index - tentative_parent->index).norm();
                    tentative_gScore = tentative_parent->gScore + dist_to_parent;

                    if (!flag_explored)
                    {
                        neighborPtr->state = GridNode::OPENSET;
                        neighborPtr->cameFrom = tentative_parent; // 可能直接指向爷爷节点
                        neighborPtr->gScore = tentative_gScore;
                        neighborPtr->fScore = tentative_gScore + getHeu(neighborPtr, endPtr);
                        openSet_.push(neighborPtr);
                    }
                    else if (tentative_gScore < neighborPtr->gScore)
                    {
                        neighborPtr->cameFrom = tentative_parent;
                        neighborPtr->gScore = tentative_gScore;
                        neighborPtr->fScore = tentative_gScore + getHeu(neighborPtr, endPtr);
                        // 重新推入优先队列 (优化手段)
                        openSet_.push(neighborPtr); 
                    }
                }
        
        // 超时保护保持原版设定
        ros::Time time_2 = ros::Time::now();
        if ((time_2 - time_1).toSec() > 0.2) {
            ROS_WARN("Failed in Theta* path searching !!! 0.2s time limit exceeded.");
            return false;
        }
    }
    return false;
}

// ---------------- 核心 3：强制稠密插值（拯救 B 样条投影器） ----------------
vector<Vector3d> ThetaStar::getPath()
{
    vector<Vector3d> sparse_path;
    for (auto ptr : gridPath_)
        sparse_path.push_back(Index2Coord(ptr->index));

    reverse(sparse_path.begin(), sparse_path.end());

    vector<Vector3d> dense_path;
    if (sparse_path.empty()) return dense_path;

    dense_path.push_back(sparse_path[0]); // 起点

    for (size_t i = 1; i < sparse_path.size(); ++i)
    {
        Vector3d p1 = sparse_path[i - 1];
        Vector3d p2 = sparse_path[i];
        
        double dist = (p2 - p1).norm();
        // 以 step_size_ (通常是 0.1m) 为分辨率强行补点
        int num_segments = std::max(1, (int)std::ceil(dist / step_size_)); 

        for (int j = 1; j <= num_segments; ++j)
        {
            // 线性插值
            dense_path.push_back(p1 + (p2 - p1) * ((double)j / num_segments));
        }
    }

    return dense_path;
}