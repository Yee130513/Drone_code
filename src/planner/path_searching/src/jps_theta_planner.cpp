#include "path_searching/jps_theta_planner.h"

using namespace std;
using namespace Eigen;

JPSThetaPlanner::~JPSThetaPlanner()
{
    for (int i = 0; i < POOL_SIZE_(0); i++)
        for (int j = 0; j < POOL_SIZE_(1); j++)
            for (int k = 0; k < POOL_SIZE_(2); k++)
                delete GridNodeMap_[i][j][k];
}

void JPSThetaPlanner::initGridMap(GridMap::Ptr occ_map, const Eigen::Vector3i pool_size)
{
    POOL_SIZE_ = pool_size;
    CENTER_IDX_ = pool_size / 2;

    GridNodeMap_ = new GridNodePtr **[POOL_SIZE_(0)];
    for (int i = 0; i < POOL_SIZE_(0); i++) {
        GridNodeMap_[i] = new GridNodePtr *[POOL_SIZE_(1)];
        for (int j = 0; j < POOL_SIZE_(1); j++) {
            GridNodeMap_[i][j] = new GridNodePtr[POOL_SIZE_(2)];
            for (int k = 0; k < POOL_SIZE_(2); k++) {
                GridNodeMap_[i][j][k] = new GridNode;
            }
        }
    }
    grid_map_ = occ_map;
}

inline Eigen::Vector3d JPSThetaPlanner::Index2Coord(const Eigen::Vector3i &index) const {
    return ((index - CENTER_IDX_).cast<double>() * step_size_) + center_;
}

inline bool JPSThetaPlanner::Coord2Index(const Eigen::Vector3d &pt, Eigen::Vector3i &idx) const {
    idx = ((pt - center_) * inv_step_size_ + Eigen::Vector3d(0.5, 0.5, 0.5)).cast<int>() + CENTER_IDX_;
    if (idx(0) < 1 || idx(0) >= POOL_SIZE_(0) - 1 || 
        idx(1) < 1 || idx(1) >= POOL_SIZE_(1) - 1 || 
        idx(2) < 1 || idx(2) >= POOL_SIZE_(2) - 1)
        return false;
    return true;
}

// ---------------- 极致优化的 3D 视线检测 ----------------
bool JPSThetaPlanner::fastLineOfSight(const Eigen::Vector3d &pt1, const Eigen::Vector3d &pt2)
{
    Eigen::Vector3d diff = pt2 - pt1;
    double dist = diff.norm();
    if (dist < 1e-4) return true;
    
    // 采用增量加法替代原版乘法，极大地节省 CPU 耗时，步长为栅格的一半以防穿模
    int steps = std::max(1, (int)std::ceil(dist / (step_size_ * 0.5))); 
    Eigen::Vector3d step_vec = diff / (double)steps;
    Eigen::Vector3d curr_pt = pt1 + step_vec;

    for (int i = 1; i < steps; ++i) {
        if (checkOccupancy(curr_pt)) return false;
        curr_pt += step_vec;
    }
    return true;
}

// ---------------- 防死锁的起点/终点校验 ----------------
bool JPSThetaPlanner::ConvertToIndexAndAdjustStartEndPoints(Vector3d &start_pt, Vector3d &end_pt, Vector3i &start_idx, Vector3i &end_idx)
{
    // 【重要防护】防止起点终点过近导致 normalize() 产生 NaN 崩溃
    if ((start_pt - end_pt).norm() < 1e-3) return false; 

    if (!Coord2Index(start_pt, start_idx) || !Coord2Index(end_pt, end_idx)) return false;

    int max_push_iter = 50; 
    int iter = 0;
    while (checkOccupancy(Index2Coord(start_idx)) && iter++ < max_push_iter) {
        start_pt = (start_pt - end_pt).normalized() * step_size_ + start_pt;
        if (!Coord2Index(start_pt, start_idx)) return false;
    }

    iter = 0;
    while (checkOccupancy(Index2Coord(end_idx)) && iter++ < max_push_iter) {
        end_pt = (end_pt - start_pt).normalized() * step_size_ + end_pt;
        if (!Coord2Index(end_pt, end_idx)) return false;
    }

    if (iter >= max_push_iter) return false; 
    return true;
}

// 【优化】Theta* 是任意角度拉直线，用精确的欧氏距离作为 Heuristic 最合理
double JPSThetaPlanner::getOctileHeu(const Eigen::Vector3i &idx1, const Eigen::Vector3i &idx2)
{
    return (idx2 - idx1).cast<double>().norm() * 1.0001; // 1.0001 为 EGO 原版的 Tie-breaker
}

std::vector<GridNodePtr> JPSThetaPlanner::retrievePath(GridNodePtr current)
{
    std::vector<GridNodePtr> path;
    path.push_back(current);
    while (current->cameFrom != NULL) {
        current = current->cameFrom;
        path.push_back(current);
    }
    return path;
}

// ---------------- 核心：带 LoS 剪枝的极速寻路主循环 ----------------
bool JPSThetaPlanner::AstarSearch(const double step_size, Vector3d start_pt, Vector3d end_pt)
{
    ros::Time time_1 = ros::Time::now();
    ++rounds_;
    step_size_ = step_size;
    inv_step_size_ = 1.0 / step_size;
    center_ = (start_pt + end_pt) / 2.0;

    Vector3i start_idx, end_idx;
    if (!ConvertToIndexAndAdjustStartEndPoints(start_pt, end_pt, start_idx, end_idx)) {
        return false;
    }

    if (start_idx == end_idx) return false; 

    GridNodePtr startPtr = GridNodeMap_[start_idx(0)][start_idx(1)][start_idx(2)];
    GridNodePtr endPtr = GridNodeMap_[end_idx(0)][end_idx(1)][end_idx(2)];

    std::priority_queue<GridNodePtr, std::vector<GridNodePtr>, NodeComparator> empty;
    openSet_.swap(empty);

    startPtr->index = start_idx;
    startPtr->rounds = rounds_;
    startPtr->gScore = 0.0;
    startPtr->fScore = getOctileHeu(start_idx, end_idx);
    startPtr->state = GridNode::OPENSET;
    startPtr->cameFrom = NULL;
    openSet_.push(startPtr);

    endPtr->index = end_idx;

    GridNodePtr current = NULL;
    GridNodePtr neighborPtr = NULL;

    while (!openSet_.empty())
    {
        current = openSet_.top();
        openSet_.pop();

        // 【致命 Bug 修复处】：绝对不能重复扩展已经被关闭的节点旧副本！
        if (current->state == GridNode::CLOSEDSET) continue; 

        if (current->index == endPtr->index) {
            gridPath_ = retrievePath(current);
            return true;
        }
        current->state = GridNode::CLOSEDSET; // 标记为已扩展

        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                for (int dz = -1; dz <= 1; dz++) {
                    if (dx == 0 && dy == 0 && dz == 0) continue;

                    Vector3i neighborIdx = current->index + Vector3i(dx, dy, dz);
                    if (neighborIdx(0) < 1 || neighborIdx(0) >= POOL_SIZE_(0) - 1 || 
                        neighborIdx(1) < 1 || neighborIdx(1) >= POOL_SIZE_(1) - 1 || 
                        neighborIdx(2) < 1 || neighborIdx(2) >= POOL_SIZE_(2) - 1)
                        continue;

                    neighborPtr = GridNodeMap_[neighborIdx(0)][neighborIdx(1)][neighborIdx(2)];
                    
                    // 跳过已经扩展过的节点
                    if (neighborPtr->rounds == rounds_ && neighborPtr->state == GridNode::CLOSEDSET) continue;

                    neighborPtr->index = neighborIdx;
                    neighborPtr->rounds = rounds_;

                    if (checkOccupancy(Index2Coord(neighborPtr->index))) continue;

                    // === Any-Angle Theta* 松弛与前置剪枝 (大幅提速) ===
                    double static_cost = sqrt(dx*dx + dy*dy + dz*dz);
                    double tentative_gScore = current->gScore + static_cost;
                    GridNodePtr tentative_parent = current;

                    if (current->cameFrom != NULL) {
                        Vector3d parent_pos = Index2Coord(current->cameFrom->index);
                        Vector3d neighbor_pos = Index2Coord(neighborPtr->index);
                        
                        // 【提速核心】：只有当直线跨越的代价确定变小时，才去执行昂贵的射线检测！
                        double g_via_parent = current->cameFrom->gScore + (neighbor_pos - parent_pos).norm() * inv_step_size_;
                        if (g_via_parent < tentative_gScore) {
                            if (fastLineOfSight(parent_pos, neighbor_pos)) {
                                tentative_gScore = g_via_parent;
                                tentative_parent = current->cameFrom;
                            }
                        }
                    }

                    if (neighborPtr->state != GridNode::OPENSET) {
                        neighborPtr->state = GridNode::OPENSET;
                        neighborPtr->cameFrom = tentative_parent;
                        neighborPtr->gScore = tentative_gScore;
                        neighborPtr->fScore = tentative_gScore + getOctileHeu(neighborPtr->index, end_idx);
                        openSet_.push(neighborPtr);
                    }
                    else if (tentative_gScore < neighborPtr->gScore) {
                        neighborPtr->cameFrom = tentative_parent;
                        neighborPtr->gScore = tentative_gScore;
                        neighborPtr->fScore = tentative_gScore + getOctileHeu(neighborPtr->index, end_idx);
                        openSet_.push(neighborPtr); // 重新压入队列以更新优先级 (旧副本会在 pop 时被忽略)
                    }
                }
            }
        }

        // 超时保护
        if ((ros::Time::now() - time_1).toSec() > 0.2) {
            ROS_ERROR("[Theta*] Timeout! Map is too complex or entirely blocked.");
            return false;
        }
    }
    return false;
}

// ---------------- 强制稠密插值 (迎合 B-Spline 投影) ----------------
vector<Vector3d> JPSThetaPlanner::getPath()
{
    vector<Vector3d> sparse_path;
    for (auto ptr : gridPath_)
        sparse_path.push_back(Index2Coord(ptr->index));

    reverse(sparse_path.begin(), sparse_path.end());

    vector<Vector3d> dense_path;
    if (sparse_path.empty()) return dense_path;

    dense_path.push_back(sparse_path[0]); 

    for (size_t i = 1; i < sparse_path.size(); ++i) {
        Vector3d p1 = sparse_path[i - 1];
        Vector3d p2 = sparse_path[i];
        
        double dist = (p2 - p1).norm();
        // 强制 0.1m 稠密插值
        int num_segments = std::max(1, (int)std::ceil(dist / step_size_)); 

        for (int j = 1; j <= num_segments; ++j) {
            dense_path.push_back(p1 + (p2 - p1) * ((double)j / num_segments));
        }
    }
    return dense_path;
}