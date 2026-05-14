#include "path_searching/dyn_a_star.h"

using namespace std;
using namespace Eigen;

AStar::~AStar()
{
    for (int i = 0; i < POOL_SIZE_(0); i++)
        for (int j = 0; j < POOL_SIZE_(1); j++)
            for (int k = 0; k < POOL_SIZE_(2); k++)
                delete GridNodeMap_[i][j][k];
}

void AStar::initGridMap(GridMap::Ptr occ_map, const Eigen::Vector3i pool_size)
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

// ==================== dyn_a_star.cpp ====================
// 在合适位置插入此内联函数

inline int AStar::getClearancePenalty(const Eigen::Vector3d &pt)
{
    int penalty = 0;
    // 物理探测距离
    double d = probe_dist_ * step_size_; 
    
    // 预定义 6 个方向的物理偏移量 (避免循环内动态创建 Eigen 对象)
    const double offsets[6][3] = {
        {d, 0, 0}, {-d, 0, 0}, 
        {0, d, 0}, {0, -d, 0}, 
        {0, 0, d}, {0, 0, -d}
    };

    Eigen::Vector3i probe_idx;
    for (int i = 0; i < 6; ++i)
    {
        Eigen::Vector3d probe_pt(pt(0) + offsets[i][0], pt(1) + offsets[i][1], pt(2) + offsets[i][2]);
        
        // 🚨 架构师约束：越界检查！
        // Coord2Index 内部已包含边界检查，如果不合法 (越界)，则当做安全区域 (不增加惩罚)
        if (!Coord2Index(probe_pt, probe_idx)) {
            continue; 
        }

        // 探测膨胀障碍物图
        if (checkOccupancy(probe_pt)) {
            penalty++;
        }
    }
    return penalty;
}

bool AStar::lineOfSight(const Eigen::Vector3i &idx1, const Eigen::Vector3i &idx2)
{
    Eigen::Vector3d p1 = Index2Coord(idx1);
    Eigen::Vector3d p2 = Index2Coord(idx2);
    Eigen::Vector3d dir = p2 - p1;
    double dist = dir.norm();
    
    if (dist < 1e-4) return true;
    dir /= dist;
    
    // 采用 0.5 倍栅格分辨率步进，保证绝对不会漏检薄墙
    double step = step_size_ / 2.0; 
    for (double d = step; d < dist; d += step)
    {
        if (checkOccupancy(p1 + dir * d)) 
            return false; // 中间有障碍物，视线被阻挡
    }
    return true;
}

double AStar::getDiagHeu(GridNodePtr node1, GridNodePtr node2)
{
    double dx = abs(node1->index(0) - node2->index(0));
    double dy = abs(node1->index(1) - node2->index(1));
    double dz = abs(node1->index(2) - node2->index(2));

    double h = 0.0;
    int diag = min(min(dx, dy), dz);
    dx -= diag;
    dy -= diag;
    dz -= diag;

    if (dx == 0)
    {
        h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dy, dz) + 1.0 * abs(dy - dz);
    }
    if (dy == 0)
    {
        h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dx, dz) + 1.0 * abs(dx - dz);
    }
    if (dz == 0)
    {
        h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dx, dy) + 1.0 * abs(dx - dy);
    }
    return h;
}

double AStar::getManhHeu(GridNodePtr node1, GridNodePtr node2)
{
    double dx = abs(node1->index(0) - node2->index(0));
    double dy = abs(node1->index(1) - node2->index(1));
    double dz = abs(node1->index(2) - node2->index(2));

    return dx + dy + dz;
}

double AStar::getEuclHeu(GridNodePtr node1, GridNodePtr node2)
{
    return (node2->index - node1->index).norm();
}

vector<GridNodePtr> AStar::retrievePath(GridNodePtr current)
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

bool AStar::ConvertToIndexAndAdjustStartEndPoints(Vector3d start_pt, Vector3d end_pt, Vector3i &start_idx, Vector3i &end_idx)
{
    if (!Coord2Index(start_pt, start_idx) || !Coord2Index(end_pt, end_idx))
        return false;

    if (checkOccupancy(Index2Coord(start_idx)))
    {
        //ROS_WARN("Start point is insdide an obstacle.");
        do
        {
            start_pt = (start_pt - end_pt).normalized() * step_size_ + start_pt;
            if (!Coord2Index(start_pt, start_idx))
                return false;
        } while (checkOccupancy(Index2Coord(start_idx)));
    }

    if (checkOccupancy(Index2Coord(end_idx)))
    {
        //ROS_WARN("End point is insdide an obstacle.");
        do
        {
            end_pt = (end_pt - start_pt).normalized() * step_size_ + end_pt;
            if (!Coord2Index(end_pt, end_idx))
                return false;
        } while (checkOccupancy(Index2Coord(end_idx)));
    }

    return true;
}
bool AStar::AstarSearch(const double step_size, Vector3d start_pt, Vector3d end_pt)
{
    expanded_nodes_ = 0; // 每次搜索前清零
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

    // 🚨 修复点 1：补全所有指针的声明！
    GridNodePtr startPtr = GridNodeMap_[start_idx(0)][start_idx(1)][start_idx(2)];
    GridNodePtr endPtr = GridNodeMap_[end_idx(0)][end_idx(1)][end_idx(2)];
    GridNodePtr current = NULL;
    GridNodePtr neighborPtr = NULL;

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
    // ... [前置初始化代码保持原样] ...
    // 稀疏 A* 正常绝不会超过 2000 次迭代。一旦超过，说明陷入死区或目标不可达，立刻熔断，避免引发内存溢出和 0.2s 线程卡死。
    

    int num_iter = 0;
    while (!openSet_.empty())
    {
        expanded_nodes_++; // 埋点 1：记录真实出队的节点数
        num_iter++;

        if (num_iter > 3000) {
        ROS_WARN("[SCA-A*] Triggered Kill-switch! Goal unreachable or stuck in local minima.");
        return false;
    }
        current = openSet_.top();
        openSet_.pop();

        // 终点放宽判定 (不变)
        double dx_goal = (double)(current->index(0) - endPtr->index(0));
        double dy_goal = (double)(current->index(1) - endPtr->index(1));
        double dz_goal = (double)(current->index(2) - endPtr->index(2));
        double dist_to_goal = sqrt(dx_goal*dx_goal + dy_goal*dy_goal + dz_goal*dz_goal);

        if (dist_to_goal <= expand_step_ * 1.74) 
        {
            // O(1) 终点防穿墙检查
            Vector3i mid_to_goal = (current->index + endPtr->index) / 2;
            if (!checkOccupancy(Index2Coord(mid_to_goal))) 
            {
                endPtr->cameFrom = current;
                gridPath_ = retrievePath(endPtr);
                return true;
            }
        }

        current->state = GridNode::CLOSEDSET;

        // =================================================================
        // 🚀 极速引擎 2：提取父节点方向动量，用于无掉头剪枝
        // =================================================================
        Vector3d dir_from_parent = Vector3d::Zero();
        if (current->cameFrom != NULL) {
            dir_from_parent = (current->index - current->cameFrom->index).cast<double>().normalized();
        }

        // =================================================================
        // 🚀 极速引擎 3：计算动态启发式权重 (离目标越近，贪婪度越低)
        // =================================================================
        double dynamic_epsilon = 1.0 + (epsilon_ - 1.0) * (dist_to_goal / (POOL_SIZE_(0) / 2.0));
        dynamic_epsilon = std::max(1.1, std::min(dynamic_epsilon, epsilon_));

        for (int dx_i = -expand_step_; dx_i <= expand_step_; dx_i += expand_step_)
            for (int dy_i = -expand_step_; dy_i <= expand_step_; dy_i += expand_step_)
                for (int dz_i = -expand_step_; dz_i <= expand_step_; dz_i += expand_step_)
                {
                    if (dx_i == 0 && dy_i == 0 && dz_i == 0) continue;

                    // 🚀 极速引擎 2 (续)：无掉头剪枝 (No U-Turn Pruning)
                    // 如果当前有运动方向，且尝试扩展的方向与来时方向相反（点积 < -0.5，夹角大于120度），直接剪枝！
                    Vector3d dir_to_neighbor = Vector3d(dx_i, dy_i, dz_i).normalized();
                    if (current->cameFrom != NULL && dir_from_parent.dot(dir_to_neighbor) < -0.5) {
                        continue; 
                    }

                    Vector3i neighborIdx(current->index(0) + dx_i, current->index(1) + dy_i, current->index(2) + dz_i);

                    // 边界检查
                    if (neighborIdx(0) < 2 || neighborIdx(0) >= POOL_SIZE_(0) - 2 || 
                        neighborIdx(1) < 2 || neighborIdx(1) >= POOL_SIZE_(1) - 2 || 
                        neighborIdx(2) < 2 || neighborIdx(2) >= POOL_SIZE_(2) - 2)
                        continue;

                    neighborPtr = GridNodeMap_[neighborIdx(0)][neighborIdx(1)][neighborIdx(2)];
                    neighborPtr->index = neighborIdx;

                    bool flag_explored = neighborPtr->rounds == rounds_;
                    if (flag_explored && neighborPtr->state == GridNode::CLOSEDSET) continue;

                    neighborPtr->rounds = rounds_;
                    Vector3d neighbor_pos = Index2Coord(neighborPtr->index);
                    
                    if (checkOccupancy(neighbor_pos)) continue;

                    // =================================================================
                    // 🚀 极速引擎 1：O(1) 极简防穿墙检测 (Midpoint Check)
                    // 完全废弃昂贵的 lineOfSight 射线遍历，直接查两点的数学中点！
                    // =================================================================
                    Vector3i mid_idx = (current->index + neighborIdx) / 2;
                    if (checkOccupancy(Index2Coord(mid_idx))) continue; 

                    // 机制 3: 间隙探测
                    int c_penalty = getClearancePenalty(neighbor_pos);
                    
                    double static_cost = sqrt(dx_i * dx_i + dy_i * dy_i + dz_i * dz_i) * step_size_;
                    double tentative_gScore = current->gScore + static_cost;

                    if (!flag_explored || tentative_gScore < neighborPtr->gScore)
                    {
                        neighborPtr->state = GridNode::OPENSET;
                        neighborPtr->cameFrom = current;
                        neighborPtr->gScore = tentative_gScore;
                        
                        // 使用动态 Epsilon 进行加速
                        neighborPtr->fScore = tentative_gScore + 
                                              dynamic_epsilon * getHeu(neighborPtr, endPtr) + 
                                              weight_clearance_ * c_penalty;
                        openSet_.push(neighborPtr);
                    }
                }
        
        // 超时保护
        ros::Time time_2 = ros::Time::now();
        if ((time_2 - time_1).toSec() > 0.2)
        {
            ROS_WARN("Failed in Fast-SCA-A* !!! 0.2s time limit exceeded. Iters: %d", num_iter);
            return false;
        }
    }
    return false;
}
/*bool AStar::AstarSearch(const double step_size, Vector3d start_pt, Vector3d end_pt)
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

    // if ( start_pt(0) > -1 && start_pt(0) < 0 )
    //     cout << "start_pt=" << start_pt.transpose() << " end_pt=" << end_pt.transpose() << endl;

    GridNodePtr startPtr = GridNodeMap_[start_idx(0)][start_idx(1)][start_idx(2)];
    GridNodePtr endPtr = GridNodeMap_[end_idx(0)][end_idx(1)][end_idx(2)];

    std::priority_queue<GridNodePtr, std::vector<GridNodePtr>, NodeComparator> empty;
    openSet_.swap(empty);

    GridNodePtr neighborPtr = NULL;
    GridNodePtr current = NULL;

    startPtr->index = start_idx;
    startPtr->rounds = rounds_;
    startPtr->gScore = 0;
    startPtr->fScore = getHeu(startPtr, endPtr);
    startPtr->state = GridNode::OPENSET; //put start node in open set
    startPtr->cameFrom = NULL;
    openSet_.push(startPtr); //put start in open set

    endPtr->index = end_idx;

    double tentative_gScore;

    int num_iter = 0;
    while (!openSet_.empty())
    {
        num_iter++;
        current = openSet_.top();
        openSet_.pop();

        // if ( num_iter < 10000 )
        //     cout << "current=" << current->index.transpose() << endl;

        if (current->index(0) == endPtr->index(0) && current->index(1) == endPtr->index(1) && current->index(2) == endPtr->index(2))
        {
            // ros::Time time_2 = ros::Time::now();
            // printf("\033[34mA star iter:%d, time:%.3f\033[0m\n",num_iter, (time_2 - time_1).toSec()*1000);
            // if((time_2 - time_1).toSec() > 0.1)
            //     ROS_WARN("Time consume in A star path finding is %f", (time_2 - time_1).toSec() );
            gridPath_ = retrievePath(current);
            return true;
        }
        current->state = GridNode::CLOSEDSET; //move current node from open set to closed set.

        for (int dx = -1; dx <= 1; dx++)
            for (int dy = -1; dy <= 1; dy++)
                for (int dz = -1; dz <= 1; dz++)
                {
                    if (dx == 0 && dy == 0 && dz == 0)
                        continue;

                    Vector3i neighborIdx;
                    neighborIdx(0) = (current->index)(0) + dx;
                    neighborIdx(1) = (current->index)(1) + dy;
                    neighborIdx(2) = (current->index)(2) + dz;

                    if (neighborIdx(0) < 1 || neighborIdx(0) >= POOL_SIZE_(0) - 1 || neighborIdx(1) < 1 || neighborIdx(1) >= POOL_SIZE_(1) - 1 || neighborIdx(2) < 1 || neighborIdx(2) >= POOL_SIZE_(2) - 1)
                    {
                        continue;
                    }

                    neighborPtr = GridNodeMap_[neighborIdx(0)][neighborIdx(1)][neighborIdx(2)];
                    neighborPtr->index = neighborIdx;

                    bool flag_explored = neighborPtr->rounds == rounds_;

                    if (flag_explored && neighborPtr->state == GridNode::CLOSEDSET)
                    {
                        continue; //in closed set.
                    }

                    neighborPtr->rounds = rounds_;

                    if (checkOccupancy(Index2Coord(neighborPtr->index)))
                    {
                        continue;
                    }

                    double static_cost = sqrt(dx * dx + dy * dy + dz * dz);
                    tentative_gScore = current->gScore + static_cost;

                    if (!flag_explored)
                    {
                        //discover a new node
                        neighborPtr->state = GridNode::OPENSET;
                        neighborPtr->cameFrom = current;
                        neighborPtr->gScore = tentative_gScore;
                        neighborPtr->fScore = tentative_gScore + getHeu(neighborPtr, endPtr);
                        openSet_.push(neighborPtr); //put neighbor in open set and record it.
                    }
                    else if (tentative_gScore < neighborPtr->gScore)
                    { //in open set and need update
                        neighborPtr->cameFrom = current;
                        neighborPtr->gScore = tentative_gScore;
                        neighborPtr->fScore = tentative_gScore + getHeu(neighborPtr, endPtr);
                    }
                }
        ros::Time time_2 = ros::Time::now();
        if ((time_2 - time_1).toSec() > 0.2)
        {
            ROS_WARN("Failed in A star path searching !!! 0.2 seconds time limit exceeded.");
            return false;
        }
    }

    ros::Time time_2 = ros::Time::now();

    if ((time_2 - time_1).toSec() > 0.1)
        ROS_WARN("Time consume in A star path finding is %.3fs, iter=%d", (time_2 - time_1).toSec(), num_iter);

    return false;
}*/

vector<Vector3d> AStar::getPath()
{
    vector<Vector3d> raw_path;
    for (auto ptr : gridPath_)
        raw_path.push_back(Index2Coord(ptr->index));
    reverse(raw_path.begin(), raw_path.end());

    if (raw_path.size() < 2) return raw_path;

    // =================================================================
    // 🚀 新增极简视线剪枝 (LoS Pruning) - 消除锯齿，拯救 Jerk！
    // =================================================================
    vector<Vector3d> pruned_path;
    pruned_path.push_back(raw_path.front());

    int current_idx = 0;
    while (current_idx < (int)raw_path.size() - 1)
    {
        int next_idx = raw_path.size() - 1;
        // 从终点往回找，看 current_idx 能不能直接直线看到 next_idx
        while (next_idx > current_idx + 1)
        {
            Eigen::Vector3i idx1, idx2;
            Coord2Index(raw_path[current_idx], idx1);
            Coord2Index(raw_path[next_idx], idx2);
            
            // 如果两点之间视线无遮挡，直接拉直线！
            if (lineOfSight(idx1, idx2)) {
                break;
            }
            next_idx--;
        }
        pruned_path.push_back(raw_path[next_idx]);
        current_idx = next_idx;
    }

    // =================================================================
    // 保持原有的稠密插值 (对剪枝后极其平滑的拉扯线进行 0.1m 插值)
    // =================================================================
    vector<Vector3d> dense_path;
    dense_path.push_back(pruned_path[0]); 

    for (size_t i = 1; i < pruned_path.size(); ++i)
    {
        Vector3d p1 = pruned_path[i - 1];
        Vector3d p2 = pruned_path[i];
        
        double dist = (p2 - p1).norm();
        int num_segments = std::max(1, (int)std::ceil(dist / step_size_)); 

        for (int j = 1; j <= num_segments; ++j)
        {
            dense_path.push_back(p1 + (p2 - p1) * ((double)j / num_segments));
        }
    }

    path_length_ = 0.0;
    if (raw_path.size() > 1) {
        for (size_t i = 0; i < raw_path.size() - 1; ++i) {
            path_length_ += (raw_path[i+1] - raw_path[i]).norm();
        }
    }

    return dense_path;
}