#include "path_searching/dyn_a_star_sca.h"

using namespace std;
using namespace Eigen;

AStarSCA::~AStarSCA()
{
    if (GridNodeMap_ == nullptr)
        return;

    for (int i = 0; i < POOL_SIZE_(0); i++)
        for (int j = 0; j < POOL_SIZE_(1); j++)
            for (int k = 0; k < POOL_SIZE_(2); k++)
                delete GridNodeMap_[i][j][k];
}

void AStarSCA::initGridMap(GridMap::Ptr occ_map, const Eigen::Vector3i pool_size)
{
    if (GridNodeMap_ != nullptr)
    {
        grid_map_ = occ_map;
        return;
    }

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

double AStarSCA::estimateClearance(const Eigen::Vector3d &pt)
{
    if (grid_map_ == nullptr) return 0.0;
    if (!grid_map_->isInMap(pt)) return 0.0;
    if (checkOccupancy(pt)) return 0.0;

    const double res = grid_map_->getResolution();
    static const int dirs[6][3] = {
        {1, 0, 0}, {-1, 0, 0},
        {0, 1, 0}, {0, -1, 0},
        {0, 0, 1}, {0, 0, -1}
    };

    for (int r = 1; r <= max_probe_steps_; ++r)
    {
        const double d = r * res;

        for (int i = 0; i < 6; ++i)
        {
            Eigen::Vector3d probe_pt(
                pt(0) + dirs[i][0] * d,
                pt(1) + dirs[i][1] * d,
                pt(2) + dirs[i][2] * d);

            if (!grid_map_->isInMap(probe_pt))
            {
                return d;
            }

            if (checkOccupancy(probe_pt))
            {
                return d;
            }

            if (reject_unknown_ && grid_map_->isUnknown(probe_pt))
            {
                return d;
            }
        }
    }

    return (max_probe_steps_ + 1) * res;
}

double AStarSCA::getClearanceCost(const Eigen::Vector3d &pt)
{
    const double res = grid_map_ ? grid_map_->getResolution() : step_size_;
    const double eps = std::max(0.5 * res, 1e-3);

    double clearance = estimateClearance(pt);
    if (clearance <= 1e-6)
    {
        return 1.0 / eps;
    }

    return 1.0 / std::max(clearance, eps);
}

double AStarSCA::getTurnPenalty(GridNodePtr current,
                                const Eigen::Vector3i &neighbor_idx) const
{
    if (current == nullptr || current->cameFrom == nullptr)
    {
        return 0.0;
    }

    Eigen::Vector3d a = (current->index - current->cameFrom->index).cast<double>();
    Eigen::Vector3d b = (neighbor_idx - current->index).cast<double>();

    double na = a.norm();
    double nb = b.norm();
    if (na < 1e-6 || nb < 1e-6)
    {
        return 0.0;
    }

    double dot = a.dot(b) / (na * nb);
    dot = std::max(-1.0, std::min(1.0, dot));

    return 1.0 - dot;
}

bool AStarSCA::segmentFreeSampled(const Eigen::Vector3d &p1,
                                  const Eigen::Vector3d &p2)
{
    if (grid_map_ == nullptr) return false;

    Eigen::Vector3d diff = p2 - p1;
    double dist = diff.norm();

    if (dist < 1e-6)
    {
        return grid_map_->isInMap(p1) && !checkOccupancy(p1);
    }

    const double res = grid_map_->getResolution();
    const double sample_step = std::max(0.5 * res, 1e-3);
    const int steps = std::max(1, static_cast<int>(std::ceil(dist / sample_step)));

    for (int i = 0; i <= steps; ++i)
    {
        double s = static_cast<double>(i) / static_cast<double>(steps);
        Eigen::Vector3d p = p1 + s * diff;

        if (!grid_map_->isInMap(p))
        {
            return false;
        }

        if (checkOccupancy(p))
        {
            return false;
        }

        if (reject_unknown_ && grid_map_->isUnknown(p))
        {
            return false;
        }
    }

    return true;
}

bool AStarSCA::edgeFeasible(const Eigen::Vector3i &idx1,
                            const Eigen::Vector3i &idx2,
                            bool conservative_check)
{
    if (grid_map_ == nullptr) return false;

    Eigen::Vector3d p1 = Index2Coord(idx1);
    Eigen::Vector3d p2 = Index2Coord(idx2);

    if (!grid_map_->isInMap(p1) || !grid_map_->isInMap(p2))
    {
        return false;
    }

    if (checkOccupancy(p1) || checkOccupancy(p2))
    {
        return false;
    }

    Eigen::Vector3d mid_pt = 0.5 * (p1 + p2);
    if (!grid_map_->isInMap(mid_pt))
    {
        return false;
    }

    if (checkOccupancy(mid_pt))
    {
        return false;
    }

    if (reject_unknown_ &&
        (grid_map_->isUnknown(p1) || grid_map_->isUnknown(mid_pt) || grid_map_->isUnknown(p2)))
    {
        return false;
    }

    if (conservative_check)
    {
        return segmentFreeSampled(p1, p2);
    }

    return true;
}

bool AStarSCA::lineOfSight(const Eigen::Vector3i &idx1,
                           const Eigen::Vector3i &idx2)
{
    return edgeFeasible(idx1, idx2, true);
}

bool AStarSCA::clearancePreservingLineOfSight(const Eigen::Vector3d &p1,
                                              const Eigen::Vector3d &p2,
                                              double min_required_clearance)
{
    if (!segmentFreeSampled(p1, p2))
    {
        return false;
    }

    if (min_required_clearance <= 1e-6 || grid_map_ == nullptr)
    {
        return true;
    }

    const Eigen::Vector3d diff = p2 - p1;
    const double dist = diff.norm();
    const double res = grid_map_->getResolution();
    const double sample_step = std::max(res, 1e-3);
    const int steps = std::max(1, static_cast<int>(std::ceil(dist / sample_step)));

    for (int i = 1; i < steps; ++i)
    {
        const double s = static_cast<double>(i) / static_cast<double>(steps);
        const Eigen::Vector3d p = p1 + s * diff;
        if (estimateClearance(p) + 1e-6 < min_required_clearance)
        {
            return false;
        }
    }

    return true;
}

double AStarSCA::getDiagHeu(GridNodePtr node1, GridNodePtr node2)
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

double AStarSCA::getManhHeu(GridNodePtr node1, GridNodePtr node2)
{
    double dx = abs(node1->index(0) - node2->index(0));
    double dy = abs(node1->index(1) - node2->index(1));
    double dz = abs(node1->index(2) - node2->index(2));

    return dx + dy + dz;
}

double AStarSCA::getEuclHeu(GridNodePtr node1, GridNodePtr node2)
{
    return (node2->index - node1->index).norm();
}

vector<GridNodePtr> AStarSCA::retrievePath(GridNodePtr current)
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

bool AStarSCA::validateGridPath(const std::vector<GridNodePtr> &path,
                                double *min_clearance)
{
    if (path.size() < 2)
    {
        last_fallback_reason_ = "PATH_TOO_SHORT";
        return false;
    }

    double min_c = 1e9;

    for (size_t i = 0; i < path.size(); ++i)
    {
        Eigen::Vector3d p = Index2Coord(path[i]->index);

        if (!grid_map_->isInMap(p) || checkOccupancy(p))
        {
            last_fallback_reason_ = "PATH_NODE_IN_COLLISION";
            return false;
        }

        min_c = std::min(min_c, estimateClearance(p));
    }

    for (size_t i = 1; i < path.size(); ++i)
    {
        if (!edgeFeasible(path[i - 1]->index, path[i]->index, true))
        {
            ++segment_check_fail_count_;
            last_fallback_reason_ = "PATH_EDGE_COLLISION";
            return false;
        }
    }

    if (min_clearance != nullptr)
    {
        *min_clearance = min_c;
    }

    if (min_path_clearance_limit_ > 1e-6 && min_c < min_path_clearance_limit_)
    {
        last_fallback_reason_ = "PATH_CLEARANCE_TOO_LOW";
        return false;
    }

    return true;
}

double AStarSCA::computePathLength(const std::vector<Eigen::Vector3d> &path) const
{
    if (path.size() < 2)
    {
        return 0.0;
    }

    double length = 0.0;
    for (size_t i = 1; i < path.size(); ++i)
    {
        length += (path[i] - path[i - 1]).norm();
    }

    return length;
}

bool AStarSCA::ConvertToIndexAndAdjustStartEndPoints(Vector3d start_pt,
                                                     Vector3d end_pt,
                                                     Vector3i &start_idx,
                                                     Vector3i &end_idx)
{
    if (!Coord2Index(start_pt, start_idx) || !Coord2Index(end_pt, end_idx))
    {
        return false;
    }

    const int max_adjust_steps = std::max(5, static_cast<int>(std::ceil(2.0 / std::max(step_size_, 1e-3))));

    if (checkOccupancy(Index2Coord(start_idx)))
    {
        Eigen::Vector3d dir = start_pt - end_pt;
        if (dir.norm() < 1e-6)
        {
            return false;
        }

        dir.normalize();

        int adjust_cnt = 0;
        do
        {
            start_pt = start_pt + dir * step_size_;

            if (!Coord2Index(start_pt, start_idx))
            {
                return false;
            }

            ++adjust_cnt;
            if (adjust_cnt > max_adjust_steps)
            {
                ROS_WARN("[A*] Failed to adjust start point out of obstacle.");
                return false;
            }
        } while (checkOccupancy(Index2Coord(start_idx)));
    }

    if (checkOccupancy(Index2Coord(end_idx)))
    {
        Eigen::Vector3d dir = end_pt - start_pt;
        if (dir.norm() < 1e-6)
        {
            return false;
        }

        dir.normalize();

        int adjust_cnt = 0;
        do
        {
            end_pt = end_pt + dir * step_size_;

            if (!Coord2Index(end_pt, end_idx))
            {
                return false;
            }

            ++adjust_cnt;
            if (adjust_cnt > max_adjust_steps)
            {
                ROS_WARN("[A*] Failed to adjust end point out of obstacle.");
                return false;
            }
        } while (checkOccupancy(Index2Coord(end_idx)));
    }

    return true;
}

bool AStarSCA::AstarSearch(const double step_size,
                           Vector3d start_pt,
                           Vector3d end_pt)
{
    const ros::Time search_start_time = ros::Time::now();

    auto finishSearch = [&](bool ok) -> bool
    {
        last_search_time_ms_ = (ros::Time::now() - search_start_time).toSec() * 1000.0;
        return ok;
    };

    expanded_nodes_ = 0;
    path_length_ = 0.0;
    segment_check_fail_count_ = 0;
    last_search_time_ms_ = 0.0;
    last_min_clearance_ = 0.0;
    used_fallback_ = false;
    last_fallback_reason_ = "NONE";
    gridPath_.clear();

    if (step_size <= 1e-6)
    {
        last_fallback_reason_ = "INVALID_STEP_SIZE";
        ROS_ERROR("[A*] Invalid step_size.");
        return finishSearch(false);
    }

    step_size_ = step_size;
    inv_step_size_ = 1.0 / step_size_;
    center_ = 0.5 * (start_pt + end_pt);

    bool sca_ok = searchOnce(start_pt, end_pt, true);
    if (sca_ok)
    {
        if (validateGridPath(gridPath_, &last_min_clearance_))
        {
            used_fallback_ = false;
            return finishSearch(true);
        }
    }

    std::string sca_fail_reason = last_fallback_reason_;
    ROS_WARN_THROTTLE(0.5, "[SCA-A*] Fallback to dense A*. reason=%s", sca_fail_reason.c_str());

    gridPath_.clear();

    bool dense_ok = searchOnce(start_pt, end_pt, false);
    if (dense_ok)
    {
        if (validateGridPath(gridPath_, &last_min_clearance_))
        {
            used_fallback_ = true;
            last_fallback_reason_ = sca_fail_reason;
            return finishSearch(true);
        }
    }

    ROS_WARN_THROTTLE(0.5, "[A*] Dense fallback also failed. reason=%s", last_fallback_reason_.c_str());

    gridPath_.clear();
    used_fallback_ = true;
    return finishSearch(false);
}

bool AStarSCA::searchOnce(const Eigen::Vector3d &start_pt,
                          const Eigen::Vector3d &end_pt,
                          bool use_sca)
{
    ros::Time time_1 = ros::Time::now();

    ++rounds_;

    Eigen::Vector3i start_idx, end_idx;
    if (!ConvertToIndexAndAdjustStartEndPoints(start_pt, end_pt, start_idx, end_idx))
    {
        last_fallback_reason_ = use_sca ? "SCA_START_END_INVALID" : "DENSE_START_END_INVALID";
        return false;
    }

    GridNodePtr startPtr = GridNodeMap_[start_idx(0)][start_idx(1)][start_idx(2)];
    GridNodePtr endPtr = GridNodeMap_[end_idx(0)][end_idx(1)][end_idx(2)];

    std::priority_queue<GridNodePtr, std::vector<GridNodePtr>, NodeComparator> empty;
    openSet_.swap(empty);

    startPtr->index = start_idx;
    startPtr->rounds = rounds_;
    startPtr->gScore = 0.0;
    startPtr->fScore = getHeu(startPtr, endPtr);
    startPtr->state = GridNode::OPENSET;
    startPtr->cameFrom = NULL;
    openSet_.push(startPtr);

    endPtr->index = end_idx;

    const int expand_step = use_sca ? expand_step_ : dense_expand_step_;
    const int max_expand_count = use_sca ? max_expand_count_sca_ : max_expand_count_fallback_;
    const double max_search_time = use_sca ? max_search_time_sca_ : max_search_time_fallback_;

    int local_expand_count = 0;

    while (!openSet_.empty())
    {
        ros::Time now = ros::Time::now();
        double elapsed = (now - time_1).toSec();

        if (elapsed > max_search_time)
        {
            last_search_time_ms_ = elapsed * 1000.0;
            last_fallback_reason_ = use_sca ? "SCA_TIMEOUT" : "DENSE_TIMEOUT";
            return false;
        }

        GridNodePtr current = openSet_.top();
        openSet_.pop();

        if (current->state == GridNode::CLOSEDSET)
        {
            continue;
        }

        ++expanded_nodes_;
        ++local_expand_count;

        if (local_expand_count > max_expand_count)
        {
            last_search_time_ms_ = elapsed * 1000.0;
            last_fallback_reason_ = use_sca ? "SCA_MAX_EXPAND" : "DENSE_MAX_EXPAND";
            return false;
        }

        double dx_goal = static_cast<double>(current->index(0) - endPtr->index(0));
        double dy_goal = static_cast<double>(current->index(1) - endPtr->index(1));
        double dz_goal = static_cast<double>(current->index(2) - endPtr->index(2));
        double dist_to_goal = std::sqrt(dx_goal * dx_goal + dy_goal * dy_goal + dz_goal * dz_goal);

        bool reach_goal = false;
        if (use_sca)
        {
            reach_goal = dist_to_goal <= expand_step * std::sqrt(3.0);
        }
        else
        {
            reach_goal = (current->index == endPtr->index);
        }

        if (reach_goal)
        {
            if (current->index == endPtr->index || edgeFeasible(current->index, endPtr->index, true))
            {
                if (current->index != endPtr->index)
                {
                    endPtr->rounds = rounds_;
                    endPtr->state = GridNode::OPENSET;
                    endPtr->cameFrom = current;
                    endPtr->gScore = current->gScore + (Index2Coord(endPtr->index) - Index2Coord(current->index)).norm();
                    endPtr->fScore = endPtr->gScore;
                    gridPath_ = retrievePath(endPtr);
                }
                else
                {
                    gridPath_ = retrievePath(current);
                }

                last_search_time_ms_ = (ros::Time::now() - time_1).toSec() * 1000.0;
                return true;
            }
        }

        current->state = GridNode::CLOSEDSET;

        double dynamic_epsilon = 1.0;
        if (use_sca)
        {
            dynamic_epsilon = 1.0 + (epsilon_ - 1.0) * (dist_to_goal / (POOL_SIZE_(0) / 2.0));
            dynamic_epsilon = std::max(1.0, std::min(dynamic_epsilon, epsilon_));
        }

        for (int dx_i = -expand_step; dx_i <= expand_step; dx_i += expand_step)
            for (int dy_i = -expand_step; dy_i <= expand_step; dy_i += expand_step)
                for (int dz_i = -expand_step; dz_i <= expand_step; dz_i += expand_step)
                {
                    if (dx_i == 0 && dy_i == 0 && dz_i == 0)
                    {
                        continue;
                    }

                    Eigen::Vector3i neighborIdx(
                        current->index(0) + dx_i,
                        current->index(1) + dy_i,
                        current->index(2) + dz_i);

                    int margin = std::max(1, expand_step);
                    if (neighborIdx(0) < margin || neighborIdx(0) >= POOL_SIZE_(0) - margin ||
                        neighborIdx(1) < margin || neighborIdx(1) >= POOL_SIZE_(1) - margin ||
                        neighborIdx(2) < margin || neighborIdx(2) >= POOL_SIZE_(2) - margin)
                    {
                        continue;
                    }

                    GridNodePtr neighborPtr = GridNodeMap_[neighborIdx(0)][neighborIdx(1)][neighborIdx(2)];
                    neighborPtr->index = neighborIdx;

                    bool flag_explored = neighborPtr->rounds == rounds_;
                    if (flag_explored && neighborPtr->state == GridNode::CLOSEDSET)
                    {
                        continue;
                    }

                    Eigen::Vector3d neighbor_pos = Index2Coord(neighborIdx);
                    if (!grid_map_->isInMap(neighbor_pos))
                    {
                        continue;
                    }

                    if (checkOccupancy(neighbor_pos))
                    {
                        continue;
                    }

                    if (reject_unknown_ && grid_map_->isUnknown(neighbor_pos))
                    {
                        continue;
                    }

                    if (!edgeFeasible(current->index, neighborIdx, conservative_edge_check_))
                    {
                        ++segment_check_fail_count_;
                        continue;
                    }

                    neighborPtr->rounds = rounds_;

                    const double move_cost = std::sqrt(
                        static_cast<double>(dx_i * dx_i + dy_i * dy_i + dz_i * dz_i)) * step_size_;

                    const double clearance_weight = use_sca ? weight_clearance_ : 0.5 * weight_clearance_;
                    const double turn_weight = use_sca ? weight_turn_ : 0.5 * weight_turn_;
                    const double clearance_cost = clearance_weight * getClearanceCost(neighbor_pos);
                    const double turn_cost = turn_weight * getTurnPenalty(current, neighborIdx);

                    double tentative_gScore = current->gScore + move_cost + clearance_cost + turn_cost;

                    if (!flag_explored || tentative_gScore < neighborPtr->gScore)
                    {
                        neighborPtr->state = GridNode::OPENSET;
                        neighborPtr->cameFrom = current;
                        neighborPtr->gScore = tentative_gScore;
                        neighborPtr->fScore = tentative_gScore + dynamic_epsilon * getHeu(neighborPtr, endPtr);
                        openSet_.push(neighborPtr);
                    }
                }
    }

    last_search_time_ms_ = (ros::Time::now() - time_1).toSec() * 1000.0;
    last_fallback_reason_ = use_sca ? "SCA_OPENSET_EMPTY" : "DENSE_OPENSET_EMPTY";
    return false;
}

vector<Vector3d> AStarSCA::getPath()
{
    vector<Vector3d> raw_path;

    for (auto ptr : gridPath_)
    {
        raw_path.push_back(Index2Coord(ptr->index));
    }

    reverse(raw_path.begin(), raw_path.end());

    if (raw_path.size() < 2)
    {
        path_length_ = computePathLength(raw_path);
        return raw_path;
    }

    double tmp_min_clearance = 0.0;
    if (!validateGridPath(gridPath_, &tmp_min_clearance))
    {
        path_length_ = 0.0;
        ROS_WARN_THROTTLE(0.5, "[A*] getPath rejected invalid gridPath. reason=%s", last_fallback_reason_.c_str());
        return vector<Vector3d>();
    }

    // Precompute axial clearance of raw_path nodes to preserve high-clearance corridors during LOS pruning
    vector<double> raw_clearances(raw_path.size(), 0.0);
    for (size_t k = 0; k < raw_path.size(); ++k)
    {
        raw_clearances[k] = estimateClearance(raw_path[k]);
    }

    vector<Vector3d> pruned_path;
    pruned_path.push_back(raw_path.front());

    int current_idx = 0;
    while (current_idx < static_cast<int>(raw_path.size()) - 1)
    {
        int next_idx = static_cast<int>(raw_path.size()) - 1;
        bool found = false;

        while (next_idx > current_idx + 1)
        {
            double shortcut_dist = (raw_path[next_idx] - raw_path[current_idx]).norm();
            if (shortcut_dist > los_prune_max_dist_)
            {
                --next_idx;
                continue;
            }

            double subpath_min_c = raw_clearances[current_idx];
            for (int k = current_idx + 1; k <= next_idx; ++k)
            {
                subpath_min_c = std::min(subpath_min_c, raw_clearances[k]);
            }

            // Reject shortcuts that cut into the obstacle-hugging 1-voxel layer (< 2*res),
            // while still allowing smooth pruning when clearance stays >= 2*res.
            const double res_prune = grid_map_ ? grid_map_->getResolution() : step_size_;
            const double min_safe_prune_c = std::min(subpath_min_c, 2.0 * res_prune);
            if (clearancePreservingLineOfSight(raw_path[current_idx], raw_path[next_idx], min_safe_prune_c))
            {
                found = true;
                break;
            }

            --next_idx;
        }

        if (!found)
        {
            next_idx = current_idx + 1;
        }

        pruned_path.push_back(raw_path[next_idx]);
        current_idx = next_idx;
    }

    vector<Vector3d> dense_path;
    dense_path.push_back(pruned_path.front());
    double final_min_clearance = estimateClearance(pruned_path.front());

    const double res = grid_map_ ? grid_map_->getResolution() : step_size_;
    const double dense_step = std::max(0.5 * res, dense_sample_resolution_ratio_ * step_size_);

    for (size_t i = 1; i < pruned_path.size(); ++i)
    {
        Vector3d p1 = pruned_path[i - 1];
        Vector3d p2 = pruned_path[i];

        double dist = (p2 - p1).norm();
        int num_segments = std::max(1, static_cast<int>(std::ceil(dist / dense_step)));

        for (int j = 1; j <= num_segments; ++j)
        {
            Vector3d p = p1 + (p2 - p1) * (static_cast<double>(j) / static_cast<double>(num_segments));

            if (!grid_map_->isInMap(p) || checkOccupancy(p))
            {
                ++segment_check_fail_count_;
                path_length_ = computePathLength(raw_path);
                ROS_WARN_THROTTLE(0.5, "[A*] Dense path interpolation hits obstacle. Return raw path instead.");
                return raw_path;
            }

            final_min_clearance = std::min(final_min_clearance, estimateClearance(p));
            dense_path.push_back(p);
        }
    }

    last_min_clearance_ = final_min_clearance;
    path_length_ = computePathLength(dense_path);
    return dense_path;
}
