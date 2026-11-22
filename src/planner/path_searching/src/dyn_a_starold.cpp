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

double AStar::computePathLength(const std::vector<Eigen::Vector3d> &path) const
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

double AStar::estimateClearanceMetricOnly(const Eigen::Vector3d &pt)
{
    if (grid_map_ == nullptr)
    {
        return 0.0;
    }
    if (!grid_map_->isInMap(pt))
    {
        return 0.0;
    }
    if (checkOccupancy(pt))
    {
        return 0.0;
    }

    const double res = grid_map_->getResolution();
    static const int dirs[6][3] = {
        {1, 0, 0}, {-1, 0, 0},
        {0, 1, 0}, {0, -1, 0},
        {0, 0, 1}, {0, 0, -1}};

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
            if (reject_unknown_for_metric_ && grid_map_->isUnknown(probe_pt))
            {
                return d;
            }
        }
    }

    return (max_probe_steps_ + 1) * res;
}

bool AStar::segmentFreeSampledMetricOnly(const Eigen::Vector3d &p1,
                                         const Eigen::Vector3d &p2,
                                         bool ignore_start_occupancy,
                                         bool ignore_end_occupancy)
{
    if (grid_map_ == nullptr)
    {
        return false;
    }

    Eigen::Vector3d diff = p2 - p1;
    const double dist = diff.norm();

    if (dist < 1e-6)
    {
        if (ignore_start_occupancy || ignore_end_occupancy)
        {
            return grid_map_->isInMap(p1);
        }
        return grid_map_->isInMap(p1) && !checkOccupancy(p1);
    }

    const double res = grid_map_->getResolution();
    const double sample_step = std::max(0.5 * res, 1e-3);
    const int steps = std::max(1, static_cast<int>(std::ceil(dist / sample_step)));

    for (int i = 0; i <= steps; ++i)
    {
        const double s = static_cast<double>(i) / static_cast<double>(steps);
        const bool is_start_sample = (i == 0);
        const bool is_end_sample = (i == steps);

        if ((is_start_sample && ignore_start_occupancy) ||
            (is_end_sample && ignore_end_occupancy))
        {
            continue;
        }

        Eigen::Vector3d p = p1 + s * diff;

        if (!grid_map_->isInMap(p))
        {
            return false;
        }
        if (checkOccupancy(p))
        {
            return false;
        }
        if (reject_unknown_for_metric_ && grid_map_->isUnknown(p))
        {
            return false;
        }
    }

    return true;
}

void AStar::updatePathMetricsFromPath(const std::vector<Eigen::Vector3d> &path)
{
    last_min_clearance_ = 0.0;
    path_length_ = 0.0;
    output_edge_fail_count_ = 0;

    if (path.empty())
    {
        segment_check_fail_count_ = edge_reject_count_;
        return;
    }

    path_length_ = computePathLength(path);

    double min_clearance = 1e9;
    for (const auto &p : path)
    {
        if (!p.allFinite())
        {
            min_clearance = 0.0;
            continue;
        }
        min_clearance = std::min(min_clearance, estimateClearanceMetricOnly(p));
    }

    if (min_clearance < 1e8)
    {
        last_min_clearance_ = min_clearance;
    }
    else
    {
        last_min_clearance_ = 0.0;
    }

    for (size_t i = 1; i < path.size(); ++i)
    {
        const bool ignore_start = (i == 1);
        const bool ignore_end = (i + 1 == path.size());

        if (!segmentFreeSampledMetricOnly(path[i - 1], path[i], ignore_start, ignore_end))
        {
            ++output_edge_fail_count_;
        }
    }

    segment_check_fail_count_ = edge_reject_count_ + output_edge_fail_count_;
}

void AStar::updatePathMetricsFromGridPath()
{
    std::vector<Eigen::Vector3d> path;
    path.reserve(gridPath_.size());

    for (auto ptr : gridPath_)
    {
        if (ptr == nullptr)
        {
            continue;
        }
        path.push_back(Index2Coord(ptr->index));
    }

    std::reverse(path.begin(), path.end());
    updatePathMetricsFromPath(path);
}

std::vector<Eigen::Vector3d> AStar::buildFinalPathFromGridPath() const
{
    std::vector<Eigen::Vector3d> path;

    if (!has_last_query_ || gridPath_.empty())
    {
        return path;
    }

    path.push_back(last_query_start_pt_);

    std::vector<Eigen::Vector3d> grid_centers;
    grid_centers.reserve(gridPath_.size());
    for (auto ptr : gridPath_)
    {
        if (ptr != nullptr)
        {
            grid_centers.push_back(Index2Coord(ptr->index));
        }
    }

    std::reverse(grid_centers.begin(), grid_centers.end());
    for (const auto &p : grid_centers)
    {
        if (path.empty() || (p - path.back()).norm() > 1e-4)
        {
            path.push_back(p);
        }
    }

    if (path.empty() || (last_query_end_pt_ - path.back()).norm() > 1e-4)
    {
        path.push_back(last_query_end_pt_);
    }

    return path;
}

bool AStar::validateFinalPathAndUpdateMetrics(const std::vector<Eigen::Vector3d> &path,
                                              std::string *reason)
{
    path_length_ = 0.0;
    last_min_clearance_ = 0.0;
    output_edge_fail_count_ = 0;
    segment_check_fail_count_ = edge_reject_count_;

    if (reason != nullptr)
    {
        *reason = "NONE";
    }

    if (grid_map_ == nullptr)
    {
        if (reason != nullptr)
        {
            *reason = "PATH_NO_MAP";
        }
        return false;
    }

    if (path.size() < 3)
    {
        if (reason != nullptr)
        {
            *reason = "PATH_TOO_SHORT";
        }
        return false;
    }

    double min_clearance = 1e9;

    for (size_t i = 0; i < path.size(); ++i)
    {
        if (!path[i].allFinite())
        {
            if (reason != nullptr)
            {
                *reason = "PATH_POINT_NAN";
            }
            return false;
        }

        if (!grid_map_->isInMap(path[i]))
        {
            if (reason != nullptr)
            {
                *reason = "PATH_POINT_OUT_OF_MAP";
            }
            return false;
        }

        const bool is_endpoint = (i == 0 || i + 1 == path.size());
        if (!is_endpoint && checkOccupancy(path[i]))
        {
            if (reason != nullptr)
            {
                *reason = "PATH_POINT_IN_COLLISION";
            }
            last_min_clearance_ = 0.0;
            return false;
        }

        if (!is_endpoint)
        {
            min_clearance = std::min(min_clearance, estimateClearanceMetricOnly(path[i]));
        }
    }

    for (size_t i = 1; i < path.size(); ++i)
    {
        const bool ignore_start = (i == 1);
        const bool ignore_end = (i + 1 == path.size());

        if (!segmentFreeSampledMetricOnly(path[i - 1], path[i], ignore_start, ignore_end))
        {
            ++output_edge_fail_count_;
            segment_check_fail_count_ = edge_reject_count_ + output_edge_fail_count_;
            if (reason != nullptr)
            {
                *reason = "PATH_EDGE_COLLISION";
            }
            last_min_clearance_ = 0.0;
            return false;
        }

        path_length_ += (path[i] - path[i - 1]).norm();
    }

    if (min_clearance < 1e8)
    {
        last_min_clearance_ = min_clearance;
    }
    else
    {
        last_min_clearance_ = 0.0;
    }

    segment_check_fail_count_ = edge_reject_count_ + output_edge_fail_count_;
    return true;
}

bool AStar::ConvertToIndexAndAdjustStartEndPoints(Vector3d start_pt,
                                                  Vector3d end_pt,
                                                  Vector3i &start_idx,
                                                  Vector3i &end_idx,
                                                  Vector3d *adjusted_start_pt,
                                                  Vector3d *adjusted_end_pt)
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

    if (adjusted_start_pt != nullptr)
    {
        *adjusted_start_pt = start_pt;
    }
    if (adjusted_end_pt != nullptr)
    {
        *adjusted_end_pt = end_pt;
    }

    return true;
}

bool AStar::AstarSearch(const double step_size, Vector3d start_pt, Vector3d end_pt)
{
    ros::Time time_1 = ros::Time::now();

    expanded_nodes_ = 0;
    path_length_ = 0.0;
    last_search_time_ms_ = 0.0;
    last_min_clearance_ = 0.0;
    segment_check_fail_count_ = 0;
    edge_reject_count_ = 0;
    output_edge_fail_count_ = 0;
    used_fallback_ = false;
    last_fallback_reason_ = "NONE";
    gridPath_.clear();
    final_path_.clear();
    last_query_start_pt_ = start_pt;
    last_query_end_pt_ = end_pt;
    last_adjusted_start_pt_ = start_pt;
    last_adjusted_end_pt_ = end_pt;
    has_last_query_ = true;
    last_path_invalid_reason_ = "NONE";

    auto finishSearch = [&](bool ok) -> bool
    {
        last_search_time_ms_ = (ros::Time::now() - time_1).toSec() * 1000.0;
        return ok;
    };

    ++rounds_;

    step_size_ = step_size;
    inv_step_size_ = 1 / step_size;
    center_ = (start_pt + end_pt) / 2;

    Vector3i start_idx, end_idx;
    if (!ConvertToIndexAndAdjustStartEndPoints(start_pt, end_pt,
                                               start_idx, end_idx,
                                               &last_adjusted_start_pt_,
                                               &last_adjusted_end_pt_))
    {
        last_fallback_reason_ = "ASTAR_START_END_INVALID";
        ROS_ERROR("Unable to handle the initial or end point, force return!");
        return finishSearch(false);
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
        expanded_nodes_ = num_iter;

        // if ( num_iter < 10000 )
        //     cout << "current=" << current->index.transpose() << endl;

        if (current->index(0) == endPtr->index(0) && current->index(1) == endPtr->index(1) && current->index(2) == endPtr->index(2))
        {
            // ros::Time time_2 = ros::Time::now();
            // printf("\033[34mA star iter:%d, time:%.3f\033[0m\n",num_iter, (time_2 - time_1).toSec()*1000);
            // if((time_2 - time_1).toSec() > 0.1)
            //     ROS_WARN("Time consume in A star path finding is %f", (time_2 - time_1).toSec() );
            gridPath_ = retrievePath(current);
            std::vector<Eigen::Vector3d> candidate_path = buildFinalPathFromGridPath();
            std::string invalid_reason = "NONE";
            bool candidate_valid = false;

            if (candidate_path.size() < 3 && has_last_query_)
            {
                std::vector<Eigen::Vector3d> short_path;
                short_path.push_back(last_query_start_pt_);
                short_path.push_back(0.5 * (last_query_start_pt_ + last_query_end_pt_));
                short_path.push_back(last_query_end_pt_);

                candidate_valid = validateFinalPathAndUpdateMetrics(short_path, &invalid_reason);
                if (candidate_valid)
                {
                    candidate_path.swap(short_path);
                }
            }
            else
            {
                candidate_valid = validateFinalPathAndUpdateMetrics(candidate_path, &invalid_reason);
            }

            if (!candidate_valid)
            {
                final_path_.clear();
                last_path_invalid_reason_ = invalid_reason;
                last_fallback_reason_ = invalid_reason;

                ROS_WARN("[AStarBaseline] search found grid path but final guide invalid. reason=%s, grid_nodes=%zu, start=(%.2f %.2f %.2f), end=(%.2f %.2f %.2f)",
                         invalid_reason.c_str(), gridPath_.size(),
                         last_query_start_pt_.x(), last_query_start_pt_.y(), last_query_start_pt_.z(),
                         last_query_end_pt_.x(), last_query_end_pt_.y(), last_query_end_pt_.z());

                return finishSearch(false);
            }

            final_path_ = candidate_path;
            last_path_invalid_reason_ = "NONE";
            last_fallback_reason_ = "NONE";

            ROS_WARN_THROTTLE(0.5,
                              "[AStarBaseline] valid path size=%zu, grid_nodes=%zu, length=%.3f, min_clearance=%.3f, seg_reject=%d",
                              final_path_.size(), gridPath_.size(), path_length_, last_min_clearance_, segment_check_fail_count_);

            return finishSearch(true);
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

                    const Eigen::Vector3d current_coord = Index2Coord(current->index);
                    const Eigen::Vector3d neighbor_coord = Index2Coord(neighborPtr->index);

                    if (!segmentFreeSampledMetricOnly(current_coord, neighbor_coord))
                    {
                        ++edge_reject_count_;
                        segment_check_fail_count_ = edge_reject_count_ + output_edge_fail_count_;
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
            last_fallback_reason_ = "ASTAR_TIMEOUT_0P2S";
            ROS_WARN("Failed in A star path searching !!! 0.2 seconds time limit exceeded.");
            return finishSearch(false);
        }
    }

    ros::Time time_2 = ros::Time::now();

    if ((time_2 - time_1).toSec() > 0.1)
        ROS_WARN("Time consume in A star path finding is %.3fs, iter=%d", (time_2 - time_1).toSec(), num_iter);

    last_fallback_reason_ = "ASTAR_OPENSET_EMPTY";
    return finishSearch(false);
}

vector<Vector3d> AStar::getPath()
{
    if (!final_path_.empty())
    {
        return final_path_;
    }

    ROS_WARN_THROTTLE(0.5,
                      "[AStarBaseline] getPath called but final_path_ is empty. last_reason=%s",
                      last_path_invalid_reason_.c_str());
    return vector<Vector3d>();
}
