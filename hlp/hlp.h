// High-level planner glue between the scenario message and RrtStar
// (04_Planning_and_Control_Overview, "General planning & control framework": HLP block).
//
// Simulator object IDs (my_PyDrivingSim/pydrivingsim/agent.py):
//   1 = cone (obstacle), 2 = coin (waypoint, ignored here), 3 = target, 4 = off-road obstacle.
// Obstacle footprint radius = max(ObjLen, ObjWidth) / 2.

#ifndef HLP_H
#define HLP_H

#include "interfaces_data_structs.h"
#include "rrt_star.h"
#include "pose_estimator.h"

#define OBJ_ID_CONE     1
#define OBJ_ID_TARGET   3
#define OBJ_ID_OBSTACLE 4

#define HLP_NOT_RUN   1     // planner not called yet
#define HLP_NO_TARGET 2     // no ObjID 3 in the message (e.g. the original road scenario)
// otherwise one of the RRTS_* codes of rrt_star.h

// Default RRT* parameters. TODO(dev): tune (iterations vs time, step, clearance margin).
RrtParams hlp_default_params(const input_data_str *in);

// Plans from the current pose to the target, obstacles taken from the message and moved to W0.
int hlp_plan(const input_data_str *in, const PoseEstimator &pose, const RrtParams &p,
             RrtStar &rrt, RrtPath *path);

// Writes the path ahead of the vehicle in TrajectoryPointIX/IY (vehicle frame), at most 23 points,
// starting from the path vertex closest to the vehicle. *idx keeps that vertex across cycles
// (monotonic, start with 0). Returns the number of points written.
int hlp_path_to_trajectory(const RrtPath *path, const PoseEstimator &pose, int *idx,
                           output_data_str *out);

#endif
