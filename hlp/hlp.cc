// HLP glue, see hlp.h.

#include "hlp.h"

#include <math.h>

#define TRAJ_MAX_POINTS 23      // = size of TrajectoryPointIX[] in output_data_str

RrtParams hlp_default_params(const input_data_str *in) {
    RrtParams p;
    p.max_iter   = 3000;
    p.step_max   = 2.0;                                                   // [m] S_max
    p.gamma      = 0.;                                                    // from the sampling area
    p.goal_tol   = 2.0;                                                   // [m]
    p.clearance  = 0.5 * sqrt(in->VehicleLen * in->VehicleLen +
                              in->VehicleWidth * in->VehicleWidth) + 0.3; // [m] circle around the
    // VehicleLen x VehicleWidth rectangle the simulator uses for collisions (scenarios/maps/README.md) + margin
    p.bounds_pad = 10.0;                                                  // [m]
    p.seed       = 12345u;
    return p;
}

int hlp_plan(const input_data_str *in, const PoseEstimator &pose, const RrtParams &p,
             RrtStar &rrt, RrtPath *path) {
    RrtObstacle obs[RRTS_MAX_OBST];
    int n_obs = 0;
    bool has_target = false;
    double xg = 0., yg = 0.;

    int n = in->NrObjs;
    if (n > RRTS_MAX_OBST) n = RRTS_MAX_OBST;
    for (int i = 0; i < n; i++) {
        double xw, yw;
        pose.vehicle_to_world(in->ObjX[i], in->ObjY[i], &xw, &yw);
        if (in->ObjID[i] == OBJ_ID_TARGET) {
            xg = xw; yg = yw; has_target = true;
        } else if (in->ObjID[i] == OBJ_ID_OBSTACLE || in->ObjID[i] == OBJ_ID_CONE) {
            obs[n_obs].x = xw;
            obs[n_obs].y = yw;
            obs[n_obs].r = 0.5 * fmax(in->ObjLen[i], in->ObjWidth[i]);
            n_obs++;
        }
    }
    if (!has_target) return HLP_NO_TARGET;
    return rrt.plan(pose.x, pose.y, xg, yg, obs, n_obs, p, path);
}

int hlp_path_to_trajectory(const RrtPath *path, const PoseEstimator &pose, int *idx,
                           output_data_str *out) {
    out->NTrajectoryPoints = 0;
    if (path == 0 || path->n <= 0) return 0;
    if (*idx < 0 || *idx >= path->n) *idx = 0;

    // advance while the next vertex is closer than the current one
    while (*idx + 1 < path->n) {
        double d0 = hypot(path->x[*idx] - pose.x, path->y[*idx] - pose.y);
        double d1 = hypot(path->x[*idx + 1] - pose.x, path->y[*idx + 1] - pose.y);
        if (d1 > d0) break;
        (*idx)++;
    }

    int k = 0;
    for (int i = *idx; i < path->n && k < TRAJ_MAX_POINTS; i++, k++) {
        double xv, yv;      // locals: no pointers into the packed struct
        pose.world_to_vehicle(path->x[i], path->y[i], &xv, &yv);
        out->TrajectoryPointIX[k] = xv;
        out->TrajectoryPointIY[k] = yv;
        out->TrajectoryPointITime[k] = 0.;      // geometric path, no timing yet
    }
    out->NTrajectoryPoints = k;
    return k;
}
