// Lateral control: kinematic Stanley controller on the RRT* path.
//
// 05_Low_Level_MotionPlanning deck (PDF page / printed slide number):
//   - "Stanley lateral controller - kinematic model - 1/2" (p. 43 / 37): P_F = centre of the front
//     wheel, P_P = closest path point to P_F, e = path tracking error
//   - "Stanley lateral controller - kinematic model - 2/2" (p. 44 / 38): Dpsi = psi_r - psi,
//     delta = Dpsi + atan(k_e e / v_f), delta saturated at delta_min / delta_max
//   - "Stanley lateral controller - discussion" (p. 50 / 44): the path must be dense or continuous
//     -> P_P is the projection on the polyline segments, not the closest vertex
// Reference cited in the slides: Hoffmann, Tomlin, Montemerlo, Thrun, ACC 2007.
//
// Frames and signs: W0 and pose as in pose_estimator.h (x forward, y left, psi CCW). The pose is the
// CoM (my_PyDrivingSim vehicle_double_track_model.py X[0], X[1]), so P_F = CoM + L_f (cos psi, sin psi).
// e is signed: e > 0 when the path is to the LEFT of P_F, so a positive (CCW, left) delta reduces it.
// With this sign the slide error dynamics de/dt = v_f sin(Dpsi - delta) hold as written.
// delta is the front-wheel angle; the simulator wants RequestedSteerWhlAg = VEH_STEER_RATIO * delta.

#ifndef LATERAL_H
#define LATERAL_H

#include "interfaces_data_structs.h"
#include "rrt_star.h"
#include "pose_estimator.h"

#define LAT_OK      0   // delta computed
#define LAT_NO_PATH 1   // path with fewer than 2 vertices: delta = 0

struct LatParams {
    double k_e;         // [1/s]  Stanley gain (p. 44), k_e e / v_f is dimensionless
    double v_min;       // [m/s]  lower bound on v_f in the atan term (guard against v_f -> 0)
    double delta_max;   // [rad]  > 0, front-wheel angle saturation (p. 44: delta_max, delta_min = -delta_max)
    double lf;          // [m]    CoM to front axle
    int    search_ahead;// [-]    segments searched ahead of the current one for the projection
};

LatParams lat_default_params();

// Controller state, kept across cycles (static in the caller).
struct LatState {
    int    seg;         // current path segment [seg, seg+1], monotonic (never moves back)
    double e;           // [m]   last signed tracking error (for logging)
    double dpsi;        // [rad] last heading error psi_r - psi, wrapped to (-pi, pi]
    double delta;       // [rad] last front-wheel angle after saturation
};

void lat_reset(LatState *s);

// Front-wheel speed from the CoM longitudinal speed u = VLgtFild and the measured front-wheel angle
// SteerWhlAg: no-slip kinematic model (p. 16, "Kinematic single track model - 2/2"), the rear wheel moves
// at v_f cos(delta), and u equals the rear-wheel speed when slip is neglected. cos(delta) is guarded.
double lat_front_speed(double u, double delta_meas);

// One Stanley step. v_f = front-wheel speed [m/s] from lat_front_speed.
// Writes the front-wheel angle in *delta. Returns LAT_*.
int lat_stanley_step(LatState *s, const RrtPath *path, const PoseEstimator &pose, double v_f,
                     const LatParams &p, double *delta);

#endif
