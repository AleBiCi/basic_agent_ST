// Longitudinal control: action priming, action selection and low-level control.
//
// project_info deck:
//   - motor primitives: Stop (slides 47-50), Pass / FreeFlow (slides 52-60, 104 line 3)
//   - action selection: main logic slides 104-106, with the RRT* target used as the stop point
//     in place of the traffic-light stop line (local adaptation, see lgt_select)
//   - low-level control: integrated jerk, trapezoidal, internal a0, saturation (slide 68),
//     PI on the pedal (slide 69)
//
// Coefficient vector convention (coef_list_fun): m = [c0 = 0, c1 = v0, c2 = a0, c3 = j0, c4, c5],
// so m[3] is j(0) (slide 76).

#ifndef LONGITUDINAL_H
#define LONGITUDINAL_H

#include "interfaces_data_structs.h"
#include "rrt_star.h"
#include "pose_estimator.h"

#define LGT_NONE     0   // no valid primitive (e.g. standing still at the target)
#define LGT_FREEFLOW 1   // Pass(v0, a0, lookahead, vr, vr, 0, 0), slide 104 line 3
#define LGT_STOP     2   // Stop(v0, a0, s_stop), slide 50

struct LgtParams {
    double vr;              // [m/s]  cruising speed for FreeFlow
    double lookahead_min;   // [m]    slide 104 line 4: lookahead = max(lookahead_min, lookahead_time * v0)
    double lookahead_time;  // [s]
    double stop_margin;     // [m]    stop this far before the target point
    double a_min;           // [m/s^2] > 0, max deceleration (slide 68 saturation: -a_min)
    double a_max;           // [m/s^2] > 0, max acceleration
    double kp, ki;          // [-]    PI gains on the acceleration error (slide 69)
    double ped_min, ped_max;// [-]    pedal range accepted by the simulator
    double v_stop;          // [m/s]  below this, with no primitive, the vehicle is considered stopped
    double ped_hold;        // [-]    pedal while stopped (negative = brake)
};

LgtParams lgt_default_params(const input_data_str *in);

// Remaining length of the RRT* path from the current pose to the goal [m].
// idx = vertex index kept by hlp_path_to_trajectory.
double lgt_path_remaining(const RrtPath *path, const PoseEstimator &pose, int idx);

// Straight-line distance to the first ObjID 3 in the message [m], -1 if there is none.
double lgt_target_distance(const input_data_str *in);

// Action priming + selection. s_target < 0 means "no stop point" (pure FreeFlow).
// Writes the chosen primitive in m[6] and its duration in *T. Returns LGT_*.
int lgt_select(double v0, double a0, double s_target, const LgtParams &p, double m[6], double *T);

// Low-level controller state, kept across cycles (static in the caller).
struct LgtLlc {
    bool   init;
    double a_int;   // internal a0 (slide 67/68)
    double e_int;   // integral of the acceleration error (slide 69)
    double a_req;   // last requested acceleration (for logging)
};

void lgt_llc_reset(LgtLlc *s);

// One LLC step. m = chosen primitive, or nullptr if lgt_select returned LGT_NONE.
// Returns the requested pedal (RequestedAcc in the manoeuvre message).
double lgt_llc_step(LgtLlc *s, const double m[6], double v_meas, double a_meas, double dt,
                    const LgtParams &p);

#endif
