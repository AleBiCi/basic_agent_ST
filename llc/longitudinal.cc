// Longitudinal control, see longitudinal.h.

#include "longitudinal.h"
#include "hlp.h"            // OBJ_ID_TARGET
#include "primitives.h"     // MATLAB Coder: student_*_primitive, coeffs_j_opt

#include <math.h>

LgtParams lgt_default_params(const input_data_str *in) {
    LgtParams p;
    // TODO(dev): off-road cruising speed. The simulator sends RequestedCruisingSpeed = 20 m/s
    // (pydrivingsim/agent.py), probably too fast off-road; cap it or derive it from path curvature.
    p.vr             = in->RequestedCruisingSpeed;
    p.lookahead_min  = 50.0;    // slide 104 line 4
    p.lookahead_time = 5.0;     // slide 104 line 4
    p.stop_margin    = 0.0;     // TODO(dev): where to stop with respect to the target point
    p.a_min          = 5.0;     // TODO(dev): saturation limits (slide 68), from the vehicle/tyres
    p.a_max          = 3.0;     // TODO(dev)
    p.kp             = 0.2;     // TODO(dev): PI gains (slide 69), values from the previous code
    p.ki             = 0.1;     // TODO(dev)
    p.ped_min        = -1.0;    // the simulator does not clamp the pedal (vehicle_double_track_model.py)
    p.ped_max        = 1.0;
    p.v_stop         = 0.1;     // TODO(dev)
    p.ped_hold       = -0.3;    // TODO(dev): brake pedal to hold the vehicle once stopped
    return p;
}

double lgt_path_remaining(const RrtPath *path, const PoseEstimator &pose, int idx) {
    if (path->n <= 0) return -1.;
    if (idx < 0) idx = 0;
    if (idx >= path->n) idx = path->n - 1;
    // From the vehicle to the vertex after the closest one, then along the path to the goal
    int k = (idx + 1 < path->n) ? idx + 1 : idx;
    double s = hypot(path->x[k] - pose.x, path->y[k] - pose.y);
    for (int i = k; i + 1 < path->n; i++)
        s += hypot(path->x[i + 1] - path->x[i], path->y[i + 1] - path->y[i]);
    return s;
}

double lgt_target_distance(const input_data_str *in) {
    for (int i = 0; i < in->NrObjs && i < 20; i++)
        if (in->ObjID[i] == OBJ_ID_TARGET) return hypot(in->ObjX[i], in->ObjY[i]);
    return -1.;
}

static void zero6(double m[6]) { for (int i = 0; i < 6; i++) m[i] = 0.; }

int lgt_select(double v0, double a0, double s_target, const LgtParams &p, double m[6], double *T) {
    zero6(m);
    *T = 0.;
    double lookahead = fmax(p.lookahead_min, p.lookahead_time * v0);   // slide 104 line 4

    // Slide 104 line 16 / 105 line 17: free flow when there is no stop point within the lookahead.
    // TODO(dev): off-road policy on top of this (speed limit from curvature, obstacles close to the
    // path, mirroring between FreeFlow and Stop as in slide 106 lines 37-44).
    if (s_target < 0. || s_target >= lookahead) {
        double m2[6], v1, T1, v2, T2;
        // FreeFlow(v0, a0, xf, vr) = Pass(v0, a0, xf, vr, vr, 0, 0); output order of the .m:
        // [coeffsT2, v2, T2, coeffsT1, v1, T1], with vr == vr both primitives coincide.
        student_pass_primitive(v0, a0, lookahead, p.vr, p.vr, 0., 0., m2, &v2, &T2, m, &v1, &T1);
        if (T1 > 0.) { *T = T1; return LGT_FREEFLOW; }
        zero6(m);
        // No free-flow solution (a0 < 0 with v* > vr, slide 58: "there is no solution, we need to
        // stop"). TODO(dev): Stop or PassJ0 here; for now no primitive -> constant velocity.
        return LGT_NONE;
    }

    // Target within the lookahead: the target plays the role of a red light (slide 106 line 35).
    double s_stop = s_target - p.stop_margin;
    double smax, tf;
    if (s_stop <= 0.) return LGT_NONE;
    student_stop_primitive(v0, a0, s_stop, m, &smax, &tf);
    if (tf > 0.) { *T = tf; return LGT_STOP; }
    zero6(m);
    return LGT_NONE;
}

void lgt_llc_reset(LgtLlc *s) {
    s->init  = false;
    s->a_int = 0.;
    s->e_int = 0.;
    s->a_req = 0.;
}

double lgt_llc_step(LgtLlc *s, const double m[6], double v_meas, double a_meas, double dt,
                    const LgtParams &p) {
    if (!s->init) {                       // slide 68: initialisation a0 = a0 (measured)
        s->a_int = a_meas;
        s->e_int = 0.;
        s->init  = true;
    }

    if (m == nullptr) {
        if (v_meas < p.v_stop) {          // stopped (target reached): hold with the brake
            s->a_int = 0.;
            s->e_int = 0.;
            s->a_req = 0.;
            return p.ped_hold;
        }
        // No primitive while moving: constant velocity (slide 70, "Test Stopping Primitive" line 7)
        s->a_int = 0.;
        s->a_req = 0.;
    } else {
        // Integrated jerk, trapezoidal, internal a0 (slide 68)
        double a_req = s->a_int + 0.5 * dt * (coeffs_j_opt(0., m) + coeffs_j_opt(dt, m));
        a_req = fmin(fmax(a_req, -p.a_min), p.a_max);    // saturation (slide 68)
        s->a_int = a_req;
        s->a_req = a_req;
    }

    // PI on the pedal (slide 69)
    double e = s->a_req - a_meas;
    double e_int_new = s->e_int + e * dt;
    double ped = p.kp * e + p.ki * e_int_new;
    // Pedal clamp + conditional integration (anti-windup). NOT IN SLIDES (slide 69 has no clamp);
    // needed because the simulator applies any pedal value.
    if (ped > p.ped_max)      ped = p.ped_max;
    else if (ped < p.ped_min) ped = p.ped_min;
    else                      s->e_int = e_int_new;
    return ped;
}
