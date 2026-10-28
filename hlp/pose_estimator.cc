// Pose estimator, see pose_estimator.h.

#include "pose_estimator.h"

#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double wrap_pi(double a) {
    while (a > M_PI) a -= 2. * M_PI;
    while (a <= -M_PI) a += 2. * M_PI;
    return a;
}

PoseEstimator::PoseEstimator()
    : x(0.), y(0.), psi(0.), landmarks_used(0), init_done(false), lane_heading0(0.), n_lm(0) {}

void PoseEstimator::init(const input_data_str *in) {
    x = 0.; y = 0.; psi = 0.;
    lane_heading0 = in->LaneHeading;
    n_lm = in->NrObjs;
    if (n_lm < 0) n_lm = 0;
    if (n_lm > POSE_MAX_LANDMARKS) n_lm = POSE_MAX_LANDMARKS;
    for (int i = 0; i < n_lm; i++) {
        lm_id[i] = in->ObjID[i];
        lm_x[i] = in->ObjX[i];      // W0 == vehicle frame at this cycle
        lm_y[i] = in->ObjY[i];
    }
    landmarks_used = n_lm;
    init_done = true;
}

void PoseEstimator::update(const input_data_str *in, double dt) {
    psi = wrap_pi(lane_heading0 - in->LaneHeading);
    double c = cos(psi), s = sin(psi);

    int n = in->NrObjs < n_lm ? in->NrObjs : n_lm;
    double sx = 0., sy = 0.;
    int used = 0;
    for (int i = 0; i < n; i++) {
        if (in->ObjID[i] != lm_id[i] || in->ObjVel[i] != 0.) continue;   // only static, same object
        sx += lm_x[i] - (c * in->ObjX[i] - s * in->ObjY[i]);
        sy += lm_y[i] - (s * in->ObjX[i] + c * in->ObjY[i]);
        used++;
    }
    if (used > 0) {
        x = sx / used;
        y = sy / used;
    } else {
        x += in->VLgtFild * c * dt;
        y += in->VLgtFild * s * dt;
    }
    landmarks_used = used;
}

void PoseEstimator::world_to_vehicle(double xw, double yw, double *xv, double *yv) const {
    double dx = xw - x, dy = yw - y, c = cos(psi), s = sin(psi);
    *xv = c * dx + s * dy;
    *yv = -s * dx + c * dy;
}

void PoseEstimator::vehicle_to_world(double xv, double yv, double *xw, double *yw) const {
    double c = cos(psi), s = sin(psi);
    *xw = x + c * xv - s * yv;
    *yw = y + s * xv + c * yv;
}
