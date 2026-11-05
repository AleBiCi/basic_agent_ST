// Offline closed-loop check of the kinematic Stanley controller (no simulator).
// Plant: kinematic single track model of 05 p. 16 (front-wheel point, no slip), DT = 0.05 s, plus the
// simulator's first-order steering lag tau_H = 0.03 s (vehicle_double_track_model.py delta_Dot).
// Path: an RRT*-like polyline with 2 m edges and sharp vertices. Prints max/final |e| per case.

#include "lateral.h"
#include "vehicle_params.h"

#include <math.h>
#include <stdio.h>

static const double DT_TEST = 0.05;     // [s] agent cycle
static const double TAU_H   = 0.03;     // [s] SteeringSystemData.tau_H

static void make_path(RrtPath *p) {
    // straight 20 m, 45 deg left over 20 m, then straight along +x, 2 m spacing
    p->n = 0;
    double x = 0., y = 0.;
    const double head[3] = {0., M_PI / 4., 0.};
    const int    nseg[3] = {10, 10, 15};
    p->x[p->n] = x; p->y[p->n] = y; p->n++;
    for (int k = 0; k < 3; k++)
        for (int i = 0; i < nseg[k]; i++) {
            x += 2. * cos(head[k]); y += 2. * sin(head[k]);
            p->x[p->n] = x; p->y[p->n] = y; p->n++;
        }
    p->length = 70.;
}

static int run(double y0, double psi0, double u) {
    RrtPath path;
    make_path(&path);
    LatParams lp = lat_default_params();
    LatState st;
    lat_reset(&st);
    PoseEstimator pose;   // only x, y, psi are used
    pose.x = 0.; pose.y = y0; pose.psi = psi0;
    double delta = 0., e_max = 0.;
    const int steps = (int)(path.length / u / DT_TEST) + 40;
    const int sub = 50;
    for (int k = 0; k < steps; k++) {
        double d_req;
        lat_stanley_step(&st, &path, pose, lat_front_speed(u, delta), lp, &d_req);
        if (fabs(st.e) > e_max && k * DT_TEST * u > 5.) e_max = fabs(st.e);
        for (int j = 0; j < sub; j++) {                 // plant integration inside one cycle
            const double h = DT_TEST / sub;
            const double vf = u / cos(delta);
            // integrate the front point, then move the CoM back by Lf
            double xf = pose.x + VEH_LF * cos(pose.psi), yf = pose.y + VEH_LF * sin(pose.psi);
            xf += vf * cos(pose.psi + delta) * h;
            yf += vf * sin(pose.psi + delta) * h;
            pose.psi += vf / VEH_L * sin(delta) * h;
            pose.x = xf - VEH_LF * cos(pose.psi);
            pose.y = yf - VEH_LF * sin(pose.psi);
            delta += (d_req - delta) / TAU_H * h;
        }
    }
    printf("y0=%5.2f psi0=%5.2f u=%4.1f  max|e| after 5 m = %.3f m  final e = %.4f m  final dpsi = %.4f rad  seg=%d/%d\n",
           y0, psi0, u, e_max, st.e, st.dpsi, st.seg, path.n - 2);
    return (fabs(st.e) < 0.05 && fabs(st.dpsi) < 0.02) ? 0 : 1;
}

int main() {
    int fail = 0;
    fail += run(0.0, 0.0, 2.0);
    fail += run(1.5, 0.0, 2.0);      // start 1.5 m left of the path
    fail += run(-2.0, 0.3, 4.0);     // right of the path, heading away
    fail += run(0.5, -0.4, 6.0);
    printf(fail ? "FAIL (%d)\n" : "all converged\n", fail);
    return fail;
}
