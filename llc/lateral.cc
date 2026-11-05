// Lateral control, see lateral.h.

#include "lateral.h"
#include "vehicle_params.h"

#include <math.h>

static const double LAT_DELTA_MAX_FALLBACK = 0.5;   // [rad] TODO(dev): Chimera set has no steering limit

static double wrap_pi(double a) {
    while (a > M_PI)   a -= 2. * M_PI;
    while (a <= -M_PI) a += 2. * M_PI;
    return a;
}

LatParams lat_default_params() {
    LatParams p;
    p.k_e          = 1.0;   // TODO(dev): tune in simulation (p. 44 "k_e is a tunable coefficient")
    p.v_min        = 0.5;   // TODO(dev): guard only; NOT IN SLIDES as such (closest: p. 49, k_v term
                            // "to improve the low-speed behavior")
    p.delta_max    = VEH_MAX_WHL_AG > 0. ? VEH_MAX_WHL_AG : LAT_DELTA_MAX_FALLBACK;
    p.lf           = VEH_LF;
    p.search_ahead = 5;     // TODO(dev): >= max segments travelled in one cycle (RRT* S_max = 2 m)
    return p;
}

double lat_front_speed(double u, double delta_meas) {
    const double cd = cos(delta_meas);
    return (cd > 0.1) ? u / cd : u / 0.1;   // |delta| < 84 deg in practice; guard only
}

void lat_reset(LatState *s) {
    s->seg   = 0;
    s->e     = 0.;
    s->dpsi  = 0.;
    s->delta = 0.;
}

// Projection of (px, py) on segment i of the path. On the last segment t is not clamped above 1,
// so past the goal the segment is extended and psi_r stays defined.
static double project(const RrtPath *path, int i, double px, double py, double *t_out) {
    double ax = path->x[i], ay = path->y[i];
    double dx = path->x[i + 1] - ax, dy = path->y[i + 1] - ay;
    double l2 = dx * dx + dy * dy;
    double t = (l2 > 1e-12) ? ((px - ax) * dx + (py - ay) * dy) / l2 : 0.;
    if (t < 0.) t = 0.;
    if (t > 1. && i + 2 < path->n) t = 1.;
    *t_out = t;
    return hypot(ax + t * dx - px, ay + t * dy - py);
}

int lat_stanley_step(LatState *s, const RrtPath *path, const PoseEstimator &pose, double v_f,
                     const LatParams &p, double *delta) {
    *delta = 0.;
    if (path->n < 2) { s->delta = 0.; return LAT_NO_PATH; }
    if (s->seg < 0) s->seg = 0;
    if (s->seg > path->n - 2) s->seg = path->n - 2;

    // P_F, centre of the front wheel (p. 43)
    const double c = cos(pose.psi), sn = sin(pose.psi);
    const double pfx = pose.x + p.lf * c, pfy = pose.y + p.lf * sn;

    // P_P: closest point on the segments [seg, seg + search_ahead] (monotonic progress)
    int best = s->seg;
    double t_best = 0.;
    double d_best = project(path, best, pfx, pfy, &t_best);
    for (int i = s->seg + 1; i <= s->seg + p.search_ahead && i <= path->n - 2; i++) {
        double t;
        double d = project(path, i, pfx, pfy, &t);
        if (d < d_best) { d_best = d; best = i; t_best = t; }
    }
    s->seg = best;

    const double dx = path->x[best + 1] - path->x[best], dy = path->y[best + 1] - path->y[best];
    const double psi_r = atan2(dy, dx);                               // path heading at P_P
    const double ppx = path->x[best] + t_best * dx, ppy = path->y[best] + t_best * dy;

    // Signed error: lateral coordinate of P_P in the vehicle frame (> 0 = path on the left)
    const double e = -(ppx - pfx) * sn + (ppy - pfy) * c;
    const double dpsi = wrap_pi(psi_r - pose.psi);                    // p. 44: Dpsi = psi_r - psi

    // p. 44: delta = Dpsi + atan(k_e e / v_f), saturated
    const double vf = fmax(fabs(v_f), p.v_min);
    double d = dpsi + atan(p.k_e * e / vf);
    if (d >  p.delta_max) d =  p.delta_max;
    if (d < -p.delta_max) d = -p.delta_max;
    // TODO(dev): psi_r jumps at the RRT* vertices -> steering steps. Options: clothoid smoothing of the
    // path (06 "Clothoid curve", Bertolazzi & Frego) or the k_s damping term of p. 48.

    s->e = e;
    s->dpsi = dpsi;
    s->delta = d;
    *delta = d;
    return LAT_OK;
}
