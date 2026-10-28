// RRT* planner, see rrt_star.h for the references.

#include "rrt_star.h"

#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

RrtStar::RrtStar()
    : xmin(0.), xmax(0.), ymin(0.), ymax(0.), n_nodes(0), n_obs(0), rng_state(1u) {}

// xorshift32: small, deterministic across compilers (std::rand is not)
double RrtStar::uniform(double a, double b) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    return a + (b - a) * (rng_state / 4294967296.0);
}

bool RrtStar::point_free(double x, double y) const {
    for (int k = 0; k < n_obs; k++) {
        double dx = x - ox[k], dy = y - oy[k];
        if (dx * dx + dy * dy <= orad[k] * orad[k]) return false;
    }
    return true;
}

// ObstacleFree in Alg. 2: segment vs every inflated circle
bool RrtStar::obstacle_free(double x1, double y1, double x2, double y2) const {
    double sx = x2 - x1, sy = y2 - y1;
    double len2 = sx * sx + sy * sy;
    for (int k = 0; k < n_obs; k++) {
        double t = 0.;
        if (len2 > 1e-12) {
            t = ((ox[k] - x1) * sx + (oy[k] - y1) * sy) / len2;
            if (t < 0.) t = 0.;
            if (t > 1.) t = 1.;
        }
        double dx = x1 + t * sx - ox[k], dy = y1 + t * sy - oy[k];
        if (dx * dx + dy * dy <= orad[k] * orad[k]) return false;
    }
    return true;
}

int RrtStar::nearest(double x, double y) const {
    int best = 0;
    double best_d2 = 1e300;
    for (int i = 0; i < n_nodes; i++) {
        double dx = nx[i] - x, dy = ny[i] - y;
        double d2 = dx * dx + dy * dy;
        if (d2 < best_d2) { best_d2 = d2; best = i; }
    }
    return best;
}

int RrtStar::near_vertices(double x, double y, double r) {
    int n = 0;
    double r2 = r * r;
    for (int i = 0; i < n_nodes; i++) {
        double dx = nx[i] - x, dy = ny[i] - y;
        if (dx * dx + dy * dy <= r2) near_idx[n++] = i;
    }
    return n;
}

void RrtStar::attach(int child, int par) {
    parent[child] = par;
    next_sibling[child] = first_child[par];
    first_child[par] = child;
}

void RrtStar::detach(int child) {
    int par = parent[child];
    if (par < 0) return;
    if (first_child[par] == child) {
        first_child[par] = next_sibling[child];
    } else {
        for (int c = first_child[par]; c >= 0; c = next_sibling[c]) {
            if (next_sibling[c] == child) { next_sibling[c] = next_sibling[child]; break; }
        }
    }
    parent[child] = -1;
    next_sibling[child] = -1;
}

int RrtStar::add_node(double x, double y, int par, double c) {
    int i = n_nodes++;
    nx[i] = x; ny[i] = y; cost[i] = c;
    parent[i] = -1; first_child[i] = -1; next_sibling[i] = -1;
    if (par >= 0) attach(i, par);
    return i;
}

// After rewiring, every descendant of `root` changes cost by the same delta.
// Iterative DFS on the child lists, using near_idx as the stack is not allowed (in use),
// so walk with parent pointers instead.
void RrtStar::propagate_cost(int root, double delta) {
    int c = first_child[root];
    while (c >= 0) {
        cost[c] += delta;
        if (first_child[c] >= 0) { c = first_child[c]; continue; }
        while (c >= 0 && c != root && next_sibling[c] < 0) c = parent[c];
        if (c < 0 || c == root) break;
        c = next_sibling[c];
    }
}

int RrtStar::plan(double x0, double y0, double xg, double yg,
                  const RrtObstacle *obs, int n_obs_in, const RrtParams &p, RrtPath *out) {
    if (out == 0) return RRTS_BAD_INPUT;
    out->n = 0;
    out->length = 0.;
    if (p.max_iter < 1 || p.max_iter > RRTS_MAX_NODES - 1 || p.step_max <= 0. || p.goal_tol <= 0. ||
        n_obs_in < 0 || n_obs_in > RRTS_MAX_OBST || (n_obs_in > 0 && obs == 0))
        return RRTS_BAD_INPUT;

    // Inflated obstacles (configuration space of a point robot) and sampling box
    n_obs = n_obs_in;
    xmin = fmin(x0, xg); xmax = fmax(x0, xg);
    ymin = fmin(y0, yg); ymax = fmax(y0, yg);
    for (int k = 0; k < n_obs; k++) {
        ox[k] = obs[k].x; oy[k] = obs[k].y;
        orad[k] = obs[k].r + p.clearance;
        xmin = fmin(xmin, ox[k] - orad[k]); xmax = fmax(xmax, ox[k] + orad[k]);
        ymin = fmin(ymin, oy[k] - orad[k]); ymax = fmax(ymax, oy[k] + orad[k]);
    }
    xmin -= p.bounds_pad; xmax += p.bounds_pad;
    ymin -= p.bounds_pad; ymax += p.bounds_pad;

    if (!point_free(x0, y0)) return RRTS_START_BLOCKED;
    if (!point_free(xg, yg)) return RRTS_GOAL_BLOCKED;

    // NearVertices radius r = min(gamma * (log n / n)^(1/d), S_max), d = 2. The slides only say "a ball
    // with radius r_b" (slide 4/7); the shrinking radius is NearVertices(G, z, |V|) of [2].
    // UNVERIFIED: default gamma = 2 (1 + 1/d)^(1/d) (mu(X)/zeta_d)^(1/d), mu(X) = sampling-box area
    // >= mu(X_free), zeta_2 = pi; check the bound stated in [2] (06 deck, References).
    // TODO(dev): tune gamma, or replace r with a fixed r_b.
    double gamma = p.gamma;
    if (gamma <= 0.) gamma = 2. * sqrt(1.5) * sqrt((xmax - xmin) * (ymax - ymin) / M_PI);

    rng_state = (p.seed != 0u) ? p.seed : 1u;
    n_nodes = 0;
    add_node(x0, y0, -1, 0.);

    for (int it = 0; it < p.max_iter && n_nodes < RRTS_MAX_NODES; it++) {
        // z_rand <- Sample(i); z_nearest <- Nearest(G, z_rand)
        double xr = uniform(xmin, xmax), yr = uniform(ymin, ymax);
        int i_nearest = nearest(xr, yr);

        // z_new <- Steer(z_nearest, z_rand): move at most S_max towards z_rand (slide 3/7)
        double dx = xr - nx[i_nearest], dy = yr - ny[i_nearest];
        double d = sqrt(dx * dx + dy * dy);
        if (d < 1e-9) continue;
        if (d > p.step_max) { dx *= p.step_max / d; dy *= p.step_max / d; d = p.step_max; }
        double xn = nx[i_nearest] + dx, yn = ny[i_nearest] + dy;

        if (!obstacle_free(nx[i_nearest], ny[i_nearest], xn, yn)) continue;

        // Z_nearby <- NearVertices(G, z_new, |V|)
        double nv = (double)n_nodes + 1.;
        double r = fmin(gamma * sqrt(log(nv) / nv), p.step_max);
        int n_near = near_vertices(xn, yn, r);

        // Choose parent z_min: min Cost(z_near) + J(edge) with a free edge (slides 4/7, 5/7)
        int i_min = i_nearest;
        double c_min = cost[i_nearest] + d;
        for (int k = 0; k < n_near; k++) {
            int j = near_idx[k];
            double ex = xn - nx[j], ey = yn - ny[j];
            double c = cost[j] + sqrt(ex * ex + ey * ey);
            if (c < c_min && obstacle_free(nx[j], ny[j], xn, yn)) { i_min = j; c_min = c; }
        }
        int i_new = add_node(xn, yn, i_min, c_min);

        // Rewire: z_new becomes the parent of z_near if that lowers its cost (slide 6/7)
        for (int k = 0; k < n_near; k++) {
            int j = near_idx[k];
            if (j == i_min) continue;
            double ex = xn - nx[j], ey = yn - ny[j];
            double c = c_min + sqrt(ex * ex + ey * ey);
            if (c < cost[j] && obstacle_free(xn, yn, nx[j], ny[j])) {
                double delta = c - cost[j];
                detach(j);
                attach(j, i_new);
                cost[j] = c;
                propagate_cost(j, delta);
            }
        }
    }

    // Goal connection: cheapest node within goal_tol with a free segment to the goal
    int i_best = -1;
    double c_best = 1e300;
    for (int i = 0; i < n_nodes; i++) {
        double ex = xg - nx[i], ey = yg - ny[i];
        double dg = sqrt(ex * ex + ey * ey);
        if (dg <= p.goal_tol && cost[i] + dg < c_best && obstacle_free(nx[i], ny[i], xg, yg)) {
            i_best = i;
            c_best = cost[i] + dg;
        }
    }
    if (i_best < 0) return RRTS_NO_PATH;

    // Backtrack goal -> start, then reverse
    int n = 0;
    out->x[n] = xg; out->y[n] = yg; n++;
    for (int i = i_best; i >= 0 && n < RRTS_MAX_PATH; i = parent[i]) {
        out->x[n] = nx[i]; out->y[n] = ny[i]; n++;
    }
    if (n >= RRTS_MAX_PATH && out->x[n - 1] != x0) return RRTS_NO_PATH;  // path longer than the buffer
    for (int a = 0, b = n - 1; a < b; a++, b--) {
        double tx = out->x[a], ty = out->y[a];
        out->x[a] = out->x[b]; out->y[a] = out->y[b];
        out->x[b] = tx; out->y[b] = ty;
    }
    out->n = n;
    out->length = c_best;
    return RRTS_OK;
}
