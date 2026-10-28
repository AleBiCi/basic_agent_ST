// Offline RRT* check, no simulator needed.
// stdin:  line 1 "x0 y0 xg yg", then one obstacle per line "x y d" (d = diameter) [m]
// stdout: CSV, "path,x,y" rows then "edge,x_parent,y_parent,x,y" rows (whole tree)
// Map file -> stdin, from my_PyDrivingSim (world frame; the plan does not depend on the start yaw):
//   python -c "import json,sys;m=json.load(open(sys.argv[1]));print(m['start']['x'],m['start']['y'],m['target']['x'],m['target']['y']);[print(o['x'],o['y'],o['d']) for o in m['obstacles']]" scenarios/maps/offroad_slalom.json | ./rrt_star_test

#include <stdio.h>
#include <math.h>
#include "rrt_star.h"

static RrtStar rrt;
static RrtPath path;

int main() {
    double x0, y0, xg, yg;
    if (scanf("%lf %lf %lf %lf", &x0, &y0, &xg, &yg) != 4) {
        fprintf(stderr, "expected 'x0 y0 xg yg' on the first line\n");
        return 1;
    }
    RrtObstacle obs[RRTS_MAX_OBST];
    int n_obs = 0;
    double ox, oy, od;
    while (n_obs < RRTS_MAX_OBST && scanf("%lf %lf %lf", &ox, &oy, &od) == 3) {
        obs[n_obs].x = ox; obs[n_obs].y = oy; obs[n_obs].r = od / 2.;
        n_obs++;
    }

    // Same values as hlp_default_params for the simulator vehicle (L = 1.54 m, Wf = 1.27 m)
    RrtParams p;
    p.max_iter = 3000; p.step_max = 2.0; p.gamma = 0.; p.goal_tol = 2.0;
    p.clearance = 0.5 * sqrt(1.54 * 1.54 + 1.27 * 1.27) + 0.3;
    p.bounds_pad = 10.0; p.seed = 12345u;

    int status = rrt.plan(x0, y0, xg, yg, obs, n_obs, p, &path);
    fprintf(stderr, "status %d, nodes %d, path vertices %d, length %.2f m, straight line %.2f m\n",
            status, rrt.node_count(), path.n, path.length, hypot(xg - x0, yg - y0));
    for (int i = 0; i < path.n; i++) printf("path,%f,%f\n", path.x[i], path.y[i]);
    for (int i = 1; i < rrt.node_count(); i++) {
        int j = rrt.node_parent(i);
        printf("edge,%f,%f,%f,%f\n", rrt.node_x(j), rrt.node_y(j), rrt.node_x(i), rrt.node_y(i));
    }
    return status == RRTS_OK ? 0 : 2;
}
