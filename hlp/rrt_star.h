// RRT* path planner in the plane (point robot, circular obstacles inflated by a clearance).
//
// Algorithm: 06_Route_and_HighLevel_Planning, "RRT* algorithm 1/7 - 7/7" and
// "RRT* algorithm - pseudo code" (slide 26), i.e. Karaman & Frazzoli, CDC 2010 [2],
// Algorithm 1 (main loop) and Algorithm 2 (Extend: Nearest, Steer, NearVertices,
// choose parent, rewire). Cost = path length (slide "RRT* algorithm - example of use - 5/7").
// Steer is a straight segment of at most S_max (slide 3/7): the path is a polyline, its curvature
// is discontinuous at the vertices (slide 27 "pros and cons", cons 2); smoothing with clothoids
// (same deck, "Clothoid curve") is left to the path-following stage.
//
// Fixed-size storage: no heap allocation, the whole object can be `static`.

#ifndef RRT_STAR_H
#define RRT_STAR_H

#define RRTS_MAX_NODES 4000   // max tree size (N in Alg. 1 cannot exceed this)
#define RRTS_MAX_OBST  20     // = size of ObjID[] in input_data_str
#define RRTS_MAX_PATH  512    // max vertices of the returned path

// Return codes of RrtStar::plan
#define RRTS_OK             0
#define RRTS_NO_PATH       -1   // no node reached the goal region within max_iter
#define RRTS_START_BLOCKED -2   // start inside an inflated obstacle
#define RRTS_GOAL_BLOCKED  -3   // goal inside an inflated obstacle
#define RRTS_BAD_INPUT     -4   // invalid parameters

struct RrtObstacle {
    double x, y;    // centre [m]
    double r;       // footprint radius [m] (clearance is added by the planner)
};

struct RrtParams {
    int      max_iter;      // number of samples N (Alg. 1); <= RRTS_MAX_NODES - 1
    double   step_max;      // S_max [m], max edge length (slide "RRT* algorithm 3/7")
    double   gamma;         // NearVertices constant [m]; <= 0 -> computed from the sampling area
    double   goal_tol;      // [m] nodes within this distance may connect to the goal
    double   clearance;     // [m] added to every obstacle radius (vehicle radius + margin)
    double   bounds_pad;    // [m] sampling box = bounding box of start/goal/obstacles + this pad
    unsigned seed;          // PRNG seed (fixed seed -> repeatable plans)
};

struct RrtPath {
    int    n;                   // number of vertices, start first, goal last
    double x[RRTS_MAX_PATH];    // [m]
    double y[RRTS_MAX_PATH];    // [m]
    double length;              // [m]
};

class RrtStar {
public:
    RrtStar();

    // Plans from (x0, y0) to (xg, yg). Returns one of the RRTS_* codes; on RRTS_OK `out` holds the path.
    int plan(double x0, double y0, double xg, double yg,
             const RrtObstacle *obs, int n_obs, const RrtParams &p, RrtPath *out);

    // Tree access (for logging/plots), valid after plan()
    int    node_count() const { return n_nodes; }
    double node_x(int i) const { return nx[i]; }
    double node_y(int i) const { return ny[i]; }
    int    node_parent(int i) const { return parent[i]; }

    // Sampling box used by the last plan() call
    double xmin, xmax, ymin, ymax;

private:
    // Tree: parent pointer + first-child / next-sibling lists (for cost propagation after rewiring)
    int    n_nodes;
    double nx[RRTS_MAX_NODES], ny[RRTS_MAX_NODES], cost[RRTS_MAX_NODES];
    int    parent[RRTS_MAX_NODES], first_child[RRTS_MAX_NODES], next_sibling[RRTS_MAX_NODES];
    int    near_idx[RRTS_MAX_NODES];

    // Obstacles, already inflated by the clearance
    int    n_obs;
    double ox[RRTS_MAX_OBST], oy[RRTS_MAX_OBST], orad[RRTS_MAX_OBST];

    unsigned rng_state;

    double uniform(double a, double b);
    bool   point_free(double x, double y) const;
    bool   obstacle_free(double x1, double y1, double x2, double y2) const;
    int    nearest(double x, double y) const;
    int    near_vertices(double x, double y, double r);
    int    add_node(double x, double y, int par, double c);
    void   attach(int child, int par);
    void   detach(int child);
    void   propagate_cost(int root, double delta);
};

#endif
