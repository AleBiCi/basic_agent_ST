// Vehicle pose in the planning frame W0 for the high-level planner.
//
// The scenario message has no absolute pose. W0 is the vehicle frame at the first received
// cycle (origin on the vehicle position, x forward, y left), so at that cycle the objects'
// ObjX/ObjY are already their W0 coordinates. Later cycles:
//   - heading psi: LaneHeading0 - LaneHeading (the simulator sends LaneHeading = -yaw, see
//     my_PyDrivingSim/pydrivingsim/agent.py), wrapped to (-pi, pi];
//   - position: static objects are landmarks. With L_i its W0 position and r_i = (ObjX, ObjY) now,
//     p = L_i - R(psi) r_i; the estimate is the mean over the objects that kept index and ObjID;
//   - fallback with no usable landmark: dead reckoning x += VLgtFild cos(psi) DT,
//     y += VLgtFild sin(psi) DT (kinematic model without side slip).
// NOT IN SLIDES: landmark localisation. Closest: 04_Planning_and_Control_Overview, framework block
// "state estimation"; Perception Systems deck. The simulator objects are noise-free and always
// visible, so the estimate is exact; a real sensor would need filtering.

#ifndef POSE_ESTIMATOR_H
#define POSE_ESTIMATOR_H

#include "interfaces_data_structs.h"

#define POSE_MAX_LANDMARKS 20   // = size of ObjID[] in input_data_str

class PoseEstimator {
public:
    PoseEstimator();

    bool initialised() const { return init_done; }
    void init(const input_data_str *in);                 // first cycle: defines W0
    void update(const input_data_str *in, double dt);    // every later cycle

    // Point conversions between W0 and the current vehicle frame
    void world_to_vehicle(double xw, double yw, double *xv, double *yv) const;
    void vehicle_to_world(double xv, double yv, double *xw, double *yw) const;

    double x, y, psi;       // pose in W0 [m, m, rad]
    int landmarks_used;     // 0 -> last update was dead reckoning

private:
    bool   init_done;
    double lane_heading0;
    int    n_lm;
    int    lm_id[POSE_MAX_LANDMARKS];
    double lm_x[POSE_MAX_LANDMARKS], lm_y[POSE_MAX_LANDMARKS];
};

#endif
