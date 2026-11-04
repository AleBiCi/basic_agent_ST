// vehicle_params.h
// Vehicle constants the agent needs (lateral control, saturations). VEHICLE_SET must match the param_set the
// simulator runs (my_PyDrivingSim/vehicle_model/param_sets.py): 0 = "chimera", 1 = "polaris_utv", 2 = "cherokee_37in",
// 3 = "uaz469". Sets 1-3 are NOT IN SLIDES: Project Chrono data (BSD-3), see rock_crawler_params/README.md.
#ifndef VEHICLE_PARAMS_H
#define VEHICLE_PARAMS_H

#define VEHICLE_SET 0   // TODO(dev): pick the set, same as the simulator

#if VEHICLE_SET == 1        // "polaris_utv"
static const double VEH_LF         = 1.5515;        // [m] CoM to front axle
static const double VEH_LR         = 1.1638;        // [m] CoM to rear axle
static const double VEH_MASS       = 1360.1;        // [kg]
static const double VEH_MAX_WHL_AG = 0.4886921906;  // [rad] = 28 deg, Polaris_PitmanArm.json "Maximum Angle (deg)"
static const double VEH_VMAX       = 6.3;           // [m/s] top speed in the fixed low gear (open-loop sim run)
#elif VEHICLE_SET == 2      // "cherokee_37in"
static const double VEH_LF         = 1.1338;        // [m]
static const double VEH_LR         = 1.4442;        // [m]
static const double VEH_MASS       = 2155.4;        // [kg]
static const double VEH_MAX_WHL_AG = 0.4886921906;  // [rad] = 28 deg, Cherokee_Steering.cpp m_maxAngle
static const double VEH_VMAX       = 4.4;           // [m/s] top speed in the fixed low gear (open-loop sim run)
#elif VEHICLE_SET == 3      // "uaz469"
static const double VEH_LF         = 1.1888;        // [m]
static const double VEH_LR         = 1.1912;        // [m]
static const double VEH_MASS       = 2785.0;        // [kg]
static const double VEH_MAX_WHL_AG = 0.4712388980;  // [rad] = 27 deg, UAZ469_Vehicle.json "Maximum Steering Angle (deg)"
static const double VEH_VMAX       = 4.3;           // [m/s] top speed in the fixed 1st gear (open-loop sim run)
#else                       // "chimera" (simulator default)
static const double VEH_LF         = 0.904;         // [m] vehicle_params.py VehicleData
static const double VEH_LR         = 0.636;         // [m]
static const double VEH_MASS       = 315.0;         // [kg]
static const double VEH_MAX_WHL_AG = 0.0;           // TODO(dev): no steering limit in the Chimera set
static const double VEH_VMAX       = 30.0;          // [m/s] open-loop full-throttle run
#endif
static const double VEH_L          = VEH_LF + VEH_LR;   // [m] wheelbase
// The simulator treats RequestedSteerWhlAg as a steering-wheel angle: steady-state front-wheel angle =
// RequestedSteerWhlAg / tau_D (vehicle_double_track_model.py, delta_Dot), while SteerWhlAg feedback is the
// front-wheel angle (agent.py reads state[10]). tau_D = 10 in all sets.
static const double VEH_STEER_RATIO = 10.0;            // [-] SteeringSystemData.tau_D

#endif // VEHICLE_PARAMS_H
