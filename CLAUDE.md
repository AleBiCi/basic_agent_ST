# CLAUDE.md (basic_agent_ST)

C++ agent for PyDrivingSim. Workspace rules, references and modes are in `../../CLAUDE.md`.

## Build and run

- Copy the communication library for the OS from `lib/<os>/` to `lib/` (README). On Windows:
  WSL (`lib/win_wsl`), MinGW (`lib/win_mingw`) or Visual Studio (`lib/win_visual_studio`).
- Build: `cmake .` then `make` (or `cmake -G "MSYS Makefiles"` + `cmake --build .` on MinGW).
  `build/` contains a stale macOS configuration; do not reuse it on Windows.
- Run order: start the agent first (UDP server on `127.0.0.1:30000`), then `python simulator.py`
  in `../my_PyDrivingSim`. Ctrl-C stops the agent.
- Logs: `logger.log_var(file, var, value)` + `logger.write_line(file)` write CSV to
  `log_internal/<file>.csv`. MATLAB import/plot: `matlab/importfile_vel_acc_coefs.m`, `matlab/plots.m`.

## Code map

- `starting_point.cc`: main loop at `DT = 0.05` s (matches the simulator agent `dt`).
  Receives `input_data_str`, fills `output_data_str`:
  `RequestedAcc` is the **requested pedal** (not an acceleration). `RequestedSteerWhlAg` [rad, +CCW]
  is read as a **steering-wheel angle** despite the header comment: front-wheel angle = request / tau_D
  (tau_D = 10, `vehicle_double_track_model.py` delta_Dot), while `SteerWhlAg` feedback is the front-wheel
  angle. Send `VEH_STEER_RATIO * delta` (`vehicle_params.h`).
- `lib/interfaces_data_structs.h`: message structs. DO NOT MODIFY (precompiled lib).
- `matlab/student_def_optimal_control.m`: OCP for min ∫j² (project_info slides 23–39),
  exports `s/v/a/j_opt`, `coef_list_fun`, `coeffs_*_opt`, `total_cost`.
- `matlab/student_def_{stop,pass}_primitive{,_j0}.m`: derivations (slides 47–61, 86–92).
- `matlab/student_{stop,pass}_primitive{,_j0}.m`: `%#codegen` primitives, compiled through
  `primitives.prj` into `matlab/codegen/lib/primitives/`.
- `hlp/`: high-level planner. `rrt_star.{h,cc}` RRT* (06 deck slides 19-26, Karaman & Frazzoli 2010),
  fixed 4000-node pool, circle obstacles inflated by the vehicle circle + 0.3 m; `pose_estimator.{h,cc}`
  pose in W0 (= vehicle frame at the first cycle) from static objects as landmarks + `LaneHeading`;
  `hlp.{h,cc}` message (ObjID 3 target, 1/4 obstacles) -> plan, path -> `TrajectoryPointIX/IY` (23 pts).
  Called from `starting_point.cc` (plans once on the first cycle). Offline test: target `rrt_star_test`.
  Pre-run map preview (2026-11-04): scenario `AutomationLevel < 0` = preview request for page -(AutomationLevel+1);
  `starting_point.cc` answers with `hlp_path_page` (whole path, 22-vertex stride), `ManoeuverType = HLP_PREVIEW_ACK + status`,
  and skips control and logs. Protocol in my_PyDrivingSim `scenarios/maps/README.md`.
- `llc/longitudinal.{h,cc}`: action priming/selection (FreeFlow, Stop at the RRT* goal) and LLC
  (trapezoidal integrated jerk with internal a0 + saturation, PI on the pedal with clamp/anti-windup).
- `llc/lateral.{h,cc}`: kinematic Stanley (05 deck p. 43-44 / printed 37-38, Hoffmann et al. 2007),
  chosen by Ale 2026-11-04. P_F = CoM + L_f heading; P_P = projection on the RRT* polyline segments
  (monotonic segment index, `search_ahead` window); signed e > 0 = path on the left; v_f = VLgtFild / cos(SteerWhlAg).
  TODO(dev): `k_e`, `v_min`, Chimera steering limit (fallback 0.5 rad), steering steps at RRT* vertices
  (clothoid smoothing or the k_s term of p. 48). Logs `log_internal/lat.csv`. Offline test: target `lateral_test`
  (kinematic plant of p. 16 + tau_H lag, 4 start offsets, converges).
- `matlab/generate_primitives.m`: runs all `student_def_*` scripts and regenerates `codegen/lib/primitives`
  (replaces the GUI `primitives.prj`, whose paths point to old folders). Run it after any `.m` change.
- Coefficient vector convention (from `coef_list_fun`): 6 entries `[0, c1=v0, c2=a0, c3=j0, c4, c5]`,
  so in C `coefs[3]` is c3 = j(0) (used by mirroring, slide 81, and main logic slide 106).

## What the scenario message gives (and does not)

- Ego: `VLgtFild`, `ALgtFild`, `YawRateFild`, `SteerWhlAg`, `VehicleLen`, `VehicleWidth`,
  `RequestedCruisingSpeed`.
- Road: `LaneWidth`, `LatOffsLineR/L`, `LaneHeading` (= −yaw in the simulator).
  `LaneCrvt` and `AdasisCurvature*` are not filled by the current simulator.
- Objects: up to 20, `ObjX/ObjY` (vehicle frame), `ObjLen/ObjWidth`, `ObjID` (1 = cone, 2 = coin).
- Traffic light: `NrTrfLights`, `TrfLightDist`, state and next three change times (slides 95–97).
- **No absolute pose.** RRT\* and lateral control need a pose estimate: dead-reckon from
  `VLgtFild`/`YawRateFild` (and/or `LaneHeading`, lateral offsets), or plan in the vehicle frame
  each cycle. Chosen 2026-10-28: `hlp/pose_estimator` (static objects as landmarks, same index every
  cycle per the simulator, `LaneHeading` for heading, dead reckoning fallback). NOT IN SLIDES.

## Known issues (audit 2026-10-28; status updated 2026-11-04)

Fixed 2026-11-04 in the sources (needs `matlab/generate_primitives.m` run in MATLAB before building):
1. `student_stop_primitive.m` passed sf = 0 to `coef_list_fun`; now `maxsf` (slide 50 line 11).
2. `student_def_optimal_control.m` coefficient indexing (`c(k+1)` is c_k, no extra `c(6)*t^5`); Inspect section.
3. `coeffs_a_opt`/`coeffs_j_opt` missing from the C library: now codegen entry points.
4./5. LLC rewritten in `llc/longitudinal.cc`: persistent `e_int`, trapezoidal integrated jerk with internal a0
   and saturation (slide 68) instead of `aEval(DT, m*)`.
6. `student_pass_primitive` call: separate m1/m2 buffers, FreeFlow(v0, a0, lookahead, vr) (slide 104).
   `.m` now returns `zeros(1,6)` rows -> fixed-size `double[6]` in C (was variable-size `_data/_size`).
7. Logs read `m_star`, not an uninitialised array.
8. `student_def_pass_primitive.m` picks the positive root T* = sqrt(15 sf/(-a0)) for both T* and v*.
12. `<chrono>` included; unused variables removed.
New, fixed: FreeFlow (vmin == vmax) returned no primitive whenever a0 < 0 (strict `vfmin < vfmax`);
free-flow branch had T1/T2 swapped; `student_pass_primitive_j0.m` function name did not match the file,
no guard on the sqrt argument or T > 0.

Still open:
9. [RISK] Slide 59 pseudo-code v* = 1/4(...) is a typo; derivation gives slide 58's 1/8(2√15√(−a0 sf) − 7v0)
   (checked by hand: dvf/dT = 0 -> T*², substitute). `min_vel.m` is right.
10. [RISK] Free-flow branch `Tmin == Tmax == 0` is a local reading of slide 104 (intersection with [0,0] empty).
11. [RISK] Off-road action selection is FreeFlow + Stop at the RRT* goal only; no PassJ0/StopJ0/mirroring use,
    no curvature speed limit. Parameters in `lgt_default_params` are TODO(dev).
13. `.m` sources still uncommitted; `primitives.prj` paths are stale (use `generate_primitives.m`).
14. [RISK] `final_opt_time_pass_j0` root order differs from slide 90 (code tries the larger-denominator root first).
