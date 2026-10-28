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
  `RequestedAcc` is the **requested pedal** (not an acceleration), `RequestedSteerWhlAg` the
  **front-wheel angle** [rad, +CCW]. The simulator applies both directly
  (`pydrivingsim/agent.py`: `self.action = (m.RequestedAcc, m.RequestedSteerWhlAg)`).
- `lib/interfaces_data_structs.h`: message structs. DO NOT MODIFY (precompiled lib).
- `matlab/student_def_optimal_control.m`: OCP for min ∫j² (project_info slides 23–39),
  exports `s/v/a/j_opt`, `coef_list_fun`, `coeffs_*_opt`, `total_cost`.
- `matlab/student_def_{stop,pass}_primitive{,_j0}.m`: derivations (slides 47–61, 86–92).
- `matlab/student_{stop,pass}_primitive{,_j0}.m`: `%#codegen` primitives, compiled through
  `primitives.prj` into `matlab/codegen/lib/primitives/`.
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
  each cycle. This is a design decision for the developer; ask before choosing.

## Known issues (audit 2026-10-28, not yet fixed; verify before relying on them)

Longitudinal control is the "to correct" part. Findings against the slides:

1. [BUG] `matlab/student_stop_primitive.m:27`: `coef_list_fun(v0, a0, 0., 0., 0., tf)` passes
   `sf = 0`; slide 50 line 11 uses `evalPrimitiveCoeffs(v0, a0, smax, 0, 0, tf)`. Same in
   generated `primitives.c`.
2. [BUG] `student_def_optimal_control.m` builds `coeffs_s_opt` with an extra `c(6)*t^5` term and
   indexes `c(1..5)` while `coef_list_fun` returns `[0, v0, a0, c3, c4, c5]`. Result:
   `coeffs_a_opt(t, coefs)` returns ≈ v0 + a0·t, not a(t). This is likely why the logged
   accelerations needed capping in `plots.m`. Also the "Inspect" section uses
   `coeffs_v_opt_fun` for acceleration and `coeffs_a_opt_fun` for jerk.
3. [BUG] `starting_point.cc` calls `coeffs_a_opt`, which is not in the generated `primitives.h`
   (not part of the Coder project), so the current file should not link. Check `primitives.prj`.
4. [BUG] `starting_point.cc`: `double eint = 0.0;` is inside the loop, so the PI integral resets
   every cycle (slide 69 needs persistent `e_int`).
5. [BUG] LLC uses `a_req = aEval(DT, m*)`, the method slide 65 calls "trivial but wrong / not
   robust". Slides 66–68: integrated jerk, trapezoidal, with internal a0 and saturation.
6. [BUG] `student_pass_primitive(...)` in `starting_point.cc` passes the same `coefs, &fvel, &fdist`
   buffers for both outputs, so m1 overwrites m2; the third output is a time `T`, not a distance.
   It is also called with `vmin = vmax = 15` instead of `RequestedCruisingSpeed` and with the
   traffic-light distance instead of the lookahead (slide 104: `FreeFlow(v0,a0,lookahead,vr)`).
7. [BUG] Logged `c1..c5` read the uninitialised array `coef`, not `coefs`.
8. [BUG] `time_min_vel.m` returns `sqrt(15)·sqrt(−a0·sf)/a0`, which is negative for a0 < 0;
   slide 58: `T* = sqrt(15 sf / (−a0))`. The derivation picks root (1) for T* but root (2) for v*.
9. [RISK] Slide 59 pseudo-code shows `v* = 1/4(√15·√(−a0 sf) − 7v0)` while slide 58 gives
   `v* = 1/8(2√15·√(−a0 sf) − 7v0)`; `min_vel.m` matches slide 58. Confirm with the professor.
10. [RISK] Free-flow branch (`Tmin == Tmax == 0`) in `student_pass_primitive.m` is a local
    interpretation of slide 104; the slide intersection `[0,0] ∩ [Tvmax,Tvmin]` would be empty.
11. [RISK] No action selection yet (stop if `dist < 50`, else pass), no traffic-light state use
    (slides 99–106), no j0=0 primitives in the loop.
12. [STYLE] `<chrono>` not included; `s`, `v0`, `a0`, `req_vel` partly unused.
13. 21 MATLAB files have uncommitted changes and 3 are untracked; `primitives.c` dates from
    2025-12-15. Regenerate codegen after fixing the `.m` sources.
