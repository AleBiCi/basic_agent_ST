//
// Created by Gastone Pietro Rosati Papini on 10/08/22.
//

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <chrono>

extern "C" {
#include "screen_print_c.h"
}
#include "screen_print.h"
#include "server_lib.h"
#include "logvars.h"

// --- MATLAB PRIMITIVES INCLUDE ---
#include "primitives.h"
// --- MATLAB PRIMITIVES INCLUDE ---

// --- HIGH-LEVEL PLANNER (RRT*) ---
#include "hlp.h"

// --- LONGITUDINAL CONTROL (primitives + LLC) ---
#include "longitudinal.h"
// --- LATERAL CONTROL (kinematic Stanley) ---
#include "lateral.h"
#include "vehicle_params.h"

#define DEFAULT_SERVER_IP    "127.0.0.1"
#define SERVER_PORT               30000  // Server port
#define DT 0.05

// Handler for CTRL-C
#include <signal.h>
static uint32_t server_run = 1;
void intHandler(int signal) {
    server_run = 0;
}

int main(int argc, const char * argv[]) {
    logger.enable(true);

    // Messages variables
    scenario_msg_t scenario_msg;
    manoeuvre_msg_t manoeuvre_msg;
    size_t scenario_msg_size = sizeof(scenario_msg.data_buffer);
    size_t manoeuvre_msg_size = sizeof(manoeuvre_msg.data_buffer);
    uint32_t message_id = 0;

#if not defined( _MSC_VER ) and not defined( _WIN32 )
    // More portable way of supporting signals on UNIX
    struct sigaction act;
    act.sa_handler = intHandler;
    sigaction(SIGINT, &act, NULL);
#else
    signal(SIGINT, intHandler);
#endif

    server_agent_init(DEFAULT_SERVER_IP, SERVER_PORT);

    // Start server of the Agent
    printLine();
    printTable("Waiting for scenario message...", 0);
    printLine();

    while (server_run == 1) {

        // Clean the buffer
        memset(scenario_msg.data_buffer, '\0', scenario_msg_size);

        // Receive scenario message from the environment
        if (server_receive_from_client(&server_run, &message_id, &scenario_msg.data_struct) == 0) {
            // Init time
            static auto start = std::chrono::system_clock::now();
            auto time = std::chrono::system_clock::now()-start;
            double num_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(time).count()/1000.0;
            printLogTitle(message_id, "received message");

            // Data struct
            input_data_str *in = &scenario_msg.data_struct;
            manoeuvre_msg.data_struct.CycleNumber = in->CycleNumber;
            manoeuvre_msg.data_struct.Status = in->Status;

            // Example of using log
            logger.log_var("Example", "cycle", in->CycleNumber);
            logger.log_var("Example", "vel", in->VLgtFild);
            
            
        // ADD AGENT CODE HERE
            // --- HIGH-LEVEL PLANNER: pose in W0, RRT* to the target (ObjID 3), path -> trajectory points ---
            static PoseEstimator pose;
            static RrtStar rrt;             // ~170 kB node pool: static, not on the stack
            static RrtPath path;
            static int hlp_status = HLP_NOT_RUN;
            static int path_idx = 0;
            if (!pose.initialised()) pose.init(in); else pose.update(in, DT);
            if (hlp_status == HLP_NOT_RUN) {
                // Plan once from the first pose. TODO(dev): replanning policy (target moved, path blocked).
                hlp_status = hlp_plan(in, pose, hlp_default_params(in), rrt, &path);
                printLogVar(message_id, "RRT* status", hlp_status);
                printLogVar(message_id, "RRT* length", path.length);
            }
            // --- PRE-RUN MAP PREVIEW (my_PyDrivingSim MapScenario, vehicle frozen until Start) ---
            //     AutomationLevel = -(page + 1) asks for one page of the whole path (hlp_path_page);
            //     no control, no logs. ManoeuverType = HLP_PREVIEW_ACK + HLP status.
            if (in->AutomationLevel < 0) {
                manoeuvre_msg.data_struct.NTrajectoryPoints = 0;
                if (hlp_status == RRTS_OK)
                    hlp_path_page(&path, pose, -in->AutomationLevel - 1, &manoeuvre_msg.data_struct);
                manoeuvre_msg.data_struct.ManoeuverType = HLP_PREVIEW_ACK + hlp_status;
                manoeuvre_msg.data_struct.RequestedAcc = 0.;
                manoeuvre_msg.data_struct.RequestedSteerWhlAg = 0.;
                if (server_send_to_client(server_run, message_id, &manoeuvre_msg.data_struct) == -1) {
                    perror("error send_message()");
                    exit(EXIT_FAILURE);
                }
                continue;
            }
            if (hlp_status == RRTS_OK) {
                hlp_path_to_trajectory(&path, pose, &path_idx, &manoeuvre_msg.data_struct);
            } else {
                manoeuvre_msg.data_struct.NTrajectoryPoints = 0;
            }
            logger.log_var("hlp", "time", in->ECUupTime);
            logger.log_var("hlp", "x", pose.x);
            logger.log_var("hlp", "y", pose.y);
            logger.log_var("hlp", "psi", pose.psi);
            logger.log_var("hlp", "landmarks", pose.landmarks_used);
            logger.log_var("hlp", "gps", pose.gps_used ? 1 : 0);
            logger.log_var("hlp", "path_idx", path_idx);
            logger.write_line("hlp");

            // --- LONGITUDINAL: action priming + selection (project_info slides 47-60, 104-106),
            //     LLC: trapezoidal integrated jerk with internal a0 + PI on the pedal (slides 68-69) ---
            const double v0 = in->VLgtFild;     // slide 104 line 1
            const double a0 = in->ALgtFild;     // slide 104 line 2
            static LgtLlc llc = {false, 0., 0., 0.};
            const LgtParams lp = lgt_default_params(in);

            // Stop point = RRT* goal: remaining path length, else straight line to ObjID 3 (-1 = none)
            double s_target = (hlp_status == RRTS_OK) ? lgt_path_remaining(&path, pose, path_idx)
                                                      : lgt_target_distance(in);
            double m_star[6];
            double T_star;
            int behaviour = lgt_select(v0, a0, s_target, lp, m_star, &T_star);
            double req_ped = lgt_llc_step(&llc, behaviour != LGT_NONE ? m_star : nullptr,
                                          v0, a0, DT, lp);

            manoeuvre_msg.data_struct.RequestedAcc = req_ped;   // NOTE: RequestedAcc === REQUESTED PEDAL
            // --- LATERAL: kinematic Stanley on the RRT* path (05_Low_Level_MotionPlanning p. 43-44) ---
            static LatState lat = {0, 0., 0., 0.};
            const LatParams latp = lat_default_params();
            double delta_req = 0.;      // [rad] front-wheel angle
            int lat_status = LAT_NO_PATH;
            if (hlp_status == RRTS_OK)
                lat_status = lat_stanley_step(&lat, &path, pose,
                                              lat_front_speed(in->VLgtFild, in->SteerWhlAg), latp, &delta_req);
            // The simulator reads RequestedSteerWhlAg as a steering-wheel angle (front wheel = request / tau_D)
            manoeuvre_msg.data_struct.RequestedSteerWhlAg = VEH_STEER_RATIO * delta_req;

            logger.log_var("lat", "time", in->ECUupTime);
            logger.log_var("lat", "status", lat_status);
            logger.log_var("lat", "seg", lat.seg);
            logger.log_var("lat", "e", lat.e);
            logger.log_var("lat", "dpsi", lat.dpsi);
            logger.log_var("lat", "delta_req", delta_req);
            logger.log_var("lat", "delta_meas", in->SteerWhlAg);
            logger.log_var("lat", "yaw_rate", in->YawRateFild);
            logger.write_line("lat");

            // LOG FILE FOR LLC TESTING
            logger.log_var("acc_vel_test", "time", in->ECUupTime);
            logger.log_var("acc_vel_test", "acc", a0);
            logger.log_var("acc_vel_test", "req_acc", llc.a_req);
            logger.log_var("acc_vel_test", "vel", v0);
            logger.log_var("acc_vel_test", "req_ped", req_ped);
            logger.log_var("acc_vel_test", "behaviour", behaviour);
            logger.log_var("acc_vel_test", "s_target", s_target);
            logger.log_var("acc_vel_test", "T", T_star);
            logger.log_var("acc_vel_test", "c1", m_star[1]);
            logger.log_var("acc_vel_test", "c2", m_star[2]);
            logger.log_var("acc_vel_test", "c3", m_star[3]);
            logger.log_var("acc_vel_test", "c4", m_star[4]);
            logger.log_var("acc_vel_test", "c5", m_star[5]);
            logger.write_line("acc_vel_test");  // Writes the actual log file

            // Write log
            // logger.write_line("Example");

            // Screen print
            printLogVar(message_id, "Time", num_seconds);
            printLogVar(message_id, "Status", in->Status);
            printLogVar(message_id, "acc", a0);  // Plot acceleration on console log
            printLogVar(message_id, "vel", v0);  // Plot velocity on console log
            printLogVar(message_id, "behaviour", behaviour);
            printLogVar(message_id, "CycleNumber", in->CycleNumber);

            // Send manoeuvre message to the environment
            if (server_send_to_client(server_run, message_id, &manoeuvre_msg.data_struct) == -1) {
                perror("error send_message()");
                exit(EXIT_FAILURE);
            } else {
                printLogTitle(message_id, "sent message");
            }
        }
    }

    // Close the server of the agent
    server_agent_close();
    return 0;
}