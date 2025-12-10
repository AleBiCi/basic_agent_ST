//
// Created by Gastone Pietro Rosati Papini on 10/08/22.
//

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <algorithm>

extern "C" {
#include "screen_print_c.h"
}
#include "screen_print.h"
#include "server_lib.h"
#include "logvars.h"

// --- MATLAB PRIMITIVES INCLUDE ---
#include "primitives.h"
// --- MATLAB PRIMITIVES INCLUDE ---

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

    /* VALUES FOR PID LLC */
    // static double req_acc = 4.0;   // Placeholder for requested acceleration for testing
    double req_vel = 0.0; // Placeholder for requested velocity

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
            double v0 = 0.;
            double a0 = 0.;
            double dist = in->TrfLightDist;
            static double init_dist = dist;
            double s = init_dist - dist;

            double acc = in->ALgtFild;  // Current vehicle acceleration
            double vel = in->VLgtFild;  // Current vehicle velocity
            double t_curr = in->ECUupTime;  // Current simulation time
            double coef[6]; // Optimal control coefficients
                     
        /*  PRIMITIVES TEST   */
            // double final_time = 20.;
            // double req_acc;
            // // OPTIMAL CONTROL OF req_vel and req_acc
            // req_vel = v_opt(DT, vel, acc, dist, 20, 0, final_time - t_curr);    // COMPUTE for EACH TIMESTEP DT; tf HAS TO UPDATE @ EACH TIME STEP
            // req_acc = a_opt(DT, vel, acc, dist, 20, 0, final_time - t_curr);

            // coef_list_fun(v0, a0, dist, req_vel, req_acc, final_time, coef); // SAVE THESE LATER IN THE LOG

            // ACTUALLY COMPUTING req_acc FROM PRIMITIVES
            double coefs[6];
            // Final (output) values
            double ftime;
            double fdist;
            double fvel;

            if (dist < 50) {
                student_stop_primitive(vel, acc, dist, coefs, &fdist, &ftime);
            } else {
                student_pass_primitive(vel, acc, dist, 15.0, 15.0, 0.0, 0.0, coefs, &fvel, &fdist, coefs, &fvel, &fdist);  
                // w/ init and final time = 0. --> free-running primitive
                // If we have 2 different final velocities (fvmin != fvmax) then we pass different values for &fvel, &fdist and coefs
            }

            double req_acc = coeffs_a_opt(DT, coefs);
            // we compute req_vel based on this req_acc
        /*  PRIMITIVES TEST OVER    */

            
        // ADD LOW LEVEL CONTROL CODE HERE
            // manoeuvre_msg.data_struct.RequestedAcc = -0.3;
            manoeuvre_msg.data_struct.RequestedSteerWhlAg = 0.0;
            
            /*  Testing PID with reference ACCELERATION PROFILE */
            double P_gain = 0.2; // TODO: test PI values
            double I_gain = 0.1;
            
            double error = req_acc - acc;   // Difference between requested acceleration and actual agent accel
            double eint = 0.0;
            eint += error*DT;   // Integral of the error

            /* Requested pedal (Proportional + Integrative controller) */
            double req_ped = P_gain * error + I_gain * eint;

            req_vel += req_acc * DT;    // Recompute requested velocity at every timestep
            
            manoeuvre_msg.data_struct.RequestedAcc = req_ped;   // NOTE: RequestedAcc === REQUESTED PEDAL

        // LOG FILE FOR LLC TESTING (filename = "acc_test")
            logger.log_var("acc_vel_test", "time", in->ECUupTime);  // logs current uptime
            logger.log_var("acc_vel_test", "acc", acc); // current acc
            logger.log_var("acc_vel_test", "req_acc", req_acc); // requested acc
            logger.log_var("acc_vel_test", "vel", vel); // requested vel
            logger.log_var("acc_vel_test", "req_vel", req_vel); // requested vel
            logger.log_var("acc_vel_test", "c1", coef[1]); // resulting coefficients
            logger.log_var("acc_vel_test", "c2", coef[2]); // resulting coefficients
            logger.log_var("acc_vel_test", "c3", coef[3]); // resulting coefficients
            logger.log_var("acc_vel_test", "c4", coef[4]); // resulting coefficients
            logger.log_var("acc_vel_test", "c5", coef[5]); // resulting coefficients

            logger.write_line("acc_vel_test");  // Writes the actual log file

            // Write log
            // logger.write_line("Example");

            // Screen print
            printLogVar(message_id, "Time", num_seconds);
            printLogVar(message_id, "Status", in->Status);
            printLogVar(message_id, "acc", acc);  // Plot acceleration on console log
            printLogVar(message_id, "vel", req_vel);  // Plot velocity on console log
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