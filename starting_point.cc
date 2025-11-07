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
    static double req_acc = 4.0;   // Placeholder for requested acceleration for testing
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
            double sf = in->TrfLightDist;

            double acc = in->ALgtFild;  // Current vehicle acceleration
            double vel = in->VLgtFild;  // Current vehicle velocity
            
            // student_pass_primitive();
            
            
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