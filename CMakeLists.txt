cmake_minimum_required(VERSION 3.5)
project(basic_agent_ST)

set(CMAKE_CXX_STANDARD 11)
# If you are on mac x86
#set(CMAKE_OSX_ARCHITECTURES x86_64)

include(lib/basic.cmake)

include_directories(
        matlab/codegen/lib/primitives
        matlab
        lib
        log
        hlp
)

link_directories(
        lib
)

# ------------------------------- basic agent ------------------------------------
add_executable(basic_agent
        log/screen_print.c
        log/logvars.cc
        starting_point.cc
        matlab/codegen/lib/primitives/primitives.c
        matlab/codegen/lib/primitives/rt_nonfinite.c
        hlp/rrt_star.cc
        hlp/pose_estimator.cc
        hlp/hlp.cc
    )


target_link_libraries(basic_agent agent_lib ${lib_ws32})

target_compile_options(basic_agent
        PUBLIC -DNOT_INLINED -DSERVER
        )

set_target_properties(basic_agent PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${AGENT_FOLDER})
# ------------------------------------------------------------------------------------

# ------------------------- RRT* offline test (no simulator) -------------------------
add_executable(rrt_star_test
        hlp/rrt_star_test.cc
        hlp/rrt_star.cc
    )
set_target_properties(rrt_star_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${AGENT_FOLDER})
# ------------------------------------------------------------------------------------
