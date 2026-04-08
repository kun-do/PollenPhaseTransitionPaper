cmake_minimum_required(VERSION 3.16)

project(PollenPhaseTransitionPaper LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(POLLEN_ENABLE_SOLVER_SMOKE_TESTS "Enable long-running solver smoke tests." OFF)

find_package(Eigen3 QUIET CONFIG)

if(NOT TARGET Eigen3::Eigen)
    set(POLLEN_EIGEN3_INCLUDE_DIR "" CACHE PATH "Path to the Eigen headers (directory containing Eigen/Dense).")
    find_path(
        POLLEN_EIGEN3_INCLUDE_DIR
        NAMES Eigen/Dense
        HINTS ${POLLEN_EIGEN3_INCLUDE_DIR}
        PATH_SUFFIXES eigen3
    )

    if(NOT POLLEN_EIGEN3_INCLUDE_DIR)
        message(FATAL_ERROR "Eigen3 headers were not found. Install Eigen3 or set POLLEN_EIGEN3_INCLUDE_DIR.")
    endif()

    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(
        Eigen3::Eigen
        PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${POLLEN_EIGEN3_INCLUDE_DIR}"
    )
endif()

set(POLLEN_COMMON_DIR "${CMAKE_CURRENT_SOURCE_DIR}/PhaseDiagramCalculations/Common")
set(POLLEN_COMMON_SOURCES "${POLLEN_COMMON_DIR}/wigner_compat.cpp")

function(add_pollen_solver target source_dir)
    add_executable(
        ${target}
        "${source_dir}/main.cpp"
        ${POLLEN_COMMON_SOURCES}
    )

    target_include_directories(
        ${target}
        PRIVATE
        "${source_dir}"
        "${POLLEN_COMMON_DIR}"
    )

    target_link_libraries(${target} PRIVATE Eigen3::Eigen)
endfunction()

add_pollen_solver(
    gradient_descent_1l
    "${CMAKE_CURRENT_SOURCE_DIR}/PhaseDiagramCalculations/GradientDescent_1l"
)
add_pollen_solver(
    gradient_descent_2l
    "${CMAKE_CURRENT_SOURCE_DIR}/PhaseDiagramCalculations/GradientDescent_2l"
)
add_pollen_solver(
    simulated_annealing_1l
    "${CMAKE_CURRENT_SOURCE_DIR}/PhaseDiagramCalculations/SimulatedAnnealing_1l"
)
add_pollen_solver(
    simulated_annealing_2l
    "${CMAKE_CURRENT_SOURCE_DIR}/PhaseDiagramCalculations/SimulatedAnnealing_2l"
)

include(CTest)

if(BUILD_TESTING)
    add_executable(
        wigner_compat_test
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/wigner_compat_test.cpp"
        ${POLLEN_COMMON_SOURCES}
    )

    target_include_directories(
        wigner_compat_test
        PRIVATE
        "${POLLEN_COMMON_DIR}"
    )

    add_test(NAME wigner_compat_test COMMAND wigner_compat_test)

    if(POLLEN_ENABLE_SOLVER_SMOKE_TESTS)
        foreach(target IN ITEMS gradient_descent_1l gradient_descent_2l simulated_annealing_1l simulated_annealing_2l)
            set(test_directory "${CMAKE_CURRENT_BINARY_DIR}/smoke/${target}")
            file(MAKE_DIRECTORY "${test_directory}")
            add_test(NAME ${target}_smoke COMMAND ${target})
            set_tests_properties(
                ${target}_smoke
                PROPERTIES
                LABELS "solver;smoke"
                TIMEOUT 300
                WORKING_DIRECTORY "${test_directory}"
            )
        endforeach()
    endif()
endif()
