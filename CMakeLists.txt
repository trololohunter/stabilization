cmake_minimum_required(VERSION 3.8)
project(untitled2)
link_libraries(m)
set(CMAKE_C_STANDARD 99)

set(SOURCE_FILES main.c gas_two.c case.c gas_two.h residuals.h residuals.c functions.h functions.c case.h
        laspack/itersolv.h laspack/itersolv.c
        laspack/eigenval.h laspack/eigenval.c
        laspack/operats.h laspack/operats.c
        laspack/mlsolv.h laspack/mlsolv.c
        laspack/errhandl.h laspack/errhandl.c
        laspack/factor.h laspack/factor.c
        laspack/matrix.h laspack/matrix.c
        laspack/vector.h laspack/vector.c
        laspack/qmatrix.h laspack/qmatrix.c
        laspack/lastypes.h
        laspack/precond.h laspack/precond.c
        laspack/rtc.h laspack/rtc.c
        laspack/elcmp.h gnuploting.h gnuploting.c)
add_executable(untitled2 ${SOURCE_FILES})