cmake_minimum_required (VERSION 2.6)
project (sqct)
add_definitions("-std=c++0x")
set( CMAKE_CXX_FLAGS_RELEASE "-O3")
set( CMAKE_EXE_LINKER_FLAGS "-static-libstdc++ -static-libgcc -static")
set(Boost_USE_STATIC_LIBS ON)

add_executable(sqct
 hprhelpers.cpp
 resring.cpp
 rint.cpp
 matrix2x2.cpp
 vector2.cpp
 vector3hpr.cpp
 output.cpp

 numbersgen.cpp
 gatelibrary.cpp
 optsequencegenerator.cpp
 seqlookupcliff.cpp
 exactdecomposer.cpp

 epsilonnet.cpp
 netgenerator.cpp
 unitaryapproximator.cpp
 gcommdecomposer.cpp
 sk.cpp

 numbers-stat.cpp
 theoremverification.cpp
 toptimalitytest.cpp
 hoptimalitytest.cpp

 main.cpp
)

set(BOOST_LIBS program_options timer chrono system)
find_package(Boost COMPONENTS ${BOOST_LIBS} REQUIRED)

message("XXX ${Boost_LIBRARIES}")

target_link_libraries(sqct ${Boost_LIBRARIES})
target_link_libraries(sqct gomp pthread mpfr gmpxx gmp rt)
