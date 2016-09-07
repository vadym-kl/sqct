cmake_minimum_required (VERSION 2.6)
project (sqct)
add_definitions("-std=c++11 -fopenmp")

#INCLUDE_DIRECTORIES(/usr/local/include)
#LINK_DIRECTORIES(/opt/local/lib)
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} )


set( CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG ")
set( CMAKE_EXE_LINKER_FLAGS "-static-libstdc++ -static-libgcc -static")
set(Boost_USE_STATIC_LIBS ON)

add_executable(sqctSK
 hprhelpers.cpp
 resring.cpp
 rint.cpp
 matrix2x2.cpp
 vector2.cpp
 output.cpp
 gatelibrary.cpp
 
 es/numbersgen.cpp
 es/optsequencegenerator.cpp
 es/seqlookupcliff.cpp
 es/exactdecomposer.cpp

 sk/vector3hpr.cpp
 sk/epsilonnet.cpp
 sk/netgenerator.cpp
 sk/unitaryapproximator.cpp
 sk/gcommdecomposer.cpp
 sk/sk.cpp

 theory/numbers-stat.cpp
 theory/theoremverification.cpp
 theory/toptimalitytest.cpp
 theory/hoptimalitytest.cpp

 main.cpp
)

add_executable(sqct
 requestprocessor.cpp
 hprhelpers.cpp
 resring.cpp
 rint.cpp
 matrix2x2.cpp
 vector2.cpp
 output.cpp
 gatelibrary.cpp
 symbolic_angle.cpp
 serializers.h
 timemeasurement.h

 
 es/numbersgen.cpp
 es/optsequencegenerator.cpp
 es/seqlookupcliff.cpp
 es/exactdecomposer.cpp
 
 appr/zrot_cache.cpp
 appr/toptzrot2.cpp
 appr/normsolver.cpp
 appr/approxlist.cpp
 appr/findhalves.cpp
 appr/rcup.cpp
 appr/topt-bfs.cpp
 appr/cup.cpp
 appr/observations.cpp

 test.cpp
 factorzs2.cpp
 tcount.cpp
 solvenormequation.cpp
 fixedpoint.cpp

 mainA.cpp
)



find_package(Boost 1.54.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(BOOST_LIBS program_options timer chrono system )
find_package(Boost COMPONENTS ${BOOST_LIBS} REQUIRED)

MESSAGE("Used boost libraries:")
MESSAGE("${Boost_LIBRARIES}")

target_link_libraries(sqct ${Boost_LIBRARIES})
target_link_libraries(sqct gomp pthread pari mpfr gmpxx gmp rt )

target_link_libraries(sqctSK ${Boost_LIBRARIES})
target_link_libraries(sqctSK gomp pthread mpfr gmpxx gmp rt)

