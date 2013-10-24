cmake_minimum_required (VERSION 2.6)
project (sqct)
add_definitions("-std=c++0x")

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} )
set( CMAKE_CXX_FLAGS_RELEASE "-O3")
#set( CMAKE_EXE_LINKER_FLAGS "-static-libstdc++ -static-libgcc -static")
#set(Boost_USE_STATIC_LIBS ON)

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
 hprhelpers.cpp
 resring.cpp
 rint.cpp
 matrix2x2.cpp
 vector2.cpp
 output.cpp
 gatelibrary.cpp
 symbolic_angle.cpp
 
 es/numbersgen.cpp
 es/optsequencegenerator.cpp
 es/seqlookupcliff.cpp
 es/exactdecomposer.cpp
 
 appr/zrot_cache.cpp
 appr/toptzrot2.cpp
 appr/normsolver.cpp
 appr/approxlist.cpp

 test.cpp
 factorzs2.cpp
 tcount.cpp
 solvenormequation.cpp

 mainA.cpp
)

set(BOOST_LIBS program_options timer chrono system)
find_package(Boost COMPONENTS ${BOOST_LIBS} REQUIRED)

target_link_libraries(sqctSK ${Boost_LIBRARIES})
target_link_libraries(sqctSK gomp pthread mpfr gmpxx gmp rt)

target_link_libraries(sqct ${Boost_LIBRARIES})
target_link_libraries(sqct gomp pthread pari mpfr gmpxx gmp rt )
