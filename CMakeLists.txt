cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_C_COMPILER "gcc")

project(Gamapola)
set(CMAKE_BUILD_TYPE "Debug")

ADD_DEFINITIONS( "-Wl,--no-as-needed" )
ADD_DEFINITIONS( "`root-config --libs --cflags`" )
ADD_DEFINITIONS( "-O2 -W -Wall -Wextra -Wno-long-long -Wno-deprecated" )
ADD_DEFINITIONS( "-fno-common" )
ADD_DEFINITIONS( "-fopenmp" )
ADD_DEFINITIONS( "-std=c++11" )

#find_package(root)
FIND_PACKAGE( Boost 1.54 COMPONENTS program_options REQUIRED )
include_directories(${ROOT_INCLUDE_DIRS})

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR})
include_directories(${PROJECT_SOURCE_DIR}/include/GamPolaParser/)
include_directories(${PROJECT_SOURCE_DIR}/include/GamPolaPlotter/)
include_directories(${PROJECT_SOURCE_DIR}/include/EvtGenModels/)
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

link_directories($ENV{ROOTSYS}/lib)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cxx)
file(GLOB sourcesEvtGenModels ${PROJECT_SOURCE_DIR}/src/EvtGenModels/*.cpp)

file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB headersParser ${PROJECT_SOURCE_DIR}/include/GamPolaParser/*.h)
file(GLOB headersPlotter ${PROJECT_SOURCE_DIR}/include/GamPolaPlotter/*.h)
file(GLOB headersEvtGenModels ${PROJECT_SOURCE_DIR}/include/EvtGenModels/*.hh)

add_library(GamPola SHARED ${sources} ${sourcesEvtGenModels} ${headers} ${headersParser} ${headersPlotter} ${headersEvtGenModels})
target_link_libraries(GamPola  -lm -lcln -lginac -lstdc++ -fopenmp -pthread -I/usr/local/include/root -L/usr/lib/x86_64-linux-gnu/ -L/usr/aarch64-linux-gnu/lib -L/usr/local/lib/root -L/usr/local/lib/ -lcln -lginac -lTree -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lMinuit2 -lMinuit -lMathMore -lThread -lm -rdynamic -L/usr/local/lib -lboost_filesystem -lboost_system -ldl )

add_library(GamPolaStatic STATIC ${sources} ${sourcesEvtGenModels} ${headers} ${headersParser} ${headersPlotter} ${headersEvtGenModels})
target_link_libraries(GamPolaStatic  -lm -lcln -lginac -lstdc++ -fopenmp -pthread -I/usr/local/include/root -L/usr/lib/x86_64-linux-gnu/ -L/usr/aarch64-linux-gnu/lib -L/usr/local/lib/root -L/usr/local/lib/ -lcln -lginac -lTree -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lMinuit2 -lMinuit -lMathMore -lThread -lm -rdynamic -L/usr/local/lib -lboost_filesystem -lboost_system -ldl )

