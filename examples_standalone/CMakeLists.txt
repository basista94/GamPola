cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_C_COMPILER "gcc")

project(Gamapola)
set(CMAKE_BUILD_TYPE "Debug")

ADD_DEFINITIONS( "-Wl,--no-as-needed" )
ADD_DEFINITIONS( "`root-config --cflags --libs`" )
ADD_DEFINITIONS( "-O2 -W -Wall -Wextra -Wno-long-long -Wno-deprecated" )
ADD_DEFINITIONS( "-fno-common" )
ADD_DEFINITIONS( "-fopenmp" )
ADD_DEFINITIONS( "-std=c++11" )

include_directories(${ROOT_INCLUDE_DIRS})

include_directories(${PROJECT_SOURCE_DIR}/../include)
include_directories(${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR})

link_directories($ENV{ROOTSYS}/lib)
link_directories(${PROJECT_SOURCE_DIR}/../build)


add_executable(generator generator.cxx)
target_link_libraries(generator -lGamPola -lginac -lCore)

add_executable(calculateNormInts calculateNormInts.cxx)
target_link_libraries(calculateNormInts -lGamPola -lginac -lCore)

add_executable(fitter fitter.cxx)
target_link_libraries(fitter -lGamPola -lginac -lCore)

add_executable(plotter plotter.cxx)
target_link_libraries(plotter -lGamPola -lginac -lCore -lRIO -lHist)

add_executable(bf2 bf2.cxx)
target_link_libraries(bf2 -lGamPola -lginac -lCore -lRIO -lHist)
