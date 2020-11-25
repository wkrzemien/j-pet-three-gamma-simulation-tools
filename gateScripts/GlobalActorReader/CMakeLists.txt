cmake_minimum_required(VERSION 2.6)

set(projectName GAR)
set(projectBinary ${projectName}.x)

project(${projectName} CXX) # using only C++

add_definitions(-std=c++11 -Wall) # using C++11 standard

find_package(ROOT REQUIRED COMPONENTS  
  Hist 
  HistPainter 
  Physics 
  RIO 
  Thread 
  Gpad 
  Tree 
  Minuit 
  Html 
  Core 
  RooFit 
  ROOTDataFrame
  MultiProc
  ROOTVecOps
  RooFitCore 
  TreePlayer 
  MathMore 
  ) 
include_directories(${ROOT_INCLUDE_DIRS}) 
add_definitions(${ROOT_DEFINITIONS}) 
link_directories(${ROOT_LIBRARY_DIRS})

include(${ROOT_USE_FILE})

################################################################################
# Sources and linking 
################################################################################
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)  #to change the LD_LIBRARY_PATH and accept pcm

file(GLOB HEADERS *.hh)
file(GLOB SOURCES *.cc)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
ROOT_GENERATE_DICTIONARY(EventDict Event.h LINKDEF EventLinkDef.h)
#---Create a shared library with geneated dictionary
add_library(Event SHARED Event.cxx EventDict.cxx ${ROOT_LIBRARIES})
target_link_libraries(Event ${ROOT_LIBRARIES})


add_executable(GAR ${SOURCES} Event.cxx)
target_link_libraries(GAR ${ROOT_LIBRARIES} Event)

add_library(readTree SHARED readTree.cpp ${ROOT_LIBRARIES})
add_executable(start.exe readTree.cpp Event.cxx)
target_link_libraries(start.exe PRIVATE Event ${ROOT_LIBRARIES} readTree) 

