cmake_minimum_required(VERSION 3.0)
project(Zgamma)

if(Boost_FOUND)
include_directories(${Boost_INCLUDE_DIRS})
endif()

find_package(ROOT REQUIRED COMPONENTS RIO Net)
include(${ROOT_USE_FILE})
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include"  "Zgamma" "$ENV{HOME}/atlasstyle-00-04-02/")
set(SOURCE_EXE main.cpp)

#file (GLOB_RECURSE CPP_FILES Root/*.cpp $ENV{HOME}/atlasstyle-00-04-02/*.C)
file (GLOB_RECURSE CPP_FILES Root/*.cpp Zgamma/*.h $ENV{HOME}/atlasstyle-00-04-02/*.C)

set(SOURCE_LIB ${CPP_FILES})
add_library(Zgam STATIC ${SOURCE_LIB})
add_executable(main ${SOURCE_EXE})

if(Boost_FOUND)
target_link_libraries(main Zgam ${ROOT_LIBRARIES} ${Boost_LIBRARIES}) # Линковка программы с библиотекой
else()
target_link_libraries(main Zgam ${ROOT_LIBRARIES}) # Линковка программы с библиотекой
endif()
