cmake_minimum_required(VERSION 3.5)

project(TCPolygonMatching)

set(CMAKE_CXX_STANDARD 17)

option(CXX "enable C++ compilation" ON)

if(CXX)
  enable_language(CXX)
endif()

# Manual override: force HiGHS even if Gurobi is installed
option(FORCE_HIGHS "Force HiGHS backend even if Gurobi is found" OFF)

#activate address sanitizing
if(CMAKE_BUILD_TYPE MATCHES Debug)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=address")
endif()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Visual Studio compiler with static runtime libraries
if(MSVC AND MT)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd")
endif()

find_package(CGAL 6.0.1 EXACT REQUIRED)
set(CGAL_BUILD_TYPE Release)
message(STATUS "Found CGAL version ${CGAL_VERSION}")

find_package(Boost 1.83.0 EXACT REQUIRED COMPONENTS system thread)
if(Boost_FOUND)
  message(STATUS "Found Boost version ${Boost_VERSION}")
else()
  message(FATAL_ERROR "Boost 1.83 not found")
endif()

find_package(Threads REQUIRED)

find_package(GDAL REQUIRED)
if(GDAL_FOUND)
  message(STATUS "Found GDAL version ${GDAL_VERSION}")
else()
  message(FATAL_ERROR "GDAL not found")
endif()

# ---- Gurobi ----
set(USE_GUROBI OFF)
if(NOT FORCE_HIGHS)
  find_package(GUROBI QUIET)
  if(GUROBI_FOUND)
    set(USE_GUROBI ON)
    message(STATUS "Found Gurobi version ${GUROBI_VERSION} (using Gurobi backend)")
  else()
    message(STATUS "Gurobi not found (using HiGHS backend)")
  endif()
else()
  message(STATUS "FORCE_HIGHS=ON (using HiGHS backend)")
endif()

# ignore BOOST deprecated headers
add_definitions("-DBOOST_ALLOW_DEPRECATED_HEADERS")
add_definitions("-DBOOST_BIND_GLOBAL_PLACEHOLDERS")

# do not show 'copy' warning on build
add_compile_options(-Wno-nonnull)

# list source files here
set(sources
        src/main.cpp
        src/shapefile_io_operations.cpp
        src/localization.cpp
        src/graph_computations.cpp
        src/solution.cpp
        include/cgal_includes.h
        include/shapefile_io_operations.h
        include/localization.h
        include/graph_computations.h
        include/solution.h
        include/linear_program.h
        include/command_line_parser.h
        src/command_line_parser.cpp
        include/tree_computations.h
        src/tree_computations.cpp
        include/logger.h
        src/logger.cpp
        include/binary_io.h
        src/binary_io.cpp
        include/polygon_operations.h
        src/polygon_operations.cpp
        include/polygon_wh.h
        src/polygon_wh.cpp
        include/threading.h
        src/threading.cpp
        include/solver_context.h
        src/solver_context_singleton.cpp
        include/union_find.h
        src/cc_decomposition.cpp
        include/cc_decomposition.h
)

add_executable(${CMAKE_PROJECT_NAME} ${sources})

# Select backend implementation for linear programs
if(USE_GUROBI)
  target_sources(${CMAKE_PROJECT_NAME} PRIVATE include/solver_context_gurobi_access.h src/linear_program_gurobi.cpp src/solver_context_gurobi.cpp)
  target_compile_definitions(${CMAKE_PROJECT_NAME} PRIVATE HAVE_GUROBI=1)
else()
  target_sources(${CMAKE_PROJECT_NAME} PRIVATE src/linear_program_highs.cpp src/solver_context_highs.cpp)
  target_compile_definitions(${CMAKE_PROJECT_NAME} PRIVATE HAVE_HIGHS=1)
endif()

# ---- Linking ----

if(CXX)
  set(CMAKE_CXX_STANDARD 11)
  if(USE_GUROBI)
    target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE
            optimized ${GUROBI_CXX_LIBRARY}
            debug ${GUROBI_CXX_DEBUG_LIBRARY}
    )
  endif()
endif()

if(USE_GUROBI)
  target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE ${GUROBI_INCLUDE_DIRS})
  target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE ${GUROBI_LIBRARY})
endif()

target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE CGAL::CGAL)
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE Threads::Threads)
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE Boost::system Boost::thread)

if(NOT USE_GUROBI)
  add_subdirectory(external/HiGHS)
  target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE highs)
endif()

target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE ${GDAL_INCLUDE_DIRS})
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE ${GDAL_LIBRARIES})

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  include(FeatureSummary)
  feature_summary(WHAT ALL)
endif()
