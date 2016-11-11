SET(Group_Project_TRIBITS_CMAKE_MINIMUM_REQUIRED 2.8.11)
CMAKE_MINIMUM_REQUIRED(VERSION ${Group_Project_TRIBITS_CMAKE_MINIMUM_REQUIRED})

#
#Define project name and set up major options
#
INCLUDE(${CMAKE_SOURCE_DIR}/ProjectName.cmake)

#Indictate to TriBITS that we might use Python
SET(${PROJECT_NAME}_USES_PYTHON_TRUE)

#Automatically read in external repositories
SET(${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE "Continuous" CACHE STRING "Enabled in top level CMakeLists.txt")

#Cmake still requires the PROJ to be defined in the top level CMakeList.txt file
PROJECT(${PROJECT_NAME} NONE)

#
#Pull in TriBITS for the project
#
SET(${PROJECT_NAME}_TRIBITS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake/tribits" CACHE PATH "" )
MESSAGE(STATUS "${PROJECT_NAME}_TRIBITS_DIR - '${${PROJECT_NAME}_TRIBITS_DIR}'")
INCLUDE(${${PROJECT_NAME}_TRIBITS_DIR}/TriBITS.cmake)

#
#Process the TriBITS Project
#
TRIBITS_PROJECT()

