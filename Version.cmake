#
# Single file that needs to be chnaged on when branching
#
IF("${REPOSITORY_NAME}" STREQUAL "")
	#Need to set the repository name
	SET(REPOSITORY_NAME itssahil_Group_Project)
ENDIF()

#Version Info
SET(${REPOSITORY_NAME}_VERSION 0.0.0)
SET(${REPOSITORY_NAME}_MAJOR_VERSION 00)
SET(${REPOSITORY_NAME}_MAJOR_MINOR_VERSION 000000)
SET(${REPOSITORY_NAME}_VERSION_STRING "0.0.0 (Dev)")

