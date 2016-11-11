# Gives access to SET_DEFAULT_AND_FROM_ENV function
# This is file defined in TriBITS. Like a header file in C.
INCLUDE(SetDefaultAndFromEnv)

# Settings for CTest/CDash
if(NOT DEFINED CTEST_DROP_METHOD)
  SET_DEFAULT_AND_FROM_ENV(CTEST_DROP_METHOD "https")
ENDIF()

IF(CTEST_DROP_METHOD STREQUAL "https")
  #The normal default is ${HOST_TYPE}-${COMPILER_VERSION}-${BUILD_DIR}
  #we are over-riding it by using set.
  SET(CTEST_BUILD_NAME "Linux-Lab3-itssahil")

  #To learn more about this function look at the TriBITS documentation
  SET_DEFAULT_AND_FROM_ENV(CTEST_PROJECT_NAME "Lab3")
  SET_DEFAULT_AND_FROM_ENV(CTEST_TRIGGER_SITE "")

  #CDash server hostname
  SET_DEFAULT_AND_FROM_ENV(CTEST_DROP_SITE "cdash-ners590.aura.arc-ts.umich.edu")

  #The rest of the web address
  SET_DEFAULT_AND_FROM_ENV(CTEST_DROP_LOCATION "/submit.php?project=Lab3")

  #YES SUBMIT RESULT! I GET MY 10 POINTS!
  SET_DEFAULT_AND_FROM_ENV(CTEST_DROP_SITE_CDASH TRUE)
ENDIF()

