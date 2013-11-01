CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(HrlAnalysis)

# Add the unitTest Directory
ADD_SUBDIRECTORY(test)

SET(CMAKE_MODULE_PATH "${HRLANALYSIS_SOURCE_DIR}/CMakeModules")

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/HrlAnalysis/bin/)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/HrlAnalysis/lib)

CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
                "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
                IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
                  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")

# Make sure the Boost Library is installed
FIND_PACKAGE(Boost COMPONENTS python REQUIRED)
FIND_PACKAGE(Boost COMPONENTS serialization REQUIRED)
FIND_PACKAGE(PythonInterp 2.4)
FIND_PACKAGE(PythonLibs 2.4)

IF(PYTHONINTERP_FOUND)
    MESSAGE (STATUS "Found Python Interpreter ${PYTHON_EXECUTABLE}")
ELSE(PYTHONINTERP_FOUND)
    MESSAGE (STATUS "Could not find Python Interpreter!")
ENDIF(PYTHONINTERP_FOUND)

IF(PYTHONLIBS_FOUND)
    MESSAGE (STATUS "Found Python Libs -- ${PYTHON_LIBRARIES} -- ${PYTHON_INCLUDE_PATH} -- ${PYTHON_DEBUG_LIBRARIES}")
ELSE(PYTHONLIBS_FOUND)
    MESSAGE (STATUS "Could not find Python Libs!")
ENDIF(PYTHONLIBS_FOUND)

IF(Boost_FOUND)

    AUX_SOURCE_DIRECTORY("stlPythonWrapper" stlPythonWrapper)
    AUX_SOURCE_DIRECTORY("HrlNetworkAnalysis" HrlNetworkAnalysis)
    AUX_SOURCE_DIRECTORY("HrlNeuralAnalysis" HrlNeuralAnalysis)
        
    INCLUDE_DIRECTORIES("stlPythonWrapper"
                        "HrlNeuralAnalysis"
                        "HrlNetworkAnalysis"
                        ${Boost_INCLUDE_DIRS}
                        ${PYTHON_INCLUDE_PATH}
                        ".")

    ADD_LIBRARY(HrlAnalysis SHARED "HrlAnalysis.cpp"
                           ${stlPythonWrapper}
                           ${HrlNeuralAnalysis}
                           ${HrlNetworkAnalysis})

    SET(CFLAGS " -fPIC ")
    SET(CMAKE_BUILD_TYPE Release)

    # Set the libraries that need to be linked
    TARGET_LINK_LIBRARIES(HrlAnalysis stdc++ ${PYTHON_LIBRARIES} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
    SET_TARGET_PROPERTIES(HrlAnalysis PROPERTIES LINKER_LANGUAGE CXX)

ELSE(Boost_FOUND)
    MESSAGE(WARNING, "The Boost Library was not found.")
ENDIF(Boost_FOUND)

