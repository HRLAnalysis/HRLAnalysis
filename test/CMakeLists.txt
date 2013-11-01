ADD_SUBDIRECTORY(lib)

AUX_SOURCE_DIRECTORY("../HrlNeuralAnalysis" HrlNeuralAnalysis)
AUX_SOURCE_DIRECTORY("SpikeAnalysisAccess" SpikeAnalysisAccess)
AUX_SOURCE_DIRECTORY("AnalysisTest" AnalysisTest)
AUX_SOURCE_DIRECTORY("AnalysisData" AnalysisData)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/HrlAnalysis/test)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/HrlAnalysis/lib)

# Look for cppunit
INCLUDE(FindPkgConfig)
PKG_CHECK_MODULES(CPPUNIT REQUIRED cppunit )

MESSAGE( "CPPUNIT_INCLUDE_DIRS : ${CPPUNIT_INCLUDE_DIRS}" )
MESSAGE( "CPPUNIT_LIBRARY_DIRS : ${CPPUNIT_LIBRARY_DIRS}" )
MESSAGE( "CPPUNIT_LIBRARIES : ${CPPUNIT_LIBRARIES}" )

FIND_PACKAGE(Boost COMPONENTS serialization REQUIRED)
FIND_PACKAGE(Boost REQUIRED)

IF(Boost_FOUND)
    MESSAGE (STATUS "Found Boost ${Boost_VERSION} -- ${Boost_LIBRARIES} -- ${Boost_INCLUDE_DIRS}")
ENDIF(Boost_FOUND)       

IF( CPPUNIT_FOUND  )
    INCLUDE_DIRECTORIES("../HrlNeuralAnalysis"
                        "SpikeAnalysisAccess"
                        "AnalysisTest"
                        "AnalysisData"
                        "."
                        ${CPPUNIT_LIBRARIES}
                        ${Boost_INCLUDE_DIRS})

    ADD_EXECUTABLE(hrlAnalysisTest  "TestSuite.cpp"
    								${SpikeAnalysisAccess}
                                    ${HrlNeuralAnalysis}
                                    ${AnalysisTest}
                                    ${AnalysisData})

    #SET(CSTANDARD "-std=c++0x")
    #SET(CFLAGS " ${CSTANDARD} ")
    #SET(CMAKE_CXX_FLAGS "-std=c++0x")
    
    SET(CFLAGS " -fPIC ")
    SET(CMAKE_BUILD_TYPE Release)
    
    # Set the libraries that need to be linked
    TARGET_LINK_LIBRARIES(hrlAnalysisTest stdc++ cppunit ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
    SET_TARGET_PROPERTIES(hrlAnalysisTest PROPERTIES LINKER_LANGUAGE CXX)

ENDIF(CPPUNIT_FOUND)