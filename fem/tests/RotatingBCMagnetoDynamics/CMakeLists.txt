INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( case.sif case.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO extrude_quarter.eg mortar_quarter.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(RotatingBCMagnetoDynamics)
