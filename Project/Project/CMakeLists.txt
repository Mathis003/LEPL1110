cmake_minimum_required(VERSION 3.6)
project(myFemProject)

set(EXEC_NAME myFem)

# link the aforementioned libraries
add_executable(${EXEC_NAME} src/fem.c  src/homework.c src/main.c)

# set exec as the startup project in visual studio
if (MSVC)
    set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${EXEC_NAME})
    target_compile_options(${EXEC_NAME} PRIVATE "/W2")
endif()
if (NOT WIN32)
    find_library(MATH_LIBRARY m REQUIRED)
    target_link_libraries(${EXEC_NAME} PUBLIC ${MATH_LIBRARY})
endif()

# add a custom target so we can do make run directly
add_custom_target(run COMMAND ${EXEC_NAME} DEPENDS ${EXEC_NAME} COMMENT "./${EXEC_NAME}")