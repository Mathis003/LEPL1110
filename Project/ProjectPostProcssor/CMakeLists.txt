cmake_minimum_required(VERSION 3.6)
project(myFemProject)

set(EXEC_NAME myFem)

# OpenGL
find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR})

# GLFW
add_subdirectory(glfw)
include_directories(glfw/include)

# link the aforementioned libraries
add_executable(${EXEC_NAME} ../../Project/src/fem.c ../../ProjectPreProcessor/src/glfem.c src/main.c)
target_link_libraries(${EXEC_NAME} glfw ${OPENGL_gl_LIBRARY})

if(APPLE)
    target_compile_options(${EXEC_NAME} PRIVATE "-Wdeprecated-declarations")
    target_compile_definitions(${EXEC_NAME} PRIVATE -D GL_SILENCE_DEPRECATION)
endif()

# set exec as the startup project in visual studio
if (MSVC)
    set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${EXEC_NAME})
    target_compile_options(${EXEC_NAME} PRIVATE "/W2")
endif()

target_compile_definitions(${EXEC_NAME} PRIVATE graphic)

# add a custom target so we can do make run directly
add_custom_target(run COMMAND ${EXEC_NAME} DEPENDS ${EXEC_NAME} COMMENT "./${EXEC_NAME}")