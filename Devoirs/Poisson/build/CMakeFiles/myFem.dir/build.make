# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.25.3/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.25.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/build

# Include any dependencies generated for this target.
include CMakeFiles/myFem.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/myFem.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/myFem.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/myFem.dir/flags.make

CMakeFiles/myFem.dir/src/fem.c.o: CMakeFiles/myFem.dir/flags.make
CMakeFiles/myFem.dir/src/fem.c.o: /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/src/fem.c
CMakeFiles/myFem.dir/src/fem.c.o: CMakeFiles/myFem.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/myFem.dir/src/fem.c.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/myFem.dir/src/fem.c.o -MF CMakeFiles/myFem.dir/src/fem.c.o.d -o CMakeFiles/myFem.dir/src/fem.c.o -c /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/src/fem.c

CMakeFiles/myFem.dir/src/fem.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/myFem.dir/src/fem.c.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/src/fem.c > CMakeFiles/myFem.dir/src/fem.c.i

CMakeFiles/myFem.dir/src/fem.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/myFem.dir/src/fem.c.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/src/fem.c -o CMakeFiles/myFem.dir/src/fem.c.s

CMakeFiles/myFem.dir/src/glfem.c.o: CMakeFiles/myFem.dir/flags.make
CMakeFiles/myFem.dir/src/glfem.c.o: /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/src/glfem.c
CMakeFiles/myFem.dir/src/glfem.c.o: CMakeFiles/myFem.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/myFem.dir/src/glfem.c.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/myFem.dir/src/glfem.c.o -MF CMakeFiles/myFem.dir/src/glfem.c.o.d -o CMakeFiles/myFem.dir/src/glfem.c.o -c /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/src/glfem.c

CMakeFiles/myFem.dir/src/glfem.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/myFem.dir/src/glfem.c.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/src/glfem.c > CMakeFiles/myFem.dir/src/glfem.c.i

CMakeFiles/myFem.dir/src/glfem.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/myFem.dir/src/glfem.c.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/src/glfem.c -o CMakeFiles/myFem.dir/src/glfem.c.s

CMakeFiles/myFem.dir/src/homework.c.o: CMakeFiles/myFem.dir/flags.make
CMakeFiles/myFem.dir/src/homework.c.o: /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/src/homework.c
CMakeFiles/myFem.dir/src/homework.c.o: CMakeFiles/myFem.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/myFem.dir/src/homework.c.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/myFem.dir/src/homework.c.o -MF CMakeFiles/myFem.dir/src/homework.c.o.d -o CMakeFiles/myFem.dir/src/homework.c.o -c /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/src/homework.c

CMakeFiles/myFem.dir/src/homework.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/myFem.dir/src/homework.c.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/src/homework.c > CMakeFiles/myFem.dir/src/homework.c.i

CMakeFiles/myFem.dir/src/homework.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/myFem.dir/src/homework.c.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/src/homework.c -o CMakeFiles/myFem.dir/src/homework.c.s

CMakeFiles/myFem.dir/src/main.c.o: CMakeFiles/myFem.dir/flags.make
CMakeFiles/myFem.dir/src/main.c.o: /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/src/main.c
CMakeFiles/myFem.dir/src/main.c.o: CMakeFiles/myFem.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object CMakeFiles/myFem.dir/src/main.c.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/myFem.dir/src/main.c.o -MF CMakeFiles/myFem.dir/src/main.c.o.d -o CMakeFiles/myFem.dir/src/main.c.o -c /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/src/main.c

CMakeFiles/myFem.dir/src/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/myFem.dir/src/main.c.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/src/main.c > CMakeFiles/myFem.dir/src/main.c.i

CMakeFiles/myFem.dir/src/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/myFem.dir/src/main.c.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/src/main.c -o CMakeFiles/myFem.dir/src/main.c.s

# Object files for target myFem
myFem_OBJECTS = \
"CMakeFiles/myFem.dir/src/fem.c.o" \
"CMakeFiles/myFem.dir/src/glfem.c.o" \
"CMakeFiles/myFem.dir/src/homework.c.o" \
"CMakeFiles/myFem.dir/src/main.c.o"

# External object files for target myFem
myFem_EXTERNAL_OBJECTS =

myFem: CMakeFiles/myFem.dir/src/fem.c.o
myFem: CMakeFiles/myFem.dir/src/glfem.c.o
myFem: CMakeFiles/myFem.dir/src/homework.c.o
myFem: CMakeFiles/myFem.dir/src/main.c.o
myFem: CMakeFiles/myFem.dir/build.make
myFem: glfw/src/libglfw3.a
myFem: CMakeFiles/myFem.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking C executable myFem"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/myFem.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/myFem.dir/build: myFem
.PHONY : CMakeFiles/myFem.dir/build

CMakeFiles/myFem.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/myFem.dir/cmake_clean.cmake
.PHONY : CMakeFiles/myFem.dir/clean

CMakeFiles/myFem.dir/depend:
	cd /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/build /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/build /Users/mathisdelsart/Desktop/LEPL1110/Devoirs/Poisson/build/CMakeFiles/myFem.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/myFem.dir/depend

