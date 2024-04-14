# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /media/diehlm/Data/Linux/Stage_CY/armadillo-12.8.2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /media/diehlm/Data/Linux/Stage_CY/armadillo-12.8.2

# Include any dependencies generated for this target.
include CMakeFiles/armadillo.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/armadillo.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/armadillo.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/armadillo.dir/flags.make

CMakeFiles/armadillo.dir/src/wrapper1.cpp.o: CMakeFiles/armadillo.dir/flags.make
CMakeFiles/armadillo.dir/src/wrapper1.cpp.o: src/wrapper1.cpp
CMakeFiles/armadillo.dir/src/wrapper1.cpp.o: CMakeFiles/armadillo.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/diehlm/Data/Linux/Stage_CY/armadillo-12.8.2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/armadillo.dir/src/wrapper1.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/armadillo.dir/src/wrapper1.cpp.o -MF CMakeFiles/armadillo.dir/src/wrapper1.cpp.o.d -o CMakeFiles/armadillo.dir/src/wrapper1.cpp.o -c /media/diehlm/Data/Linux/Stage_CY/armadillo-12.8.2/src/wrapper1.cpp

CMakeFiles/armadillo.dir/src/wrapper1.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/armadillo.dir/src/wrapper1.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/diehlm/Data/Linux/Stage_CY/armadillo-12.8.2/src/wrapper1.cpp > CMakeFiles/armadillo.dir/src/wrapper1.cpp.i

CMakeFiles/armadillo.dir/src/wrapper1.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/armadillo.dir/src/wrapper1.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/diehlm/Data/Linux/Stage_CY/armadillo-12.8.2/src/wrapper1.cpp -o CMakeFiles/armadillo.dir/src/wrapper1.cpp.s

CMakeFiles/armadillo.dir/src/wrapper2.cpp.o: CMakeFiles/armadillo.dir/flags.make
CMakeFiles/armadillo.dir/src/wrapper2.cpp.o: src/wrapper2.cpp
CMakeFiles/armadillo.dir/src/wrapper2.cpp.o: CMakeFiles/armadillo.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/diehlm/Data/Linux/Stage_CY/armadillo-12.8.2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/armadillo.dir/src/wrapper2.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/armadillo.dir/src/wrapper2.cpp.o -MF CMakeFiles/armadillo.dir/src/wrapper2.cpp.o.d -o CMakeFiles/armadillo.dir/src/wrapper2.cpp.o -c /media/diehlm/Data/Linux/Stage_CY/armadillo-12.8.2/src/wrapper2.cpp

CMakeFiles/armadillo.dir/src/wrapper2.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/armadillo.dir/src/wrapper2.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/diehlm/Data/Linux/Stage_CY/armadillo-12.8.2/src/wrapper2.cpp > CMakeFiles/armadillo.dir/src/wrapper2.cpp.i

CMakeFiles/armadillo.dir/src/wrapper2.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/armadillo.dir/src/wrapper2.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/diehlm/Data/Linux/Stage_CY/armadillo-12.8.2/src/wrapper2.cpp -o CMakeFiles/armadillo.dir/src/wrapper2.cpp.s

# Object files for target armadillo
armadillo_OBJECTS = \
"CMakeFiles/armadillo.dir/src/wrapper1.cpp.o" \
"CMakeFiles/armadillo.dir/src/wrapper2.cpp.o"

# External object files for target armadillo
armadillo_EXTERNAL_OBJECTS =

libarmadillo.so.12.8.2: CMakeFiles/armadillo.dir/src/wrapper1.cpp.o
libarmadillo.so.12.8.2: CMakeFiles/armadillo.dir/src/wrapper2.cpp.o
libarmadillo.so.12.8.2: CMakeFiles/armadillo.dir/build.make
libarmadillo.so.12.8.2: /usr/lib/x86_64-linux-gnu/libopenblas.so
libarmadillo.so.12.8.2: /usr/lib/x86_64-linux-gnu/liblapack.so
libarmadillo.so.12.8.2: /usr/lib/x86_64-linux-gnu/libarpack.so
libarmadillo.so.12.8.2: /usr/lib/x86_64-linux-gnu/libsuperlu.so
libarmadillo.so.12.8.2: CMakeFiles/armadillo.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/diehlm/Data/Linux/Stage_CY/armadillo-12.8.2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX shared library libarmadillo.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/armadillo.dir/link.txt --verbose=$(VERBOSE)
	$(CMAKE_COMMAND) -E cmake_symlink_library libarmadillo.so.12.8.2 libarmadillo.so.12 libarmadillo.so

libarmadillo.so.12: libarmadillo.so.12.8.2
	@$(CMAKE_COMMAND) -E touch_nocreate libarmadillo.so.12

libarmadillo.so: libarmadillo.so.12.8.2
	@$(CMAKE_COMMAND) -E touch_nocreate libarmadillo.so

# Rule to build all files generated by this target.
CMakeFiles/armadillo.dir/build: libarmadillo.so
.PHONY : CMakeFiles/armadillo.dir/build

CMakeFiles/armadillo.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/armadillo.dir/cmake_clean.cmake
.PHONY : CMakeFiles/armadillo.dir/clean

CMakeFiles/armadillo.dir/depend:
	cd /media/diehlm/Data/Linux/Stage_CY/armadillo-12.8.2 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/diehlm/Data/Linux/Stage_CY/armadillo-12.8.2 /media/diehlm/Data/Linux/Stage_CY/armadillo-12.8.2 /media/diehlm/Data/Linux/Stage_CY/armadillo-12.8.2 /media/diehlm/Data/Linux/Stage_CY/armadillo-12.8.2 /media/diehlm/Data/Linux/Stage_CY/armadillo-12.8.2/CMakeFiles/armadillo.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/armadillo.dir/depend

