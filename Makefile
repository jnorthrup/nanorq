# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.9

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Produce verbose output by default.
VERBOSE = 1

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jim/work/narq3

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jim/work/narq3

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/local/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/local/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/jim/work/narq3/CMakeFiles /home/jim/work/narq3/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/jim/work/narq3/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named nanorq_decode

# Build rule for target.
nanorq_decode: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 nanorq_decode
.PHONY : nanorq_decode

# fast build rule for target.
nanorq_decode/fast:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/build
.PHONY : nanorq_decode/fast

#=============================================================================
# Target rules for targets named nanorq_encode

# Build rule for target.
nanorq_encode: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 nanorq_encode
.PHONY : nanorq_encode

# fast build rule for target.
nanorq_encode/fast:
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/build
.PHONY : nanorq_encode/fast

#=============================================================================
# Target rules for targets named tablegen

# Build rule for target.
tablegen: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tablegen
.PHONY : tablegen

# fast build rule for target.
tablegen/fast:
	$(MAKE) -f oblas/CMakeFiles/tablegen.dir/build.make oblas/CMakeFiles/tablegen.dir/build
.PHONY : tablegen/fast

#=============================================================================
# Target rules for targets named oblas

# Build rule for target.
oblas: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 oblas
.PHONY : oblas

# fast build rule for target.
oblas/fast:
	$(MAKE) -f oblas/CMakeFiles/oblas.dir/build.make oblas/CMakeFiles/oblas.dir/build
.PHONY : oblas/fast

bitmask.o: bitmask.c.o

.PHONY : bitmask.o

# target to build an object file
bitmask.c.o:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/bitmask.c.o
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/bitmask.c.o
.PHONY : bitmask.c.o

bitmask.i: bitmask.c.i

.PHONY : bitmask.i

# target to preprocess a source file
bitmask.c.i:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/bitmask.c.i
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/bitmask.c.i
.PHONY : bitmask.c.i

bitmask.s: bitmask.c.s

.PHONY : bitmask.s

# target to generate assembly for a file
bitmask.c.s:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/bitmask.c.s
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/bitmask.c.s
.PHONY : bitmask.c.s

chooser.o: chooser.c.o

.PHONY : chooser.o

# target to build an object file
chooser.c.o:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/chooser.c.o
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/chooser.c.o
.PHONY : chooser.c.o

chooser.i: chooser.c.i

.PHONY : chooser.i

# target to preprocess a source file
chooser.c.i:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/chooser.c.i
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/chooser.c.i
.PHONY : chooser.c.i

chooser.s: chooser.c.s

.PHONY : chooser.s

# target to generate assembly for a file
chooser.c.s:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/chooser.c.s
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/chooser.c.s
.PHONY : chooser.c.s

decode.o: decode.c.o

.PHONY : decode.o

# target to build an object file
decode.c.o:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/decode.c.o
.PHONY : decode.c.o

decode.i: decode.c.i

.PHONY : decode.i

# target to preprocess a source file
decode.c.i:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/decode.c.i
.PHONY : decode.c.i

decode.s: decode.c.s

.PHONY : decode.s

# target to generate assembly for a file
decode.c.s:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/decode.c.s
.PHONY : decode.c.s

encode.o: encode.c.o

.PHONY : encode.o

# target to build an object file
encode.c.o:
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/encode.c.o
.PHONY : encode.c.o

encode.i: encode.c.i

.PHONY : encode.i

# target to preprocess a source file
encode.c.i:
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/encode.c.i
.PHONY : encode.c.i

encode.s: encode.c.s

.PHONY : encode.s

# target to generate assembly for a file
encode.c.s:
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/encode.c.s
.PHONY : encode.c.s

graph.o: graph.c.o

.PHONY : graph.o

# target to build an object file
graph.c.o:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/graph.c.o
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/graph.c.o
.PHONY : graph.c.o

graph.i: graph.c.i

.PHONY : graph.i

# target to preprocess a source file
graph.c.i:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/graph.c.i
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/graph.c.i
.PHONY : graph.c.i

graph.s: graph.c.s

.PHONY : graph.s

# target to generate assembly for a file
graph.c.s:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/graph.c.s
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/graph.c.s
.PHONY : graph.c.s

io.o: io.c.o

.PHONY : io.o

# target to build an object file
io.c.o:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/io.c.o
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/io.c.o
.PHONY : io.c.o

io.i: io.c.i

.PHONY : io.i

# target to preprocess a source file
io.c.i:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/io.c.i
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/io.c.i
.PHONY : io.c.i

io.s: io.c.s

.PHONY : io.s

# target to generate assembly for a file
io.c.s:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/io.c.s
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/io.c.s
.PHONY : io.c.s

nanorq.o: nanorq.c.o

.PHONY : nanorq.o

# target to build an object file
nanorq.c.o:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/nanorq.c.o
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/nanorq.c.o
.PHONY : nanorq.c.o

nanorq.i: nanorq.c.i

.PHONY : nanorq.i

# target to preprocess a source file
nanorq.c.i:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/nanorq.c.i
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/nanorq.c.i
.PHONY : nanorq.c.i

nanorq.s: nanorq.c.s

.PHONY : nanorq.s

# target to generate assembly for a file
nanorq.c.s:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/nanorq.c.s
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/nanorq.c.s
.PHONY : nanorq.c.s

params.o: params.c.o

.PHONY : params.o

# target to build an object file
params.c.o:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/params.c.o
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/params.c.o
.PHONY : params.c.o

params.i: params.c.i

.PHONY : params.i

# target to preprocess a source file
params.c.i:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/params.c.i
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/params.c.i
.PHONY : params.c.i

params.s: params.c.s

.PHONY : params.s

# target to generate assembly for a file
params.c.s:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/params.c.s
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/params.c.s
.PHONY : params.c.s

precode.o: precode.c.o

.PHONY : precode.o

# target to build an object file
precode.c.o:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/precode.c.o
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/precode.c.o
.PHONY : precode.c.o

precode.i: precode.c.i

.PHONY : precode.i

# target to preprocess a source file
precode.c.i:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/precode.c.i
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/precode.c.i
.PHONY : precode.c.i

precode.s: precode.c.s

.PHONY : precode.s

# target to generate assembly for a file
precode.c.s:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/precode.c.s
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/precode.c.s
.PHONY : precode.c.s

rand.o: rand.c.o

.PHONY : rand.o

# target to build an object file
rand.c.o:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/rand.c.o
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/rand.c.o
.PHONY : rand.c.o

rand.i: rand.c.i

.PHONY : rand.i

# target to preprocess a source file
rand.c.i:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/rand.c.i
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/rand.c.i
.PHONY : rand.c.i

rand.s: rand.c.s

.PHONY : rand.s

# target to generate assembly for a file
rand.c.s:
	$(MAKE) -f CMakeFiles/nanorq_decode.dir/build.make CMakeFiles/nanorq_decode.dir/rand.c.s
	$(MAKE) -f CMakeFiles/nanorq_encode.dir/build.make CMakeFiles/nanorq_encode.dir/rand.c.s
.PHONY : rand.c.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... nanorq_decode"
	@echo "... edit_cache"
	@echo "... nanorq_encode"
	@echo "... tablegen"
	@echo "... oblas"
	@echo "... bitmask.o"
	@echo "... bitmask.i"
	@echo "... bitmask.s"
	@echo "... chooser.o"
	@echo "... chooser.i"
	@echo "... chooser.s"
	@echo "... decode.o"
	@echo "... decode.i"
	@echo "... decode.s"
	@echo "... encode.o"
	@echo "... encode.i"
	@echo "... encode.s"
	@echo "... graph.o"
	@echo "... graph.i"
	@echo "... graph.s"
	@echo "... io.o"
	@echo "... io.i"
	@echo "... io.s"
	@echo "... nanorq.o"
	@echo "... nanorq.i"
	@echo "... nanorq.s"
	@echo "... params.o"
	@echo "... params.i"
	@echo "... params.s"
	@echo "... precode.o"
	@echo "... precode.i"
	@echo "... precode.s"
	@echo "... rand.o"
	@echo "... rand.i"
	@echo "... rand.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

