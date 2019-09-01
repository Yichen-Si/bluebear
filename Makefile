# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.8

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /net/wonderland/home/ycsi/tool/bluebear

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /net/wonderland/home/ycsi/tool/bluebear

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /net/wonderland/home/ycsi/tool/bluebear/CMakeFiles /net/wonderland/home/ycsi/tool/bluebear/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /net/wonderland/home/ycsi/tool/bluebear/CMakeFiles 0
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
# Target rules for targets named bluebear

# Build rule for target.
bluebear: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 bluebear
.PHONY : bluebear

# fast build rule for target.
bluebear/fast:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/build
.PHONY : bluebear/fast

src/Error.o: src/Error.cpp.o

.PHONY : src/Error.o

# target to build an object file
src/Error.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/Error.cpp.o
.PHONY : src/Error.cpp.o

src/Error.i: src/Error.cpp.i

.PHONY : src/Error.i

# target to preprocess a source file
src/Error.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/Error.cpp.i
.PHONY : src/Error.cpp.i

src/Error.s: src/Error.cpp.s

.PHONY : src/Error.s

# target to generate assembly for a file
src/Error.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/Error.cpp.s
.PHONY : src/Error.cpp.s

src/bcf_ordered_reader.o: src/bcf_ordered_reader.cpp.o

.PHONY : src/bcf_ordered_reader.o

# target to build an object file
src/bcf_ordered_reader.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/bcf_ordered_reader.cpp.o
.PHONY : src/bcf_ordered_reader.cpp.o

src/bcf_ordered_reader.i: src/bcf_ordered_reader.cpp.i

.PHONY : src/bcf_ordered_reader.i

# target to preprocess a source file
src/bcf_ordered_reader.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/bcf_ordered_reader.cpp.i
.PHONY : src/bcf_ordered_reader.cpp.i

src/bcf_ordered_reader.s: src/bcf_ordered_reader.cpp.s

.PHONY : src/bcf_ordered_reader.s

# target to generate assembly for a file
src/bcf_ordered_reader.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/bcf_ordered_reader.cpp.s
.PHONY : src/bcf_ordered_reader.cpp.s

src/bcf_ordered_writer.o: src/bcf_ordered_writer.cpp.o

.PHONY : src/bcf_ordered_writer.o

# target to build an object file
src/bcf_ordered_writer.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/bcf_ordered_writer.cpp.o
.PHONY : src/bcf_ordered_writer.cpp.o

src/bcf_ordered_writer.i: src/bcf_ordered_writer.cpp.i

.PHONY : src/bcf_ordered_writer.i

# target to preprocess a source file
src/bcf_ordered_writer.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/bcf_ordered_writer.cpp.i
.PHONY : src/bcf_ordered_writer.cpp.i

src/bcf_ordered_writer.s: src/bcf_ordered_writer.cpp.s

.PHONY : src/bcf_ordered_writer.s

# target to generate assembly for a file
src/bcf_ordered_writer.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/bcf_ordered_writer.cpp.s
.PHONY : src/bcf_ordered_writer.cpp.s

src/bluebear.o: src/bluebear.cpp.o

.PHONY : src/bluebear.o

# target to build an object file
src/bluebear.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/bluebear.cpp.o
.PHONY : src/bluebear.cpp.o

src/bluebear.i: src/bluebear.cpp.i

.PHONY : src/bluebear.i

# target to preprocess a source file
src/bluebear.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/bluebear.cpp.i
.PHONY : src/bluebear.cpp.i

src/bluebear.s: src/bluebear.cpp.s

.PHONY : src/bluebear.s

# target to generate assembly for a file
src/bluebear.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/bluebear.cpp.s
.PHONY : src/bluebear.cpp.s

src/cmd_vcf_ibs0.o: src/cmd_vcf_ibs0.cpp.o

.PHONY : src/cmd_vcf_ibs0.o

# target to build an object file
src/cmd_vcf_ibs0.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_ibs0.cpp.o
.PHONY : src/cmd_vcf_ibs0.cpp.o

src/cmd_vcf_ibs0.i: src/cmd_vcf_ibs0.cpp.i

.PHONY : src/cmd_vcf_ibs0.i

# target to preprocess a source file
src/cmd_vcf_ibs0.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_ibs0.cpp.i
.PHONY : src/cmd_vcf_ibs0.cpp.i

src/cmd_vcf_ibs0.s: src/cmd_vcf_ibs0.cpp.s

.PHONY : src/cmd_vcf_ibs0.s

# target to generate assembly for a file
src/cmd_vcf_ibs0.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_ibs0.cpp.s
.PHONY : src/cmd_vcf_ibs0.cpp.s

src/cmd_vcf_ibs0_baseline.o: src/cmd_vcf_ibs0_baseline.cpp.o

.PHONY : src/cmd_vcf_ibs0_baseline.o

# target to build an object file
src/cmd_vcf_ibs0_baseline.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_ibs0_baseline.cpp.o
.PHONY : src/cmd_vcf_ibs0_baseline.cpp.o

src/cmd_vcf_ibs0_baseline.i: src/cmd_vcf_ibs0_baseline.cpp.i

.PHONY : src/cmd_vcf_ibs0_baseline.i

# target to preprocess a source file
src/cmd_vcf_ibs0_baseline.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_ibs0_baseline.cpp.i
.PHONY : src/cmd_vcf_ibs0_baseline.cpp.i

src/cmd_vcf_ibs0_baseline.s: src/cmd_vcf_ibs0_baseline.cpp.s

.PHONY : src/cmd_vcf_ibs0_baseline.s

# target to generate assembly for a file
src/cmd_vcf_ibs0_baseline.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_ibs0_baseline.cpp.s
.PHONY : src/cmd_vcf_ibs0_baseline.cpp.s

src/cmd_vcf_ibs0_flanking.o: src/cmd_vcf_ibs0_flanking.cpp.o

.PHONY : src/cmd_vcf_ibs0_flanking.o

# target to build an object file
src/cmd_vcf_ibs0_flanking.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_ibs0_flanking.cpp.o
.PHONY : src/cmd_vcf_ibs0_flanking.cpp.o

src/cmd_vcf_ibs0_flanking.i: src/cmd_vcf_ibs0_flanking.cpp.i

.PHONY : src/cmd_vcf_ibs0_flanking.i

# target to preprocess a source file
src/cmd_vcf_ibs0_flanking.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_ibs0_flanking.cpp.i
.PHONY : src/cmd_vcf_ibs0_flanking.cpp.i

src/cmd_vcf_ibs0_flanking.s: src/cmd_vcf_ibs0_flanking.cpp.s

.PHONY : src/cmd_vcf_ibs0_flanking.s

# target to generate assembly for a file
src/cmd_vcf_ibs0_flanking.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_ibs0_flanking.cpp.s
.PHONY : src/cmd_vcf_ibs0_flanking.cpp.s

src/cmd_vcf_ibs0_pairwise.o: src/cmd_vcf_ibs0_pairwise.cpp.o

.PHONY : src/cmd_vcf_ibs0_pairwise.o

# target to build an object file
src/cmd_vcf_ibs0_pairwise.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_ibs0_pairwise.cpp.o
.PHONY : src/cmd_vcf_ibs0_pairwise.cpp.o

src/cmd_vcf_ibs0_pairwise.i: src/cmd_vcf_ibs0_pairwise.cpp.i

.PHONY : src/cmd_vcf_ibs0_pairwise.i

# target to preprocess a source file
src/cmd_vcf_ibs0_pairwise.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_ibs0_pairwise.cpp.i
.PHONY : src/cmd_vcf_ibs0_pairwise.cpp.i

src/cmd_vcf_ibs0_pairwise.s: src/cmd_vcf_ibs0_pairwise.cpp.s

.PHONY : src/cmd_vcf_ibs0_pairwise.s

# target to generate assembly for a file
src/cmd_vcf_ibs0_pairwise.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_ibs0_pairwise.cpp.s
.PHONY : src/cmd_vcf_ibs0_pairwise.cpp.s

src/cmd_vcf_ibs0_unconditional.o: src/cmd_vcf_ibs0_unconditional.cpp.o

.PHONY : src/cmd_vcf_ibs0_unconditional.o

# target to build an object file
src/cmd_vcf_ibs0_unconditional.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_ibs0_unconditional.cpp.o
.PHONY : src/cmd_vcf_ibs0_unconditional.cpp.o

src/cmd_vcf_ibs0_unconditional.i: src/cmd_vcf_ibs0_unconditional.cpp.i

.PHONY : src/cmd_vcf_ibs0_unconditional.i

# target to preprocess a source file
src/cmd_vcf_ibs0_unconditional.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_ibs0_unconditional.cpp.i
.PHONY : src/cmd_vcf_ibs0_unconditional.cpp.i

src/cmd_vcf_ibs0_unconditional.s: src/cmd_vcf_ibs0_unconditional.cpp.s

.PHONY : src/cmd_vcf_ibs0_unconditional.s

# target to generate assembly for a file
src/cmd_vcf_ibs0_unconditional.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_ibs0_unconditional.cpp.s
.PHONY : src/cmd_vcf_ibs0_unconditional.cpp.s

src/cmd_vcf_ibs0_view.o: src/cmd_vcf_ibs0_view.cpp.o

.PHONY : src/cmd_vcf_ibs0_view.o

# target to build an object file
src/cmd_vcf_ibs0_view.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_ibs0_view.cpp.o
.PHONY : src/cmd_vcf_ibs0_view.cpp.o

src/cmd_vcf_ibs0_view.i: src/cmd_vcf_ibs0_view.cpp.i

.PHONY : src/cmd_vcf_ibs0_view.i

# target to preprocess a source file
src/cmd_vcf_ibs0_view.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_ibs0_view.cpp.i
.PHONY : src/cmd_vcf_ibs0_view.cpp.i

src/cmd_vcf_ibs0_view.s: src/cmd_vcf_ibs0_view.cpp.s

.PHONY : src/cmd_vcf_ibs0_view.s

# target to generate assembly for a file
src/cmd_vcf_ibs0_view.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_ibs0_view.cpp.s
.PHONY : src/cmd_vcf_ibs0_view.cpp.s

src/cmd_vcf_raresharing.o: src/cmd_vcf_raresharing.cpp.o

.PHONY : src/cmd_vcf_raresharing.o

# target to build an object file
src/cmd_vcf_raresharing.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_raresharing.cpp.o
.PHONY : src/cmd_vcf_raresharing.cpp.o

src/cmd_vcf_raresharing.i: src/cmd_vcf_raresharing.cpp.i

.PHONY : src/cmd_vcf_raresharing.i

# target to preprocess a source file
src/cmd_vcf_raresharing.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_raresharing.cpp.i
.PHONY : src/cmd_vcf_raresharing.cpp.i

src/cmd_vcf_raresharing.s: src/cmd_vcf_raresharing.cpp.s

.PHONY : src/cmd_vcf_raresharing.s

# target to generate assembly for a file
src/cmd_vcf_raresharing.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_raresharing.cpp.s
.PHONY : src/cmd_vcf_raresharing.cpp.s

src/cmd_vcf_sfs.o: src/cmd_vcf_sfs.cpp.o

.PHONY : src/cmd_vcf_sfs.o

# target to build an object file
src/cmd_vcf_sfs.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_sfs.cpp.o
.PHONY : src/cmd_vcf_sfs.cpp.o

src/cmd_vcf_sfs.i: src/cmd_vcf_sfs.cpp.i

.PHONY : src/cmd_vcf_sfs.i

# target to preprocess a source file
src/cmd_vcf_sfs.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_sfs.cpp.i
.PHONY : src/cmd_vcf_sfs.cpp.i

src/cmd_vcf_sfs.s: src/cmd_vcf_sfs.cpp.s

.PHONY : src/cmd_vcf_sfs.s

# target to generate assembly for a file
src/cmd_vcf_sfs.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/cmd_vcf_sfs.cpp.s
.PHONY : src/cmd_vcf_sfs.cpp.s

src/commands.o: src/commands.cpp.o

.PHONY : src/commands.o

# target to build an object file
src/commands.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/commands.cpp.o
.PHONY : src/commands.cpp.o

src/commands.i: src/commands.cpp.i

.PHONY : src/commands.i

# target to preprocess a source file
src/commands.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/commands.cpp.i
.PHONY : src/commands.cpp.i

src/commands.s: src/commands.cpp.s

.PHONY : src/commands.s

# target to generate assembly for a file
src/commands.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/commands.cpp.s
.PHONY : src/commands.cpp.s

src/compact_matrix.o: src/compact_matrix.cpp.o

.PHONY : src/compact_matrix.o

# target to build an object file
src/compact_matrix.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/compact_matrix.cpp.o
.PHONY : src/compact_matrix.cpp.o

src/compact_matrix.i: src/compact_matrix.cpp.i

.PHONY : src/compact_matrix.i

# target to preprocess a source file
src/compact_matrix.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/compact_matrix.cpp.i
.PHONY : src/compact_matrix.cpp.i

src/compact_matrix.s: src/compact_matrix.cpp.s

.PHONY : src/compact_matrix.s

# target to generate assembly for a file
src/compact_matrix.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/compact_matrix.cpp.s
.PHONY : src/compact_matrix.cpp.s

src/filter.o: src/filter.cpp.o

.PHONY : src/filter.o

# target to build an object file
src/filter.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/filter.cpp.o
.PHONY : src/filter.cpp.o

src/filter.i: src/filter.cpp.i

.PHONY : src/filter.i

# target to preprocess a source file
src/filter.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/filter.cpp.i
.PHONY : src/filter.cpp.i

src/filter.s: src/filter.cpp.s

.PHONY : src/filter.s

# target to generate assembly for a file
src/filter.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/filter.cpp.s
.PHONY : src/filter.cpp.s

src/genome_interval.o: src/genome_interval.cpp.o

.PHONY : src/genome_interval.o

# target to build an object file
src/genome_interval.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/genome_interval.cpp.o
.PHONY : src/genome_interval.cpp.o

src/genome_interval.i: src/genome_interval.cpp.i

.PHONY : src/genome_interval.i

# target to preprocess a source file
src/genome_interval.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/genome_interval.cpp.i
.PHONY : src/genome_interval.cpp.i

src/genome_interval.s: src/genome_interval.cpp.s

.PHONY : src/genome_interval.s

# target to generate assembly for a file
src/genome_interval.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/genome_interval.cpp.s
.PHONY : src/genome_interval.cpp.s

src/genotype_concordance.o: src/genotype_concordance.cpp.o

.PHONY : src/genotype_concordance.o

# target to build an object file
src/genotype_concordance.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/genotype_concordance.cpp.o
.PHONY : src/genotype_concordance.cpp.o

src/genotype_concordance.i: src/genotype_concordance.cpp.i

.PHONY : src/genotype_concordance.i

# target to preprocess a source file
src/genotype_concordance.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/genotype_concordance.cpp.i
.PHONY : src/genotype_concordance.cpp.i

src/genotype_concordance.s: src/genotype_concordance.cpp.s

.PHONY : src/genotype_concordance.s

# target to generate assembly for a file
src/genotype_concordance.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/genotype_concordance.cpp.s
.PHONY : src/genotype_concordance.cpp.s

src/hap_ibd_pbwt.o: src/hap_ibd_pbwt.cpp.o

.PHONY : src/hap_ibd_pbwt.o

# target to build an object file
src/hap_ibd_pbwt.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/hap_ibd_pbwt.cpp.o
.PHONY : src/hap_ibd_pbwt.cpp.o

src/hap_ibd_pbwt.i: src/hap_ibd_pbwt.cpp.i

.PHONY : src/hap_ibd_pbwt.i

# target to preprocess a source file
src/hap_ibd_pbwt.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/hap_ibd_pbwt.cpp.i
.PHONY : src/hap_ibd_pbwt.cpp.i

src/hap_ibd_pbwt.s: src/hap_ibd_pbwt.cpp.s

.PHONY : src/hap_ibd_pbwt.s

# target to generate assembly for a file
src/hap_ibd_pbwt.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/hap_ibd_pbwt.cpp.s
.PHONY : src/hap_ibd_pbwt.cpp.s

src/hts_utils.o: src/hts_utils.cpp.o

.PHONY : src/hts_utils.o

# target to build an object file
src/hts_utils.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/hts_utils.cpp.o
.PHONY : src/hts_utils.cpp.o

src/hts_utils.i: src/hts_utils.cpp.i

.PHONY : src/hts_utils.i

# target to preprocess a source file
src/hts_utils.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/hts_utils.cpp.i
.PHONY : src/hts_utils.cpp.i

src/hts_utils.s: src/hts_utils.cpp.s

.PHONY : src/hts_utils.s

# target to generate assembly for a file
src/hts_utils.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/hts_utils.cpp.s
.PHONY : src/hts_utils.cpp.s

src/ibd_around_pt.o: src/ibd_around_pt.cpp.o

.PHONY : src/ibd_around_pt.o

# target to build an object file
src/ibd_around_pt.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/ibd_around_pt.cpp.o
.PHONY : src/ibd_around_pt.cpp.o

src/ibd_around_pt.i: src/ibd_around_pt.cpp.i

.PHONY : src/ibd_around_pt.i

# target to preprocess a source file
src/ibd_around_pt.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/ibd_around_pt.cpp.i
.PHONY : src/ibd_around_pt.cpp.i

src/ibd_around_pt.s: src/ibd_around_pt.cpp.s

.PHONY : src/ibd_around_pt.s

# target to generate assembly for a file
src/ibd_around_pt.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/ibd_around_pt.cpp.s
.PHONY : src/ibd_around_pt.cpp.s

src/ibs0.o: src/ibs0.cpp.o

.PHONY : src/ibs0.o

# target to build an object file
src/ibs0.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/ibs0.cpp.o
.PHONY : src/ibs0.cpp.o

src/ibs0.i: src/ibs0.cpp.i

.PHONY : src/ibs0.i

# target to preprocess a source file
src/ibs0.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/ibs0.cpp.i
.PHONY : src/ibs0.cpp.i

src/ibs0.s: src/ibs0.cpp.s

.PHONY : src/ibs0.s

# target to generate assembly for a file
src/ibs0.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/ibs0.cpp.s
.PHONY : src/ibs0.cpp.s

src/ibs0_phase_backward.o: src/ibs0_phase_backward.cpp.o

.PHONY : src/ibs0_phase_backward.o

# target to build an object file
src/ibs0_phase_backward.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/ibs0_phase_backward.cpp.o
.PHONY : src/ibs0_phase_backward.cpp.o

src/ibs0_phase_backward.i: src/ibs0_phase_backward.cpp.i

.PHONY : src/ibs0_phase_backward.i

# target to preprocess a source file
src/ibs0_phase_backward.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/ibs0_phase_backward.cpp.i
.PHONY : src/ibs0_phase_backward.cpp.i

src/ibs0_phase_backward.s: src/ibs0_phase_backward.cpp.s

.PHONY : src/ibs0_phase_backward.s

# target to generate assembly for a file
src/ibs0_phase_backward.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/ibs0_phase_backward.cpp.s
.PHONY : src/ibs0_phase_backward.cpp.s

src/ibs0_phase_forward.o: src/ibs0_phase_forward.cpp.o

.PHONY : src/ibs0_phase_forward.o

# target to build an object file
src/ibs0_phase_forward.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/ibs0_phase_forward.cpp.o
.PHONY : src/ibs0_phase_forward.cpp.o

src/ibs0_phase_forward.i: src/ibs0_phase_forward.cpp.i

.PHONY : src/ibs0_phase_forward.i

# target to preprocess a source file
src/ibs0_phase_forward.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/ibs0_phase_forward.cpp.i
.PHONY : src/ibs0_phase_forward.cpp.i

src/ibs0_phase_forward.s: src/ibs0_phase_forward.cpp.s

.PHONY : src/ibs0_phase_forward.s

# target to generate assembly for a file
src/ibs0_phase_forward.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/ibs0_phase_forward.cpp.s
.PHONY : src/ibs0_phase_forward.cpp.s

src/interval.o: src/interval.cpp.o

.PHONY : src/interval.o

# target to build an object file
src/interval.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/interval.cpp.o
.PHONY : src/interval.cpp.o

src/interval.i: src/interval.cpp.i

.PHONY : src/interval.i

# target to preprocess a source file
src/interval.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/interval.cpp.i
.PHONY : src/interval.cpp.i

src/interval.s: src/interval.cpp.s

.PHONY : src/interval.s

# target to generate assembly for a file
src/interval.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/interval.cpp.s
.PHONY : src/interval.cpp.s

src/interval_tree.o: src/interval_tree.cpp.o

.PHONY : src/interval_tree.o

# target to build an object file
src/interval_tree.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/interval_tree.cpp.o
.PHONY : src/interval_tree.cpp.o

src/interval_tree.i: src/interval_tree.cpp.i

.PHONY : src/interval_tree.i

# target to preprocess a source file
src/interval_tree.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/interval_tree.cpp.i
.PHONY : src/interval_tree.cpp.i

src/interval_tree.s: src/interval_tree.cpp.s

.PHONY : src/interval_tree.s

# target to generate assembly for a file
src/interval_tree.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/interval_tree.cpp.s
.PHONY : src/interval_tree.cpp.s

src/params.o: src/params.cpp.o

.PHONY : src/params.o

# target to build an object file
src/params.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/params.cpp.o
.PHONY : src/params.cpp.o

src/params.i: src/params.cpp.i

.PHONY : src/params.i

# target to preprocess a source file
src/params.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/params.cpp.i
.PHONY : src/params.cpp.i

src/params.s: src/params.cpp.s

.PHONY : src/params.s

# target to generate assembly for a file
src/params.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/params.cpp.s
.PHONY : src/params.cpp.s

src/pbwt_prefix.o: src/pbwt_prefix.cpp.o

.PHONY : src/pbwt_prefix.o

# target to build an object file
src/pbwt_prefix.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/pbwt_prefix.cpp.o
.PHONY : src/pbwt_prefix.cpp.o

src/pbwt_prefix.i: src/pbwt_prefix.cpp.i

.PHONY : src/pbwt_prefix.i

# target to preprocess a source file
src/pbwt_prefix.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/pbwt_prefix.cpp.i
.PHONY : src/pbwt_prefix.cpp.i

src/pbwt_prefix.s: src/pbwt_prefix.cpp.s

.PHONY : src/pbwt_prefix.s

# target to generate assembly for a file
src/pbwt_prefix.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/pbwt_prefix.cpp.s
.PHONY : src/pbwt_prefix.cpp.s

src/pbwt_suffix.o: src/pbwt_suffix.cpp.o

.PHONY : src/pbwt_suffix.o

# target to build an object file
src/pbwt_suffix.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/pbwt_suffix.cpp.o
.PHONY : src/pbwt_suffix.cpp.o

src/pbwt_suffix.i: src/pbwt_suffix.cpp.i

.PHONY : src/pbwt_suffix.i

# target to preprocess a source file
src/pbwt_suffix.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/pbwt_suffix.cpp.i
.PHONY : src/pbwt_suffix.cpp.i

src/pbwt_suffix.s: src/pbwt_suffix.cpp.s

.PHONY : src/pbwt_suffix.s

# target to generate assembly for a file
src/pbwt_suffix.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/pbwt_suffix.cpp.s
.PHONY : src/pbwt_suffix.cpp.s

src/pt2interval.o: src/pt2interval.cpp.o

.PHONY : src/pt2interval.o

# target to build an object file
src/pt2interval.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/pt2interval.cpp.o
.PHONY : src/pt2interval.cpp.o

src/pt2interval.i: src/pt2interval.cpp.i

.PHONY : src/pt2interval.i

# target to preprocess a source file
src/pt2interval.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/pt2interval.cpp.i
.PHONY : src/pt2interval.cpp.i

src/pt2interval.s: src/pt2interval.cpp.s

.PHONY : src/pt2interval.s

# target to generate assembly for a file
src/pt2interval.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/pt2interval.cpp.s
.PHONY : src/pt2interval.cpp.s

src/reference_sequence.o: src/reference_sequence.cpp.o

.PHONY : src/reference_sequence.o

# target to build an object file
src/reference_sequence.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/reference_sequence.cpp.o
.PHONY : src/reference_sequence.cpp.o

src/reference_sequence.i: src/reference_sequence.cpp.i

.PHONY : src/reference_sequence.i

# target to preprocess a source file
src/reference_sequence.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/reference_sequence.cpp.i
.PHONY : src/reference_sequence.cpp.i

src/reference_sequence.s: src/reference_sequence.cpp.s

.PHONY : src/reference_sequence.s

# target to generate assembly for a file
src/reference_sequence.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/reference_sequence.cpp.s
.PHONY : src/reference_sequence.cpp.s

src/test_hts.o: src/test_hts.cpp.o

.PHONY : src/test_hts.o

# target to build an object file
src/test_hts.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/test_hts.cpp.o
.PHONY : src/test_hts.cpp.o

src/test_hts.i: src/test_hts.cpp.i

.PHONY : src/test_hts.i

# target to preprocess a source file
src/test_hts.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/test_hts.cpp.i
.PHONY : src/test_hts.cpp.i

src/test_hts.s: src/test_hts.cpp.s

.PHONY : src/test_hts.s

# target to generate assembly for a file
src/test_hts.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/test_hts.cpp.s
.PHONY : src/test_hts.cpp.s

src/trio_phase.o: src/trio_phase.cpp.o

.PHONY : src/trio_phase.o

# target to build an object file
src/trio_phase.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/trio_phase.cpp.o
.PHONY : src/trio_phase.cpp.o

src/trio_phase.i: src/trio_phase.cpp.i

.PHONY : src/trio_phase.i

# target to preprocess a source file
src/trio_phase.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/trio_phase.cpp.i
.PHONY : src/trio_phase.cpp.i

src/trio_phase.s: src/trio_phase.cpp.s

.PHONY : src/trio_phase.s

# target to generate assembly for a file
src/trio_phase.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/trio_phase.cpp.s
.PHONY : src/trio_phase.cpp.s

src/tsv_reader.o: src/tsv_reader.cpp.o

.PHONY : src/tsv_reader.o

# target to build an object file
src/tsv_reader.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/tsv_reader.cpp.o
.PHONY : src/tsv_reader.cpp.o

src/tsv_reader.i: src/tsv_reader.cpp.i

.PHONY : src/tsv_reader.i

# target to preprocess a source file
src/tsv_reader.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/tsv_reader.cpp.i
.PHONY : src/tsv_reader.cpp.i

src/tsv_reader.s: src/tsv_reader.cpp.s

.PHONY : src/tsv_reader.s

# target to generate assembly for a file
src/tsv_reader.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/tsv_reader.cpp.s
.PHONY : src/tsv_reader.cpp.s

src/utils.o: src/utils.cpp.o

.PHONY : src/utils.o

# target to build an object file
src/utils.cpp.o:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/utils.cpp.o
.PHONY : src/utils.cpp.o

src/utils.i: src/utils.cpp.i

.PHONY : src/utils.i

# target to preprocess a source file
src/utils.cpp.i:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/utils.cpp.i
.PHONY : src/utils.cpp.i

src/utils.s: src/utils.cpp.s

.PHONY : src/utils.s

# target to generate assembly for a file
src/utils.cpp.s:
	$(MAKE) -f CMakeFiles/bluebear.dir/build.make CMakeFiles/bluebear.dir/src/utils.cpp.s
.PHONY : src/utils.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... edit_cache"
	@echo "... bluebear"
	@echo "... src/Error.o"
	@echo "... src/Error.i"
	@echo "... src/Error.s"
	@echo "... src/bcf_ordered_reader.o"
	@echo "... src/bcf_ordered_reader.i"
	@echo "... src/bcf_ordered_reader.s"
	@echo "... src/bcf_ordered_writer.o"
	@echo "... src/bcf_ordered_writer.i"
	@echo "... src/bcf_ordered_writer.s"
	@echo "... src/bluebear.o"
	@echo "... src/bluebear.i"
	@echo "... src/bluebear.s"
	@echo "... src/cmd_vcf_ibs0.o"
	@echo "... src/cmd_vcf_ibs0.i"
	@echo "... src/cmd_vcf_ibs0.s"
	@echo "... src/cmd_vcf_ibs0_baseline.o"
	@echo "... src/cmd_vcf_ibs0_baseline.i"
	@echo "... src/cmd_vcf_ibs0_baseline.s"
	@echo "... src/cmd_vcf_ibs0_flanking.o"
	@echo "... src/cmd_vcf_ibs0_flanking.i"
	@echo "... src/cmd_vcf_ibs0_flanking.s"
	@echo "... src/cmd_vcf_ibs0_pairwise.o"
	@echo "... src/cmd_vcf_ibs0_pairwise.i"
	@echo "... src/cmd_vcf_ibs0_pairwise.s"
	@echo "... src/cmd_vcf_ibs0_unconditional.o"
	@echo "... src/cmd_vcf_ibs0_unconditional.i"
	@echo "... src/cmd_vcf_ibs0_unconditional.s"
	@echo "... src/cmd_vcf_ibs0_view.o"
	@echo "... src/cmd_vcf_ibs0_view.i"
	@echo "... src/cmd_vcf_ibs0_view.s"
	@echo "... src/cmd_vcf_raresharing.o"
	@echo "... src/cmd_vcf_raresharing.i"
	@echo "... src/cmd_vcf_raresharing.s"
	@echo "... src/cmd_vcf_sfs.o"
	@echo "... src/cmd_vcf_sfs.i"
	@echo "... src/cmd_vcf_sfs.s"
	@echo "... src/commands.o"
	@echo "... src/commands.i"
	@echo "... src/commands.s"
	@echo "... src/compact_matrix.o"
	@echo "... src/compact_matrix.i"
	@echo "... src/compact_matrix.s"
	@echo "... src/filter.o"
	@echo "... src/filter.i"
	@echo "... src/filter.s"
	@echo "... src/genome_interval.o"
	@echo "... src/genome_interval.i"
	@echo "... src/genome_interval.s"
	@echo "... src/genotype_concordance.o"
	@echo "... src/genotype_concordance.i"
	@echo "... src/genotype_concordance.s"
	@echo "... src/hap_ibd_pbwt.o"
	@echo "... src/hap_ibd_pbwt.i"
	@echo "... src/hap_ibd_pbwt.s"
	@echo "... src/hts_utils.o"
	@echo "... src/hts_utils.i"
	@echo "... src/hts_utils.s"
	@echo "... src/ibd_around_pt.o"
	@echo "... src/ibd_around_pt.i"
	@echo "... src/ibd_around_pt.s"
	@echo "... src/ibs0.o"
	@echo "... src/ibs0.i"
	@echo "... src/ibs0.s"
	@echo "... src/ibs0_phase_backward.o"
	@echo "... src/ibs0_phase_backward.i"
	@echo "... src/ibs0_phase_backward.s"
	@echo "... src/ibs0_phase_forward.o"
	@echo "... src/ibs0_phase_forward.i"
	@echo "... src/ibs0_phase_forward.s"
	@echo "... src/interval.o"
	@echo "... src/interval.i"
	@echo "... src/interval.s"
	@echo "... src/interval_tree.o"
	@echo "... src/interval_tree.i"
	@echo "... src/interval_tree.s"
	@echo "... src/params.o"
	@echo "... src/params.i"
	@echo "... src/params.s"
	@echo "... src/pbwt_prefix.o"
	@echo "... src/pbwt_prefix.i"
	@echo "... src/pbwt_prefix.s"
	@echo "... src/pbwt_suffix.o"
	@echo "... src/pbwt_suffix.i"
	@echo "... src/pbwt_suffix.s"
	@echo "... src/pt2interval.o"
	@echo "... src/pt2interval.i"
	@echo "... src/pt2interval.s"
	@echo "... src/reference_sequence.o"
	@echo "... src/reference_sequence.i"
	@echo "... src/reference_sequence.s"
	@echo "... src/test_hts.o"
	@echo "... src/test_hts.i"
	@echo "... src/test_hts.s"
	@echo "... src/trio_phase.o"
	@echo "... src/trio_phase.i"
	@echo "... src/trio_phase.s"
	@echo "... src/tsv_reader.o"
	@echo "... src/tsv_reader.i"
	@echo "... src/tsv_reader.s"
	@echo "... src/utils.o"
	@echo "... src/utils.i"
	@echo "... src/utils.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

