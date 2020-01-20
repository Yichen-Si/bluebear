cmake_minimum_required (VERSION 2.6)
project(ibs0phase)

include_directories(${CMAKE_SOURCE_DIR}/../htslib/include)

execute_process(COMMAND date OUTPUT_VARIABLE DATE OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND whoami OUTPUT_VARIABLE USER OUTPUT_STRIP_TRAILING_WHITESPACE)
set(VERSION "1.0.0")

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pipe -D__STDC_LIMIT_MACROS -Wall -Wno-unused-local-typedefs -Wno-enum-compare -fpic -O2 -DVERSION='\"${VERSION}\"' -DUSER='\"${USER}\"' -DDATE='\"${DATE}\"' -D__ZLIB_AVAILABLE__")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_subdirectory(src)

set(ibs0phase_SOURCES src/Error.cpp
    src/commands.cpp
    src/params.cpp
    src/hts_utils.cpp
    src/filter.cpp
    src/bcf_ordered_reader.cpp
    src/bcf_ordered_writer.cpp
    src/genome_interval.cpp
    src/genotype_concordance.cpp
    src/interval.cpp
    src/interval_tree.cpp
    src/reference_sequence.cpp
    src/tsv_reader.cpp
    src/cmd_vcf_annotate_context.cpp
    src/utils.cpp
    src/compact_matrix.cpp
    src/bp2cm.cpp
    src/cmd_vcf_ibs0_baseline.cpp
    src/cmd_vcf_ibs0.cpp
    src/cmd_vcf_ibs0_flanking.cpp
    src/cmd_vcf_ibs0_pairwise.cpp
    src/cmd_vcf_ibs0_unconditional.cpp
    src/cmd_vcf_ibs0_view.cpp
    src/ibs0_pairwise_scan.cpp
    src/cmd_vcf_raresharing.cpp
    src/cmd_vcf_sfs.cpp
    src/pbwt_prefix.cpp
    src/pbwt_suffix.cpp
    src/trio_phase.cpp
    src/pt2interval.cpp
    src/ibs0.cpp
    src/ibs0_phase_backward.cpp
    src/ibs0_phase_forward.cpp
    src/hap_ibd_pbwt.cpp
    src/ibd_around_pt.cpp
    src/test_hts.cpp
    src/ibs0rare.cpp
    src/ibs0_rare_phase_backward.cpp
    src/ibs0_rare_phase_forward.cpp
    src/ibs0_rare_only_phase_backward.cpp
    src/ibs0_rare_only_phase_forward.cpp
    src/ibs0_blip_forward.cpp
    src/kmer_count.cpp
    src/cmd_vcf_triallelic.cpp
    src/cmd_ibs0_add_control.cpp
    src/ibs0_between_group_control.cpp
    src/cmd_vcf_annotate_ibs0_small_n.cpp
    src/bluebear.cpp)

set(ibs0phase_HEADERS src/cramore.h
    src/Error.h
    src/pException.h
    src/commands.h
    src/params.h
    src/bcf_filter_arg.h
    src/bcf_ordered_reader.h
    src/bcf_ordered_writer.h
    src/bcf_variant_key.h
    src/hts_utils.h
    src/gtf_interval_tree.h
    src/interval.h
    src/interval_tree.h
    src/filter.h
    src/config.h
    src/Constant.h
    src/genome_interval.h
    src/genomeLoci.h
    src/genomeMask.h
    src/genomeScore.h
    src/genotype_concordance.h
    src/reference_sequence.h
    src/tsv_reader.h
    src/utils.h
    src/cm2bp.h
    src/compact_matrix.h
    src/pbwt_build.h
    src/bp2cm.h
    src/ibs0.h
    src/rare_variant_ibs0.h)

add_executable(bluebear ${ibs0phase_SOURCES} ${ibs0phase_HEADERS})
target_link_libraries(bluebear ${CMAKE_SOURCE_DIR}/../htslib/libhts.a dl pthread lzma z bz2 curl crypto)

