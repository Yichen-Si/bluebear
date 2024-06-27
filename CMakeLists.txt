cmake_minimum_required (VERSION 3.12)
project(bluebear)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-result -O3")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

find_path(HTS_INCLUDE_DIRS htslib/bgzf.h HINTS ${CMAKE_SOURCE_DIR}/../htslib)
if(NOT HTS_INCLUDE_DIRS )
    message(FATAL_ERROR "libhts HTS_INCLUDE_DIRS not found")
endif()
message(STATUS ${HTS_INCLUDE_DIRS})

find_library(HTS_LIBRARIES libhts.a HINTS ${CMAKE_SOURCE_DIR}/../htslib)
if(NOT HTS_LIBRARIES)
    message(FATAL_ERROR "libhts HTS_LIBRARIES not found")
endif()
message(STATUS ${HTS_LIBRARIES})

find_path(QGEN_INCLUDE_DIRS qgenlib/qgen_error.h  HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ ${CMAKE_SOURCE_DIR}/../qgenlib)
if(NOT QGEN_INCLUDE_DIRS )
    message(FATAL_ERROR "libqgen QGEN_INCLUDE_DIRS not found")
endif()
message(STATUS ${QGEN_INCLUDE_DIRS})

find_library(QGEN_LIBRARIES libqgen.a  HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ ${CMAKE_SOURCE_DIR}/../qgenlib/lib/)
if(NOT QGEN_LIBRARIES)
    message(FATAL_ERROR "libqgen QGEN_LIBRARIES not found")
endif()
message(STATUS ${QGEN_LIBRARIES})

set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(${HTS_INCLUDE_DIRS}
					${QGEN_INCLUDE_DIRS}
                    ${CMAKE_CURRENT_SOURCE_DIR}/src
					${CMAKE_CURRENT_SOURCE_DIR}/src/ibs0 )

add_executable(bluebear
    src/cramore.h
    src/Error.h
	src/Error.cpp
    src/pException.h
    src/commands.h
    src/commands.cpp
    src/params.h
    src/params.cpp
	src/genomeChunk.h
	src/genomeChunk.cpp
	src/PhredHelper.h
	src/PhredHelper.cpp
	src/bcf_filter_arg.h
	src/bcf_filtered_reader.h
    src/bcf_filtered_reader.cpp
	src/bcf_chunked_reader.h
	src/bcf_chunked_reader.cpp
    src/bcf_ordered_reader.h
    src/bcf_ordered_reader.cpp
    src/bcf_ordered_writer.h
    src/bcf_ordered_writer.cpp
    src/bcf_variant_key.h
    src/hts_utils.h
    src/hts_utils.cpp
    src/gtf_interval_tree.h
    src/interval.h
	src/interval.cpp
    src/interval_tree.h
	src/interval_tree.cpp
    src/filter.h
    src/filter.cpp
    src/config.h
    src/Constant.h
    src/genome_interval.h
	src/genome_interval.cpp
    src/genomeLoci.h
    src/genomeMask.h
    src/genomeScore.h
    src/genotype_concordance.h
	src/genotype_concordance.cpp
    src/reference_sequence.h
	src/reference_sequence.cpp
    src/tsv_reader.h
	src/tsv_reader.cpp
    src/utils.h
    src/brent.hpp
	src/brent.cpp
    src/brent_obj.h
    src/fa_reader.h
    src/seq_basics.h

	src/byte_operation.h
	src/byte_operation.cpp
	src/compact_matrix.h
	src/compact_matrix.cpp
	src/haploid_hash.h
	src/binary_variant_compress.h
	src/binary_variant_compress.cpp
    src/pbwt_build.h

    src/cm2bp.h
    src/bp2cm.h
	src/bp2cm.cpp
    src/ibs0.h
	src/ibs0.cpp
	src/ibs0rare.cpp
    src/rare_variant_ibs0.h
    src/rare_variant_config.h
	src/rare_variant_config.cpp

    src/bluebear.cpp
    src/test_hts.cpp
    src/cmd_vcf_temp.cpp

    src/utils.cpp

	src/cmd_fa_cpg.cpp
	src/cmd_dge_shuffle.cpp
	src/cmd_tsv_annotate_context.cpp
    src/cmd_vcf_annotate_context.cpp
    src/cmd_vcf_sfs.cpp
    src/cmd_vcf_sfs_mr.cpp
    src/kmer_count.cpp
    src/cmd_vcf_triallelic.cpp
    src/pt2interval.cpp
    src/pbwt_prefix.cpp
    src/pbwt_suffix.cpp
    src/hap_ibd_pbwt.cpp
    src/ibd_around_pt.cpp
    src/trio_phase.cpp
	src/cmd_vcf_phase_outlier.cpp
    src/cmd_vcf_trio_denovo_mutation.cpp
    src/cmd_vcf_annotate_ibs0.cpp
    src/cmd_vcf_raresharing.cpp
	src/cmd_vcf_rare_allele_carrier.cpp
	src/cmd_vcf_shared_allele_count.cpp
	src/cmd_vcf_collapse_variant.cpp
    src/cmd_vcf_insert_pm.cpp
    src/cmd_vcf_info_cdf.cpp
    src/cmd_vcf_site_annotate_refsample_info.cpp
    src/cmd_vcf_sample_summary.cpp
	src/cmd_vcf_sample_summary_window.cpp
    src/cmd_vcf_site_annotate_pm_posterior.cpp
    src/cmd_vcf_site_annotate_age.cpp
    src/cmd_vcf_annotate_ibd.cpp
	src/cmd_vcf_annotate_ibd_small.cpp
    src/cmd_vcf_anno_haplo_group.cpp
    src/cmd_vcf_snppairingene.cpp
    src/cmd_vcf_snppairinfo.cpp
	src/cmd_vcf_md.cpp
	src/cmd_vcf_dprime.cpp
	src/blocks_in_ibs0_info.cpp

	src/ibs0/cmd_vcf_ibs0_approx_pm.cpp

    # src/ibs0/ibs0_flanking.cpp
    # src/ibs0/ibs0_pairwise.cpp
	# src/ibs0/ibs0_pairwise_scan.cpp
    # src/ibs0/ibs0_view.cpp
    # src/ibs0/ibs0_short.cpp
    # src/ibs0/ibs0_triallelic.cpp
	# src/ibs0/ibs0_rare_phase_backward.cpp
    # src/ibs0/ibs0_rare_phase_forward.cpp
    # src/ibs0/ibs0_rare_only_phase_backward.cpp
    # src/ibs0/ibs0_rare_only_phase_forward.cpp
    # src/ibs0/ibs0_blip_forward.cpp
    # src/ibs0/ibs0_add_control.cpp
    # src/ibs0/ibs0_between_group_control.cpp

    Eigen/Cholesky
	Eigen/CholmodSupport
	Eigen/Core
	Eigen/Dense
	Eigen/Eigen
	Eigen/Eigenvalues
	Eigen/Geometry
	Eigen/Householder
	Eigen/IterativeLinearSolvers
	Eigen/Jacobi
	Eigen/LU
	Eigen/MetisSupport
	Eigen/OrderingMethods
	Eigen/PaStiXSupport
	Eigen/PardisoSupport
	Eigen/QR
	Eigen/QtAlignedMalloc
	Eigen/SPQRSupport
	Eigen/SVD
	Eigen/Sparse
	Eigen/SparseCholesky
	Eigen/SparseCore
	Eigen/SparseLU
	Eigen/SparseQR
	Eigen/StdDeque
	Eigen/StdList
	Eigen/StdVector
	Eigen/SuperLUSupport
	Eigen/UmfPackSupport
	Eigen/src/Cholesky/LDLT.h
	Eigen/src/Cholesky/LLT.h
	Eigen/src/Cholesky/LLT_LAPACKE.h
	Eigen/src/CholmodSupport/CholmodSupport.h
	Eigen/src/Core/Array.h
	Eigen/src/Core/ArrayBase.h
	Eigen/src/Core/ArrayWrapper.h
	Eigen/src/Core/Assign.h
	Eigen/src/Core/AssignEvaluator.h
	Eigen/src/Core/Assign_MKL.h
	Eigen/src/Core/BandMatrix.h
	Eigen/src/Core/Block.h
	Eigen/src/Core/BooleanRedux.h
	Eigen/src/Core/CommaInitializer.h
	Eigen/src/Core/ConditionEstimator.h
	Eigen/src/Core/CoreEvaluators.h
	Eigen/src/Core/CoreIterators.h
	Eigen/src/Core/CwiseBinaryOp.h
	Eigen/src/Core/CwiseNullaryOp.h
	Eigen/src/Core/CwiseTernaryOp.h
	Eigen/src/Core/CwiseUnaryOp.h
	Eigen/src/Core/CwiseUnaryView.h
	Eigen/src/Core/DenseBase.h
	Eigen/src/Core/DenseCoeffsBase.h
	Eigen/src/Core/DenseStorage.h
	Eigen/src/Core/Diagonal.h
	Eigen/src/Core/DiagonalMatrix.h
	Eigen/src/Core/DiagonalProduct.h
	Eigen/src/Core/Dot.h
	Eigen/src/Core/EigenBase.h
	Eigen/src/Core/ForceAlignedAccess.h
	Eigen/src/Core/Fuzzy.h
	Eigen/src/Core/GeneralProduct.h
	Eigen/src/Core/GenericPacketMath.h
	Eigen/src/Core/GlobalFunctions.h
	Eigen/src/Core/IO.h
	Eigen/src/Core/Inverse.h
	Eigen/src/Core/Map.h
	Eigen/src/Core/MapBase.h
	Eigen/src/Core/MathFunctions.h
	Eigen/src/Core/MathFunctionsImpl.h
	Eigen/src/Core/Matrix.h
	Eigen/src/Core/MatrixBase.h
	Eigen/src/Core/NestByValue.h
	Eigen/src/Core/NoAlias.h
	Eigen/src/Core/NumTraits.h
	Eigen/src/Core/PermutationMatrix.h
	Eigen/src/Core/PlainObjectBase.h
	Eigen/src/Core/Product.h
	Eigen/src/Core/ProductEvaluators.h
	Eigen/src/Core/Random.h
	Eigen/src/Core/Redux.h
	Eigen/src/Core/Ref.h
	Eigen/src/Core/Replicate.h
	Eigen/src/Core/ReturnByValue.h
	Eigen/src/Core/Reverse.h
	Eigen/src/Core/Select.h
	Eigen/src/Core/SelfAdjointView.h
	Eigen/src/Core/SelfCwiseBinaryOp.h
	Eigen/src/Core/Solve.h
	Eigen/src/Core/SolveTriangular.h
	Eigen/src/Core/SolverBase.h
	Eigen/src/Core/StableNorm.h
	Eigen/src/Core/Stride.h
	Eigen/src/Core/Swap.h
	Eigen/src/Core/Transpose.h
	Eigen/src/Core/Transpositions.h
	Eigen/src/Core/TriangularMatrix.h
	Eigen/src/Core/VectorBlock.h
	Eigen/src/Core/VectorwiseOp.h
	Eigen/src/Core/Visitor.h
	Eigen/src/Core/arch/AVX/Complex.h
	Eigen/src/Core/arch/AVX/MathFunctions.h
	Eigen/src/Core/arch/AVX/PacketMath.h
	Eigen/src/Core/arch/AVX/TypeCasting.h
	Eigen/src/Core/arch/AVX512/MathFunctions.h
	Eigen/src/Core/arch/AVX512/PacketMath.h
	Eigen/src/Core/arch/AltiVec/Complex.h
	Eigen/src/Core/arch/AltiVec/MathFunctions.h
	Eigen/src/Core/arch/AltiVec/PacketMath.h
	Eigen/src/Core/arch/CUDA/Complex.h
	Eigen/src/Core/arch/CUDA/Half.h
	Eigen/src/Core/arch/CUDA/MathFunctions.h
	Eigen/src/Core/arch/CUDA/PacketMath.h
	Eigen/src/Core/arch/CUDA/PacketMathHalf.h
	Eigen/src/Core/arch/CUDA/TypeCasting.h
	Eigen/src/Core/arch/Default/Settings.h
	Eigen/src/Core/arch/NEON/Complex.h
	Eigen/src/Core/arch/NEON/MathFunctions.h
	Eigen/src/Core/arch/NEON/PacketMath.h
	Eigen/src/Core/arch/SSE/Complex.h
	Eigen/src/Core/arch/SSE/MathFunctions.h
	Eigen/src/Core/arch/SSE/PacketMath.h
	Eigen/src/Core/arch/SSE/TypeCasting.h
	Eigen/src/Core/arch/ZVector/Complex.h
	Eigen/src/Core/arch/ZVector/MathFunctions.h
	Eigen/src/Core/arch/ZVector/PacketMath.h
	Eigen/src/Core/functors/AssignmentFunctors.h
	Eigen/src/Core/functors/BinaryFunctors.h
	Eigen/src/Core/functors/NullaryFunctors.h
	Eigen/src/Core/functors/StlFunctors.h
	Eigen/src/Core/functors/TernaryFunctors.h
	Eigen/src/Core/functors/UnaryFunctors.h
	Eigen/src/Core/products/GeneralBlockPanelKernel.h
	Eigen/src/Core/products/GeneralMatrixMatrix.h
	Eigen/src/Core/products/GeneralMatrixMatrixTriangular.h
	Eigen/src/Core/products/GeneralMatrixMatrixTriangular_BLAS.h
	Eigen/src/Core/products/GeneralMatrixMatrix_BLAS.h
	Eigen/src/Core/products/GeneralMatrixVector.h
	Eigen/src/Core/products/GeneralMatrixVector_BLAS.h
	Eigen/src/Core/products/Parallelizer.h
	Eigen/src/Core/products/SelfadjointMatrixMatrix.h
	Eigen/src/Core/products/SelfadjointMatrixMatrix_BLAS.h
	Eigen/src/Core/products/SelfadjointMatrixVector.h
	Eigen/src/Core/products/SelfadjointMatrixVector_BLAS.h
	Eigen/src/Core/products/SelfadjointProduct.h
	Eigen/src/Core/products/SelfadjointRank2Update.h
	Eigen/src/Core/products/TriangularMatrixMatrix.h
	Eigen/src/Core/products/TriangularMatrixMatrix_BLAS.h
	Eigen/src/Core/products/TriangularMatrixVector.h
	Eigen/src/Core/products/TriangularMatrixVector_BLAS.h
	Eigen/src/Core/products/TriangularSolverMatrix.h
	Eigen/src/Core/products/TriangularSolverMatrix_BLAS.h
	Eigen/src/Core/products/TriangularSolverVector.h
	Eigen/src/Core/util/BlasUtil.h
	Eigen/src/Core/util/Constants.h
	Eigen/src/Core/util/DisableStupidWarnings.h
	Eigen/src/Core/util/ForwardDeclarations.h
	Eigen/src/Core/util/MKL_support.h
	Eigen/src/Core/util/Macros.h
	Eigen/src/Core/util/Memory.h
	Eigen/src/Core/util/Meta.h
	Eigen/src/Core/util/NonMPL2.h
	Eigen/src/Core/util/ReenableStupidWarnings.h
	Eigen/src/Core/util/StaticAssert.h
	Eigen/src/Core/util/XprHelper.h
	Eigen/src/Eigenvalues/ComplexEigenSolver.h
	Eigen/src/Eigenvalues/ComplexSchur.h
	Eigen/src/Eigenvalues/ComplexSchur_LAPACKE.h
	Eigen/src/Eigenvalues/EigenSolver.h
	Eigen/src/Eigenvalues/GeneralizedEigenSolver.h
	Eigen/src/Eigenvalues/GeneralizedSelfAdjointEigenSolver.h
	Eigen/src/Eigenvalues/HessenbergDecomposition.h
	Eigen/src/Eigenvalues/MatrixBaseEigenvalues.h
	Eigen/src/Eigenvalues/RealQZ.h
	Eigen/src/Eigenvalues/RealSchur.h
	Eigen/src/Eigenvalues/RealSchur_LAPACKE.h
	Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h
	Eigen/src/Eigenvalues/SelfAdjointEigenSolver_LAPACKE.h
	Eigen/src/Eigenvalues/Tridiagonalization.h
	Eigen/src/Geometry/AlignedBox.h
	Eigen/src/Geometry/AngleAxis.h
	Eigen/src/Geometry/EulerAngles.h
	Eigen/src/Geometry/Homogeneous.h
	Eigen/src/Geometry/Hyperplane.h
	Eigen/src/Geometry/OrthoMethods.h
	Eigen/src/Geometry/ParametrizedLine.h
	Eigen/src/Geometry/Quaternion.h
	Eigen/src/Geometry/Rotation2D.h
	Eigen/src/Geometry/RotationBase.h
	Eigen/src/Geometry/Scaling.h
	Eigen/src/Geometry/Transform.h
	Eigen/src/Geometry/Translation.h
	Eigen/src/Geometry/Umeyama.h
	Eigen/src/Geometry/arch/Geometry_SSE.h
	Eigen/src/Householder/BlockHouseholder.h
	Eigen/src/Householder/Householder.h
	Eigen/src/Householder/HouseholderSequence.h
	Eigen/src/IterativeLinearSolvers/BasicPreconditioners.h
	Eigen/src/IterativeLinearSolvers/BiCGSTAB.h
	Eigen/src/IterativeLinearSolvers/ConjugateGradient.h
	Eigen/src/IterativeLinearSolvers/IncompleteCholesky.h
	Eigen/src/IterativeLinearSolvers/IncompleteLUT.h
	Eigen/src/IterativeLinearSolvers/IterativeSolverBase.h
	Eigen/src/IterativeLinearSolvers/LeastSquareConjugateGradient.h
	Eigen/src/IterativeLinearSolvers/SolveWithGuess.h
	Eigen/src/Jacobi/Jacobi.h
	Eigen/src/LU/Determinant.h
	Eigen/src/LU/FullPivLU.h
	Eigen/src/LU/InverseImpl.h
	Eigen/src/LU/PartialPivLU.h
	Eigen/src/LU/PartialPivLU_LAPACKE.h
	Eigen/src/LU/arch/Inverse_SSE.h
	Eigen/src/MetisSupport/MetisSupport.h
	Eigen/src/OrderingMethods/Amd.h
	Eigen/src/OrderingMethods/Eigen_Colamd.h
	Eigen/src/OrderingMethods/Ordering.h
	Eigen/src/PaStiXSupport/PaStiXSupport.h
	Eigen/src/PardisoSupport/PardisoSupport.h
	Eigen/src/QR/ColPivHouseholderQR.h
	Eigen/src/QR/ColPivHouseholderQR_LAPACKE.h
	Eigen/src/QR/CompleteOrthogonalDecomposition.h
	Eigen/src/QR/FullPivHouseholderQR.h
	Eigen/src/QR/HouseholderQR.h
	Eigen/src/QR/HouseholderQR_LAPACKE.h
	Eigen/src/SPQRSupport/SuiteSparseQRSupport.h
	Eigen/src/SVD/BDCSVD.h
	Eigen/src/SVD/JacobiSVD.h
	Eigen/src/SVD/JacobiSVD_LAPACKE.h
	Eigen/src/SVD/SVDBase.h
	Eigen/src/SVD/UpperBidiagonalization.h
	Eigen/src/SparseCholesky/SimplicialCholesky.h
	Eigen/src/SparseCholesky/SimplicialCholesky_impl.h
	Eigen/src/SparseCore/AmbiVector.h
	Eigen/src/SparseCore/CompressedStorage.h
	Eigen/src/SparseCore/ConservativeSparseSparseProduct.h
	Eigen/src/SparseCore/MappedSparseMatrix.h
	Eigen/src/SparseCore/SparseAssign.h
	Eigen/src/SparseCore/SparseBlock.h
	Eigen/src/SparseCore/SparseColEtree.h
	Eigen/src/SparseCore/SparseCompressedBase.h
	Eigen/src/SparseCore/SparseCwiseBinaryOp.h
	Eigen/src/SparseCore/SparseCwiseUnaryOp.h
	Eigen/src/SparseCore/SparseDenseProduct.h
	Eigen/src/SparseCore/SparseDiagonalProduct.h
	Eigen/src/SparseCore/SparseDot.h
	Eigen/src/SparseCore/SparseFuzzy.h
	Eigen/src/SparseCore/SparseMap.h
	Eigen/src/SparseCore/SparseMatrix.h
	Eigen/src/SparseCore/SparseMatrixBase.h
	Eigen/src/SparseCore/SparsePermutation.h
	Eigen/src/SparseCore/SparseProduct.h
	Eigen/src/SparseCore/SparseRedux.h
	Eigen/src/SparseCore/SparseRef.h
	Eigen/src/SparseCore/SparseSelfAdjointView.h
	Eigen/src/SparseCore/SparseSolverBase.h
	Eigen/src/SparseCore/SparseSparseProductWithPruning.h
	Eigen/src/SparseCore/SparseTranspose.h
	Eigen/src/SparseCore/SparseTriangularView.h
	Eigen/src/SparseCore/SparseUtil.h
	Eigen/src/SparseCore/SparseVector.h
	Eigen/src/SparseCore/SparseView.h
	Eigen/src/SparseCore/TriangularSolver.h
	Eigen/src/SparseLU/SparseLU.h
	Eigen/src/SparseLU/SparseLUImpl.h
	Eigen/src/SparseLU/SparseLU_Memory.h
	Eigen/src/SparseLU/SparseLU_Structs.h
	Eigen/src/SparseLU/SparseLU_SupernodalMatrix.h
	Eigen/src/SparseLU/SparseLU_Utils.h
	Eigen/src/SparseLU/SparseLU_column_bmod.h
	Eigen/src/SparseLU/SparseLU_column_dfs.h
	Eigen/src/SparseLU/SparseLU_copy_to_ucol.h
	Eigen/src/SparseLU/SparseLU_gemm_kernel.h
	Eigen/src/SparseLU/SparseLU_heap_relax_snode.h
	Eigen/src/SparseLU/SparseLU_kernel_bmod.h
	Eigen/src/SparseLU/SparseLU_panel_bmod.h
	Eigen/src/SparseLU/SparseLU_panel_dfs.h
	Eigen/src/SparseLU/SparseLU_pivotL.h
	Eigen/src/SparseLU/SparseLU_pruneL.h
	Eigen/src/SparseLU/SparseLU_relax_snode.h
	Eigen/src/SparseQR/SparseQR.h
	Eigen/src/StlSupport/StdDeque.h
	Eigen/src/StlSupport/StdList.h
	Eigen/src/StlSupport/StdVector.h
	Eigen/src/StlSupport/details.h
	Eigen/src/SuperLUSupport/SuperLUSupport.h
	Eigen/src/UmfPackSupport/UmfPackSupport.h
	Eigen/src/misc/Image.h
	Eigen/src/misc/Kernel.h
	Eigen/src/misc/RealSvd2x2.h
	Eigen/src/misc/blas.h
	Eigen/src/misc/lapack.h
	Eigen/src/misc/lapacke.h
	Eigen/src/misc/lapacke_mangling.h
	Eigen/src/plugins/ArrayCwiseBinaryOps.h
	Eigen/src/plugins/ArrayCwiseUnaryOps.h
	Eigen/src/plugins/BlockMethods.h
	Eigen/src/plugins/CommonCwiseBinaryOps.h
	Eigen/src/plugins/CommonCwiseUnaryOps.h
	Eigen/src/plugins/MatrixCwiseBinaryOps.h
	Eigen/src/plugins/MatrixCwiseUnaryOps.h
)

find_library(ZLIB z HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT ZLIB)
    message(FATAL_ERROR "libz library not found")
endif()

find_library(CURLLIB curl HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT CURLLIB)
    message(FATAL_ERROR "libcurl library not found")
endif()

find_package (BZip2)
if (NOT BZIP2_FOUND)
    message(FATAL_ERROR "libbz2 library not found")
else()
    include_directories(${BZIP2_INCLUDE_DIRS})
    set(BZIP2 ${BZIP2_LIBRARIES})
endif()

find_library(LZMA lzma HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT LZMA)
    message(WARNING "liblzma library not found, if you specified --disable-lzma when compiling libhts, please ignore this warning")
    set(LZMA "")
else()
    set(LZMA ${LZMA_LIBRARIES})
endif()

find_library(DEFLATELIB libdeflate.a deflate HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ ../libdeflate/)
if(NOT DEFLATELIB)
    message(FATAL_ERROR "libdeflate library not found")
endif()

find_library(CRYPTOLIB crypto HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT CRYPTOLIB)
    message(FATAL_ERROR "libcrypto library not found")
endif()

find_package(Threads REQUIRED)
find_package(OpenSSL REQUIRED)

target_link_libraries(bluebear ${QGEN_LIBRARIES} ${HTS_LIBRARIES} ${ZLIB} ${BZIP2} ${LZMA} ${CURLLIB} ${DEFLATELIB} ${CRYPTOLIB} ${CMAKE_THREAD_LIBS_INIT})


install(TARGETS bluebear RUNTIME DESTINATION bin)
