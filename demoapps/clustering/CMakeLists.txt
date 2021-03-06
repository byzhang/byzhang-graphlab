project(GraphLab)
if (HAS_EIGEN)
  include_directories(../../deps/eigen-eigen-3.0.2/)
  add_definitions(-DHAS_EIGEN)
endif()

set(OMP_SUPPORT 0)
  if (OMP_SUPPORT)
    add_definitions(-omp)
    add_definitions(-DOMP_SUPPORT)
  endif()

if (HAS_EIGEN)
  set( lapack_found 1) # clustering does not support eigen yet (becuase of experimental svd)
elseif(ITPP-FOUND)
  set( lapack_found 1)
endif()

if (lapack_found)
  add_graphlab_executable(glcluster clustering.cpp advanced_config.cpp 
read_matrix_market.cpp ../../libs/matrixmarket/mmio.cpp distance.cpp   
sparse_math.cpp main.cpp kmeans_plus_plus.cpp dumpcluster.cpp tfidf.cpp
lda.cpp gamma.cpp kcores.cpp itemknn.cpp)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/clive_test
    ${CMAKE_CURRENT_BINARY_DIR}/clive_test COPYONLY)
 endif()

