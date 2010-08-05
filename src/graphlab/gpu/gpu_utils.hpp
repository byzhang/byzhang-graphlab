/**
 * \file gpu_utils.hpp  GPU utilities.
 */

#ifndef GRAPHLAB_GPU_UTILS_HPP
#define GRAPHLAB_GPU_UTILS_HPP

#include <sstream>

namespace graphlab {

  //! Print to string.
  template <typename T>
  std::string to_string(const T& t) {
    std::ostringstream out;
    out << t;
    return out.str();
  }

  /**
   * Check for a CUDA error.  If one is found, throw a std::runtime_error
   * with the given error message.
   * @param e    Error code
   * @param msg  Error message
   */
  inline void cudaCheck(cudaError_t e, const std::string& msg) {
    if (e != cudaSuccess)
      throw std::runtime_error("CUDA Error (" + to_string(e) + "): " + msg);
  }

  /**
   * Allocate GPU global memory space for n objects of type T.
   * This throws a std::runtime_error if the allocation fails.
   *
   * @param n  Number of objects of type T to allocate.
   *
   * @return  Device pointer to allocated memory.
   */
  template <typename T>
  T* gpu_alloc(size_t n) {
    T* d_ptr = NULL;
    cudaCheck(cudaMalloc((void**) &d_ptr, sizeof(T)),
              "gpu_alloc failed");
    return d_ptr;
  }

  /**
   * Free GPU global memory space.
   * This throws a std::runtime_error if the free fails.
   *
   * @param  d_ptr  Device pointer to allocated memory.
   */
  template <typename T>
  void gpu_free(T* d_ptr) {
    cudaCheck(cudaFree(d_ptr),
              "gpu_free failed");
  }

  /**
   * Copy CPU memory to the GPU (n items of type T from h_ptr to d_ptr).
   * This throws a std::runtime_error if the free fails.
   *
   * @param  d_ptr  Device pointer to allocated memory.
   */
  template <typename T>
  void copy_host_to_device(T* h_ptr, size_t n, T* d_ptr) {
    cudaCheck(cudaMemcpy(d_ptr, h_ptr, n, cudaMemcpyHosttoDevice),
              "copy_host_to_device failed");
  }

} // end of namespace graphlab

#endif // GRAPHLAB_GPU_UTILS_HPP
