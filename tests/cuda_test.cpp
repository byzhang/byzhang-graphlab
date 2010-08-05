
#include <cassert>
#include <iostream>

#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include <cublas.h>

static void check(cublasStatus status, const std::string& msg) {
  if (status != CUBLAS_STATUS_SUCCESS) {
    std::cerr << msg << "\n"
              << "CUBLAS status: " << status << std::endl;
    assert(false);
  }
}

static void testCUBLAS() {
  size_t n(1000);
  assert(n > 0);

  // Init CUBLAS, and set device vector to {0,...,n-1}.
  check(cublasInit(), "CUBLAS initialization error!");
  float* d_vec = NULL;
  check(cublasAlloc(n, sizeof(float), (void**)&d_vec),
        "CUBLAS device memory allocation error!\n");
  float* h_vec = new float[n];
  for (size_t i(0); i < n; ++i)
    h_vec[i] = i;
  check(cublasSetVector(n, sizeof(float), h_vec, 1, d_vec, 1),
        "CUBLAS failed when copying host vector to device.\n");
  delete [] h_vec;
  h_vec = NULL;

  // Compute sum(device vector).
  float total(cublasSasum(n, d_vec, 1));
  check(cublasGetError(), "CUBLAS error while summing device vector!\n");
  std::cout << "CUBLAS computed sum(0,...," << (n-1) << ") = " << total
            << std::endl;

  check(cublasFree(d_vec),
        "CUBLAS error while freeing vec on device!\n");
  d_vec = NULL;
}

static void testThrust() {
  size_t n(1000);
  assert(n > 0);
  thrust::host_vector<float> h_vec(n);
  float true_sum(0);
  for (size_t i(0); i < n; ++i) {
    h_vec[i] = i;
    true_sum += i;
  }
  thrust::device_vector<float> d_vec(h_vec);
  float newval = 100;
  true_sum -= h_vec[0];
  true_sum += newval;
  d_vec[0] = newval;
  float sum =
    thrust::reduce(d_vec.begin(), d_vec.end(), (float)0, thrust::plus<float>());
  if (sum != true_sum) {
    std::cerr << "Thrust computed incorrect sum " << sum
              << " instead of the correct sum " << true_sum << std::endl;
    assert(false);
  } else {
    std::cout << "Thrust computed the correct sum " << sum << std::endl;
  }
}

int main(int argc, char** argv) {

  testCUBLAS();
  testThrust();

  return 0;

}
