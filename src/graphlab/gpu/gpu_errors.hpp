/**
 * \file gpu_error.hpp  Error codes for GPU GraphLab.
 */

#ifndef GRAPHLAB_GPU_ERRORS_HPP
#define GRAPHLAB_GPU_ERRORS_HPP

namespace graphlab {

  namespace gpu_errors {

    /**
     * Error codes for GPU GraphLab.
     *  - NO_ERROR
     *  - OUT_OF_BOUNDS
     *  - ELEMENT_NOT_FOUND
     */
    enum gpu_error_type { NO_ERROR, OUT_OF_BOUNDS, ELEMENT_NOT_FOUND };

  } // end of namespace gpu_errors

} // end of namespace graphlab

#endif // GRAPHLAB_GPU_ERRORS_HPP
