/**
 * \file gpu_stl.hpp  STL-style datastructures for GPUs.
 */

#ifndef GRAPHLAB_GPU_STL_HPP
#define GRAPHLAB_GPU_STL_HPP

namespace graphlab {

  /**
   * CUDA-compatible struct for representing a pair<T1,T2>.
   */
  template <typename T1, typename T2>
  struct pair {
    T1 first;
    T2 second;
    pair(const T1& first, const T2& second)
      : first(first), second(second) { }
  }; // struct pair

  //! Like std::make_pair
  template <typename T1, typename T2>
  make_pair(const T1& first, const T2& second) {
    return pair<T1,T2>(first, second);
  }

  /**
   * CUDA-compatible class for vectors.
   *
   * @tparam T  Type of element stored in the vector.
   */
  template <typename T>
  class gpu_vector {

    // PUBLIC TYPES
    //==========================================================================
  public:

    typedef const T* const_iterator;
    typedef T* iterator;

    // PRIVATE DATA
    //==========================================================================
  private:

    //! Number of elements.
    size_t n;

    //! Pointer to first element.
    T* data;

    // PUBLIC METHODS
    //==========================================================================
  public:

    //! Default constructor.
    gpu_vector()
      : n(0), data(NULL) { }

    //! Constructor.
    gpu_vector(size_t n, T* data)
      : n(n), data(data) { }

    /**
     * Frees the device data.
     * Note: If this is, e.g., a gpu_vector<gpu_vector<SomeType> >,
     *       then you must manually free the vectors which this vector
     *       contains before you free this vector.
     */
    void free() {
      assert(data);
      gpu_free(data);
      n = 0;
      data = NULL;
    }

    //! Number of elements.
    size_t size() const { return n; }

    //! Returns a const reference to element i.
    //! WARNING: This is NOT bound-checked.
    const T& operator[](size_t i) const { return data[i]; }

    //! Returns a mutable reference element i.
    //! WARNING: This is NOT bound-checked.
    T& operator[](size_t i) { return data[i]; }

    /** Get a const iterator to the beginning. */
    const_iterator begin() const {
      return data;
    }

    /** Get a const iterator to the end. */
    const_iterator end() const {
      return T + n;
    }

    /** Get a mutable iterator to the beginning. */
    iterator begin() {
      return data;
    }

    /** Get a mutable iterator to the end. */
    iterator end() {
      return T + n;
    }

  }; // class gpu_vector

  namespace impl {

    //! Functor which returns the first element in a pair.
    template <typename T1, typename T2>
    struct pair_first_functor {
      const T1& operator()(const pair<T1,T2>& p) { return p.first; }
    };

    //! Functor which returns the second element in a pair.
    template <typename T1, typename T2>
    struct pair_second_functor {
      const T2& operator()(const pair<T1,T2>& p) { return p.second; }
    };

  } // end of namespace graphlab::impl

} // end of namespace graphlab

#endif // GRAPHLAB_GPU_STL_HPP
