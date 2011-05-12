/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GRAPHLAB_SERIALIZE_VECTOR_HPP
#define GRAPHLAB_SERIALIZE_VECTOR_HPP
#include <vector>
#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/serialization/iterator.hpp>


namespace graphlab {
  namespace archive_detail {
    /**
     * We re-dispatch vectors because based on the contained type,
     * it is actually possible to serialize them like a POD
     */
    template <typename ArcType, typename ValueType, bool IsPOD>
    struct vector_serialize_impl {
      static void exec(ArcType &a, const ValueType& vec) {
        // really this is an assert false. But the static assert
        // must depend on a template parameter 
        BOOST_STATIC_ASSERT(sizeof(ArcType) == 0);
        assert(false);
      };
    };
    /**
     * We re-dispatch vectors because based on the contained type,
     * it is actually possible to deserialize them like a POD
     */
    template <typename ArcType, typename ValueType, bool IsPOD>
    struct vector_deserialize_impl {
      static void exec(ArcType &a, ValueType& vec) {
        // really this is an assert false. But the static assert
        // must depend on a template parameter 
        BOOST_STATIC_ASSERT(sizeof(ArcType) == 0);
        assert(false);
      };
    };
    
    /// If contained type is not a POD use the standard serializer
    template <typename ArcType, typename ValueType>
    struct vector_serialize_impl<ArcType, ValueType, false > {
      static void exec(ArcType &a, const std::vector<ValueType>& vec) {
        serialize_impl<ArcType, size_t, false>::exec(a, vec.size());
        serialize_iterator(a,vec.begin(), vec.end());
      }
    };

    /// Fast vector serialization if contained type is a POD
    template <typename ArcType, typename ValueType>
    struct vector_serialize_impl<ArcType, ValueType, true > {
      static void exec(ArcType &a, const std::vector<ValueType>& vec) {
        serialize_impl<ArcType, size_t, false>::exec(a, vec.size());
        serialize(a, &(vec[0]),sizeof(ValueType)*vec.size());
      }
    };

    /// If contained type is not a POD use the standard deserializer
    template <typename ArcType, typename ValueType>
    struct vector_deserialize_impl<ArcType, ValueType, false > {
      static void exec(ArcType& a, std::vector<ValueType>& vec){
        size_t len;
        deserialize_impl<ArcType, size_t, false>::exec(a, len);
        vec.clear(); vec.reserve(len);
        deserialize_iterator<ArcType, ValueType>(a, std::inserter(vec, vec.end()));
      }
    };

    /// Fast vector deserialization if contained type is a POD
    template <typename ArcType, typename ValueType>
    struct vector_deserialize_impl<ArcType, ValueType, true > {
      static void exec(ArcType& a, std::vector<ValueType>& vec){
        size_t len;
        deserialize_impl<ArcType, size_t, false>::exec(a, len);
        vec.clear(); vec.resize(len);
        deserialize(a, &(vec[0]), sizeof(ValueType)*vec.size());
      }
    };

    
    
    /**
       Serializes a vector */
    template <typename ArcType, typename ValueType>
    struct serialize_impl<ArcType, std::vector<ValueType>, false > {
      static void exec(ArcType &a, const std::vector<ValueType>& vec) {
        vector_serialize_impl<ArcType, ValueType, gl_is_pod<ValueType>::value>::exec(a, vec);
      }
    };
    /**
       deserializes a vector */
    template <typename ArcType, typename ValueType>
    struct deserialize_impl<ArcType, std::vector<ValueType>, false > {
      static void exec(ArcType& a, std::vector<ValueType>& vec){
        vector_deserialize_impl<ArcType, ValueType, gl_is_pod<ValueType>::value>::exec(a, vec);
      }
    };
  } // archive_detail
} // namespace graphlab

#endif 
