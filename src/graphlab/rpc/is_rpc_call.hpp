#ifndef IS_RPC_CALL_HPP
#define IS_RPC_CALL_HPP
#include <boost/type_traits/remove_pointer.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/type_traits/function_traits.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/less.hpp>
#include <boost/mpl/comparison.hpp>
#include <boost/mpl/int.hpp>

#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/function_arg_types_def.hpp>

/**
An RPC-aware call is a function of the form:
... fn(distributed_control& dc, procid_t source....)
The code here checks for that.
*/

namespace graphlab {
class distributed_control;
namespace dc_impl {

namespace is_rpc_call_detail {
/**
Whether the function has less than or equal to 2 arguments
*/
template <typename F>
struct less_than_2_args {
  typedef typename boost::mpl::bool_<FARITY < 2 >::type type;  
};


/**
Now, arg1_type and arg_2 type may not exist in function_traits if the 
number of arguments is < 2. I will need to wrap it to make it safe
*/
template <typename F, size_t nargs>
struct get_args{
 typedef NIF0 arg1_type;
 typedef NIF1 arg2_type;
};

// if 0 args. then make both void
template <typename F>
struct get_args<F, 0>{
 typedef void arg1_type;
 typedef void arg2_type;
};

// if 1 arg then make just make arg2 void
template <typename F>
struct get_args<F, 1>{
 typedef NIF0 arg1_type;
 typedef void arg2_type;
};



template <typename F>
struct check_first_arg {
  typedef typename boost::is_same<typename get_args<F,FARITY>::arg1_type, distributed_control>::type type;  
};

template <typename F>
struct check_second_arg {
  typedef typename boost::is_integral<typename get_args<F,FARITY>::arg2_type>::type type;  
};


}
template <typename F>
struct is_rpc_call {
  typedef typename boost::mpl::if_< typename is_rpc_call_detail::less_than_2_args<F>::type,
               boost::false_type,
               typename boost::mpl::and_<
                    typename is_rpc_call_detail::check_first_arg<F>::type, 
                    typename is_rpc_call_detail::check_second_arg<F>::type>::type >::type type;
               
               
};

// Varargs are all none RPC calls
#define BLOCK_VAR_ARGS(Z,N,_)  \
template <typename RetType BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
struct is_rpc_call<RetType (BOOST_PP_ENUM_PARAMS(N, T) BOOST_PP_COMMA_IF(N) ...)> { \
   typedef boost::false_type type; \
}; \
\
template <typename RetType BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
struct is_rpc_call<RetType (*)(BOOST_PP_ENUM_PARAMS(N, T) BOOST_PP_COMMA_IF(N) ...)> { \
   typedef boost::false_type type; \
};
BOOST_PP_REPEAT(6, BLOCK_VAR_ARGS, _)
#undef BLOCK_VAR_ARGS

#define GEN_GET_USER_ARG(Z,N,_)  \
template <typename F, typename BoolType>  \
struct BOOST_PP_CAT(get_cleaned_rpc_or_basic_arg, N) { \
  typedef BOOST_PP_CAT(NIF, N) arg_type;  \
};  \
template <typename F> \
struct BOOST_PP_CAT(get_cleaned_rpc_or_basic_arg, N) <F,  boost::mpl::bool_<true> > {  \
  typedef BOOST_PP_CAT(F, N) arg_type;  \
};  \
template <typename F>   \
struct BOOST_PP_CAT(get_cleaned_user_arg, N) {  \
  typedef typename BOOST_PP_CAT(get_cleaned_rpc_or_basic_arg, N)<F,typename is_rpc_call<F>::type>::arg_type arg_type; \
};

BOOST_PP_REPEAT(6, GEN_GET_USER_ARG, _)
#undef GEN_GET_USER_ARG

} // namespace dc_impl
} // namespace graphlab

#include <graphlab/rpc/function_arg_types_undef.hpp>

#endif