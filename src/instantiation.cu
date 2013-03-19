#include "c2r.h"
#include "register_ops.h"
#include "sm.h"

#ifndef INSTANTIATED_TYPE
#define INSTANTIATED_TYPE double
#endif

#ifndef WPT
#define WPT 5
#endif

#ifndef SM
#define SM sm_20
#endif

namespace inplace {
namespace detail {

template __global__ void register_col_op<SM, INSTANTIATED_TYPE, prerotator, WPT>(int, int, INSTANTIATED_TYPE*, prerotator);
template __global__ void register_row_shuffle<SM, INSTANTIATED_TYPE, WPT>(int, int, INSTANTIATED_TYPE*, shuffle);
template __global__ void register_col_op<SM, INSTANTIATED_TYPE, postpermuter, WPT>(int, int, INSTANTIATED_TYPE*, postpermuter);

}
}
