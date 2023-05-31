#ifndef TEST_HXX
#define TEST_HXX

/* test.hxx */
/* (c) Liwei Ji 2023-05-31 */
/* Produced with Julia */

#include <cctk.h>

namespace WaveToyHigherOrderX {

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE T
calc_fd_c1D2O(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  return ((1/T(12))*gf(p.I - 2*p.DI[dir]) + (2/T(3))*gf(p.I + p.DI[dir]) - (1/T(12))*gf(p.I + 2*p.DI[dir]) - (2/T(3))*gf(p.I - p.DI[dir])) / p.DX[dir];
}

} // namespace WaveToyHigherOrderX

#endif // #ifndef TEST_HXX
