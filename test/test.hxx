#ifndef TEST_HXX
#define TEST_HXX

/* test.hxx */
/* (c) Liwei Ji 2023-06-05 */
/* Produced with Julia */

#include <cctk.h>

namespace WaveToyHigherOrderX {

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE T
calc_fd_c1D2O(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  return ((1/T(2))*gf(p.I + p.DI[dir]) - (1/T(2))*gf(p.I - p.DI[dir])) / p.DX[dir];
}

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE T
calc_fd_l1D1O(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  return (gf(p.I) - gf(p.I - p.DI[dir])) / p.DX[dir];
}

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE T
calc_fd_r1D1O(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  return (gf(p.I + p.DI[dir]) - gf(p.I)) / p.DX[dir];
}

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE T
calc_fd_c2D2O(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  return (gf(p.I - p.DI[dir]) + gf(p.I + p.DI[dir]) - (2/T(1))*gf(p.I)) / (pow(p.DX[dir], 2));
}

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE T
calc_fd_l2D1O(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  return (gf(p.I - 2*p.DI[dir]) + gf(p.I) - (2/T(1))*gf(p.I - p.DI[dir])) / (pow(p.DX[dir], 2));
}

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE T
calc_fd_r2D1O(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  return (gf(p.I) + gf(p.I + 2*p.DI[dir]) - (2/T(1))*gf(p.I + p.DI[dir])) / (pow(p.DX[dir], 2));
}

} // namespace WaveToyHigherOrderX

#endif // #ifndef TEST_HXX
