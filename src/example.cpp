//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2016-19, Lawrence Livermore National Security, LLC
// and RAJA project contributors. See the RAJA/COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//////////////////////////////////////////////////////////////////////////////
#include "RAJA/RAJA.hpp"
#include <algorithm>
#include <chrono>
using namespace RAJA;
template <idx_t I0, idx_t I1, idx_t I2>
struct order_to_kpol3 {
  using Policy = KernelPolicy<
    statement::For<I0, loop_exec,
      statement::For<I1, loop_exec,
        statement::For<I2, loop_exec,
          statement::Lambda<0>
        >
      >
    >
  >;
};

template <typename Policy>
void enumerate_layouts() {
  using namespace RAJA;
  using VIEW = View<double, Layout<3>>;

  idx_t N = 128;

  VIEW a(new double[N*N*N], N,N,N);
  VIEW b(new double[N*N*N], N,N,N);

  auto lambda = [&](auto i0, auto i1, auto i2) {a(i0,i1,i2) = b(i0,i1,i2);};

  auto segs = make_tuple(RangeSegment(0,N), RangeSegment(0,N), RangeSegment(0,N));

  auto reset = make_kernel<Policy>(segs, [=](auto i0, auto i1, auto i2) {b(i0, i1, i2) = std::rand();});


  std::array<idx_t, 3> perm{{0,1,2}};

  do {
    reset();
    auto layout = make_permuted_layout({{N,N,N}}, perm);
    a.set_layout(layout);
    b.set_layout(layout);

    auto knl = make_kernel<Policy>(segs, lambda);
    
    auto access = knl.execute_symbolically().at(0);
    auto normalized = knl.normalize_access(access);

    std::cout << "(";
    std::cout << "[";
    std::string s = "";
    for(auto i : normalized) {
      std::cout << s << i;
      s = " ";
    }
    std::cout << "]";
    std::cout << " , ";
 
    auto start = std::chrono::high_resolution_clock::now();

    knl();
    auto stop = std::chrono::high_resolution_clock::now();

    auto time = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();

    std::cout << time << "),\n";
  } while(std::next_permutation(perm.begin(), perm.end()));
}



int main(int RAJA_UNUSED_ARG(argc), char** RAJA_UNUSED_ARG(argv[]))
{

  using namespace RAJA;

  std::cout << "[";
  enumerate_layouts<order_to_kpol3<0,1,2>::Policy>();
  enumerate_layouts<order_to_kpol3<0,2,1>::Policy>();
  enumerate_layouts<order_to_kpol3<1,0,2>::Policy>();
  enumerate_layouts<order_to_kpol3<1,2,0>::Policy>();
  enumerate_layouts<order_to_kpol3<2,0,1>::Policy>();
  enumerate_layouts<order_to_kpol3<2,1,0>::Policy>();

  std::cout << "]";

  return 0;
}
