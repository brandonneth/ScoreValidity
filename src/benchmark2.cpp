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
idx_t N = 512;
idx_t R = 5;
template <typename Policy, typename L1, typename L2, typename ViewType>
void enumerate_layouts(L1 & reset_lam, L2 & comp_lam, ViewType & a, ViewType & b) {

  using namespace RAJA;

  auto segs = make_tuple(RangeSegment(0,N), RangeSegment(0,N), RangeSegment(0,N));
  auto reset_knl = make_kernel<Policy>(segs,reset_lam);


  std::array<idx_t, 2> perm{{0,1}};

  do {
    std::cerr << "    layout order " << perm[0] << perm[1] << "\n";
    reset_knl();
    auto layout = make_permuted_layout({{N,N}}, perm);
    a.set_layout(layout);
    b.set_layout(layout);

    auto comp_knl = make_kernel<Policy>(segs, comp_lam);
    
    auto access = comp_knl.execute_symbolically().at(0);
    auto normalized = comp_knl.normalize_access(access);

    std::cout << "(";
    std::cout << "[";
    std::string s = "";
    for(auto i : normalized) {
      std::cout << s << i;
      s = ",";
    }
    std::cout << "]";
    std::cout << " , ";

    auto time = 0.0;
    for(int i = 0; i < R; i++) {
      reset_knl();
      auto start = std::chrono::high_resolution_clock::now();
      comp_knl();
      auto stop = std::chrono::high_resolution_clock::now();
      time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    } 

    time = time / R;
 
    std::cout << time << "," << N << "),\n";
  } while(std::next_permutation(perm.begin(), perm.end()));
}

template <typename L1, typename L2, typename ViewType>
void enumerate_policies(L1 & reset_lam, L2 & comp_lam, ViewType & a, ViewType & b) {

  std::cerr << "  policy order 012\n";
  enumerate_layouts<order_to_kpol3<0,1,2>::Policy>(reset_lam, comp_lam, a, b);
  std::cerr << "  policy order 021\n";
  enumerate_layouts<order_to_kpol3<0,2,1>::Policy>(reset_lam, comp_lam, a, b);
  std::cerr << "  policy order 102\n";
  enumerate_layouts<order_to_kpol3<1,0,2>::Policy>(reset_lam, comp_lam, a, b);
  std::cerr << "  policy order 120\n";
  enumerate_layouts<order_to_kpol3<1,2,0>::Policy>(reset_lam, comp_lam, a, b);
  std::cerr << "  policy order 201\n";
  enumerate_layouts<order_to_kpol3<2,0,1>::Policy>(reset_lam, comp_lam, a, b);
  std::cerr << "  policy order 210\n";
  enumerate_layouts<order_to_kpol3<2,1,0>::Policy>(reset_lam, comp_lam, a, b);

}

int main(int RAJA_UNUSED_ARG(argc), char** RAJA_UNUSED_ARG(argv[]))
{

  using namespace RAJA;
  using VIEW = View<double, Layout<2>>;

  std::cout << "[";
  N = 128;
  while (N <= 1024) {
  VIEW a(new double[N*N], N,N);
  VIEW b(new double[N*N], N,N);



  auto reset_lam =  [&](auto i0, auto i1, auto ) {b(i0, i1) = std::rand();};

  std::cerr << "arg order 01\n";
  auto lambda01 = [&](auto i0, auto i1, auto ) {a(i0,i1) = b(i0,i1);};
  enumerate_policies(reset_lam, lambda01, a, b);
  std::cerr << "arg order 10\n";
  auto lambda10 = [&](auto i0, auto i1, auto ) {a(i1,i0) = b(i1,i0);};
  enumerate_policies(reset_lam, lambda10, a, b);
  std::cerr << "arg order 02\n";
  auto lambda02 = [&](auto i0, auto, auto i2) {a(i0,i2) = b(i0,i2);};
  enumerate_policies(reset_lam, lambda02, a, b);
  std::cerr << "arg order 20\n";
  auto lambda20 = [&](auto i0, auto, auto i2) {a(i2,i0) = b(i2,i0);};
  enumerate_policies(reset_lam, lambda20, a, b);
  std::cerr << "arg order 21\n";
  auto lambda21 = [&](auto, auto i1, auto i2) {a(i2,i1) = b(i2,i1);};
  enumerate_policies(reset_lam, lambda21, a, b);
  std::cerr << "arg order 12\n";
  auto lambda12 = [&](auto, auto i1, auto i2) {a(i1,i2) = b(i1,i2);};
  enumerate_policies(reset_lam, lambda12, a, b);
  N = N * 2;
  }

  std::cout << "]";
  return 0;
}
