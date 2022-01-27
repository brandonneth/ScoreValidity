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
    statement::For<I0, omp_parallel_for_exec,
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

  idx_t N = 256;
  idx_t R = 5;
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
      s = ",";
    }
    std::cout << "]";
    std::cout << " , ";

    auto time = 0.0;
    for(int i = 0; i < R; i++) {
      reset();
      auto start = std::chrono::high_resolution_clock::now();
      knl();
      auto stop = std::chrono::high_resolution_clock::now();
      time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    } 

    time = time / R;
 
    std::cout << time << "),\n";
  } while(std::next_permutation(perm.begin(), perm.end()));
}



int main(int RAJA_UNUSED_ARG(argc), char** RAJA_UNUSED_ARG(argv[]))
{

  using namespace RAJA;

 
  using namespace RAJA;
  using VIEW = View<double, Layout<2>>;

  idx_t N = 1028;
  idx_t R = 10;

  std::array<idx_t,2> s{{N,N}};
  auto l_01 = make_permuted_layout(s, {0,1});
  auto l_10 = make_permuted_layout(s, {1,0});

  VIEW a(new double[N*N], l_01);
  VIEW b(new double[N*N], l_01);
  VIEW t(new double[N*N], l_01);
  VIEW c(new double[N*N], l_01);
  VIEW d(new double[N*N], l_01);


  auto lam1 = [&](auto i, auto j, auto k) {
    t(i,j) += a(i,k) * b(k,j);
  };

  auto lam2 = [&](auto i, auto j, auto k) {
    d(i,j) += t(i,k) * c(k,j);
  };


  auto segs = make_tuple(RangeSegment(0,N),RangeSegment(0,N),RangeSegment(0,N));

  auto ra = permute_view(a,l_01);
  auto rb = permute_view(b,l_01);
  auto rc = permute_view(c,l_01);
  auto rd = permute_view(d,l_01);
  auto rt = permute_view(t,l_01);

  auto reset = [&](){ra();rb();rc();rd();rt();};

  auto knl1 = make_kernel<order_to_kpol3<0,1,2>::Policy>(segs, lam1);
  auto knl2 = make_kernel<order_to_kpol3<0,1,2>::Policy>(segs, lam2);

  auto dec = format_decisions(tie(b,c), knl1, knl2);

  std::array<decltype(l_01), 2> layouts{{l_01,l_10}};


  /*for(auto a0 = 0; a0 < 2; a0++) {
  for(auto a1 = 0; a1 < 2; a1++) {
  for(auto a2 = 0; a2 < 2; a2++) {*/
  for(auto b0 = 0; b0 < 2; b0++) {
  for(auto b1 = 0; b1 < 2; b1++) {
  for(auto b2 = 0; b2 < 2; b2++) {
  for(auto c0 = 0; c0 < 2; c0++) {
  for(auto c1 = 0; c1 < 2; c1++) {
  for(auto c2 = 0; c2 < 2; c2++) {
  /*for(auto d0 = 0; d0 < 2; d0++) {
  for(auto d1 = 0; d1 < 2; d1++) {
  for(auto d2 = 0; d2 < 2; d2++) {
  for(auto t0 = 0; t0 < 2; t0++) {
  for(auto t1 = 0; t1 < 2; t1++) {
  for(auto t2 = 0; t2 < 2; t2++) {*/
/*dec.set_format_before(a,layouts[a0],knl1);
dec.set_format_before(a,layouts[a1],knl2);
dec.set_format_after(a,layouts[a2],knl2);*/
dec.set_format_before(b,layouts[b0],knl1);
dec.set_format_before(b,layouts[b1],knl2);
dec.set_format_after(b,layouts[b2],knl2);
dec.set_format_before(c,layouts[c0],knl1);
dec.set_format_before(c,layouts[c1],knl2);
dec.set_format_after(c,layouts[c2],knl2);
/*dec.set_format_before(d,layouts[d0],knl1);
dec.set_format_before(d,layouts[d1],knl2);
dec.set_format_after(d,layouts[d2],knl2);
dec.set_format_before(t,layouts[t0],knl1);
dec.set_format_before(t,layouts[t1],knl2);
dec.set_format_after(t,layouts[t2],knl2);*/

  auto knl = dec.finalize();

  dec.print_to_stream(std::cout);
 
  auto start =std::chrono::high_resolution_clock::now();
  for(auto r = 0; r < R; r++) {
    knl();
  }
  auto stop =std::chrono::high_resolution_clock::now();
 
  auto t =std::chrono::duration_cast<std::chrono::nanoseconds>(stop-start).count();
 
  std::cout << " " << t / 5 << "\n";
  //}}}
  //}}}
  //}}}
  }}}
  }}}

  //now get the same stuff for the model selected variant
  auto fresh_dec = format_decisions(tie(b,c),knl1, knl2);

  std::cout << "\nModel: ";
  
  auto modelselected = fresh_dec.finalize();
  fresh_dec.print_to_stream(std::cout);


}
