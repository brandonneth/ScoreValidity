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

  idx_t N = 2048;
  idx_t R = 10;

  std::array<idx_t,2> s{{N,N}};
  auto l_01 = make_permuted_layout(s, {0,1});
  auto l_10 = make_permuted_layout(s, {1,0});

  VIEW a(new double[N*N], l_10);
  VIEW b(new double[N*N], l_01);
  VIEW c(new double[N*N], l_10);


  auto matmul_lam = [&](auto i, auto j, auto k) {
    c(i,j) += a(i,k) * b(k,j);
  };

  auto segs = make_tuple(RangeSegment(0,N),RangeSegment(0,N),RangeSegment(0,N));

  

  auto matmul = make_kernel<order_to_kpol3<0,1,2>::Policy>(segs, matmul_lam);
  
  auto init_segs = make_tuple(RangeSegment(0,N), RangeSegment(0,N));
  auto init = make_kernel<KernelPolicy<statement::For<0,omp_parallel_for_exec,statement::For<1,loop_exec,statement::Lambda<0>>>>>(init_segs, [&](auto i, auto j) {
    a(i,j) = std::rand();
    b(i,j) = std::rand();
    c(i,j) = std::rand();
  });

  init();
  
  auto abad = permute_view(a, l_10); 
  auto bbad = permute_view(b, l_01); 
  auto cbad = permute_view(c, l_10); 

  auto agood = permute_view(a, l_01); 
  auto bgood = permute_view(b, l_10); 
  auto cgood = permute_view(c, l_01); 


  //set up the bad layouts
  abad();
  bbad();
  cbad();

  //run and time the kernel.
  auto badtime = 0.0;
  for(int i = 0; i < R; i++) {
    auto start = std::chrono::high_resolution_clock::now();
    matmul();
    auto stop = std::chrono::high_resolution_clock::now();
    badtime += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
  } 

  std::cout  << badtime << "\n";

  //run and time the conversion to good layouts, reseting them to the bad one each time.
  auto changetime = 0.0;
  for(int i = 0; i < R; i++) {
    auto start = std::chrono::high_resolution_clock::now();
    agood();
    bgood();
    cgood();
    auto stop = std::chrono::high_resolution_clock::now();
    changetime += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    abad();
    bbad();
    cbad();
  } 

  std::cout  << changetime << "\n";

  
  //run and time the comp for good layouts
  agood();
  bgood();
  cgood();
  auto goodtime = 0.0;
  for(int i = 0; i < R; i++) {
    auto start = std::chrono::high_resolution_clock::now();
    matmul();
    auto stop = std::chrono::high_resolution_clock::now();
    goodtime += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
  } 
  std::cout << goodtime << "\n";

  return 0;
}
