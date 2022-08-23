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

std::chrono::time_point<std::chrono::high_resolution_clock> _start;
void start() {
   _start = std::chrono::high_resolution_clock::now();
}

auto stop() {
  auto end = std::chrono::high_resolution_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - _start).count();
  return duration;
}

auto run(std::string variant) {
  start();
  using namespace RAJA;

 
  using namespace RAJA;
  using VIEW = View<double, Layout<2>>;

  idx_t ni = 1024;
  idx_t nj = 1024;
  idx_t nk = 1024;
  idx_t nl = 1024;
  idx_t nm = 1024;
  VIEW A(new double[ni*nk], ni,nk);
  VIEW B(new double[nj*nk], nk,nj);
  VIEW C(new double[nj*nm], nj,nm);
  VIEW D(new double[nm*nl], nm,nl);
  VIEW E(new double[ni*nj], ni,nj);
  VIEW F(new double[nj*nl], nj,nl);
  VIEW G(new double[ni*nl], ni,nl);
  
  auto lam_init_e = [&](auto i, auto j) {
    E(i,j) = 0;
  };
  auto lam_init_f = [&](auto i, auto j) {
    F(i,j) = 0;
  };
  auto lam_init_g = [&](auto i, auto j) {
    G(i,j) = 0;
  };

  auto lam_comp1 = [&](auto i, auto j, auto k) {
    E(i,j) += A(i,k) * B(k,j);
  };

  auto lam_comp2 = [&](auto i, auto j, auto k) {
    F(i,j) += C(i,k) * D(k,j);
  };
  auto lam_comp3 = [&](auto i, auto j, auto k) {
    G(i,j) += E(i,k) * F(k,j);
  };

  using namespace RAJA;
  using INIT_POL = KernelPolicy<
    statement::For<0, omp_parallel_for_exec,
      statement::For<1, loop_exec,
        statement::Lambda<0>
      >
    >
      >;
    
  using COMP_POL = KernelPolicy<
    statement::For<0, omp_parallel_for_exec,
      statement::For<1, loop_exec,
        statement::For<2, loop_exec,
          statement::Lambda<0>
        >
      >
    >
  >;
  auto eseg = make_tuple(RangeSegment(0,ni), RangeSegment(0,nj));
  auto fseg = make_tuple(RangeSegment(0,nj), RangeSegment(0,nl));
  auto gseg = make_tuple(RangeSegment(0,ni), RangeSegment(0,nl));
  
  auto seg1 = make_tuple(RangeSegment(0,ni), RangeSegment(0,nj), RangeSegment(0,nk));
  auto seg2 = make_tuple(RangeSegment(0,nj), RangeSegment(0,nl), RangeSegment(0,nm));
  auto seg3 = make_tuple(RangeSegment(0,ni), RangeSegment(0,nl), RangeSegment(0,nj));


  auto init_e = make_kernel<INIT_POL>(eseg, lam_init_e);
  auto comp1 = make_kernel<COMP_POL>(seg1, lam_comp1);
  auto init_f = make_kernel<INIT_POL>(fseg, lam_init_f);
  auto comp2 = make_kernel<COMP_POL>(seg2, lam_comp2);
  auto init_g = make_kernel<INIT_POL>(gseg, lam_init_g);
  auto comp3 = make_kernel<COMP_POL>(seg3, lam_comp3);

  std::array<idx_t, 2> bsizes{{nk,nj}};
  auto blayout = make_permuted_layout(bsizes, {{1,0}});
  std::array<idx_t, 2> dsizes{{nm,nl}};
  auto dlayout = make_permuted_layout(dsizes, {{1,0}});
  std::array<idx_t, 2> fsizes{{nj,nl}};
  auto flayout01 = make_permuted_layout(fsizes, {{0,1}});
  auto flayout10 = make_permuted_layout(fsizes, {{1,0}});
  
  

  auto dec = format_decisions(tie(B, D, F), init_e, comp1, init_f, comp2, init_g, comp3);
  
  dec.set_format_for(B, blayout, comp1);
  dec.set_format_for(D, blayout, comp2);
  dec.set_format_for(F, flayout01, init_f);
  dec.set_format_for(F, flayout10, comp3);

  if(variant == "lock") {
    dec.lock();
  }
  auto initializing_time = stop();
  start();
  auto comp = dec.finalize();
  auto solve_time = stop();
  start();
  comp();
  auto kernel_time = stop();

  auto breakdown = dec.time_execution();
  auto conversion_time = camp::get<0>(breakdown);
  auto computation_time = camp::get<1>(breakdown);


  std::cout << variant << ",Initialization," << initializing_time << "\n";
  std::cout << variant << ",Constraint Setup," << dec.setup_time << "\n";
  std::cout << variant << ",ISL Setup," << dec.isl_time << "\n";
  std::cout << variant << ",ISL Solve," << dec.solve_time << "\n";
  std::cout << variant << ",Kernels," << kernel_time << "\n";
  std::cout << variant << ",Conversion," << conversion_time << "\n";
  std::cout << variant << ",Computation," << computation_time << "\n";
  
}


int main(int RAJA_UNUSED_ARG(argc), char** RAJA_UNUSED_ARG(argv[]))
{
  std::cout << "Variant,Component,Time\n";
  run("standard");
  run("standard");
  run("lock");
  run("lock");
  

}
