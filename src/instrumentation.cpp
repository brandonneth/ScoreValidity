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

auto write_datapoint = [](auto experiment, auto dimensionality, auto problemSize, auto views, auto constraints, auto variant, auto component, auto time) {
  std::cout << experiment << "," << dimensionality << "," << problemSize << "," << views << "," << constraints << "," << variant << "," << component << "," << time << "\n";
};

std::chrono::time_point<std::chrono::high_resolution_clock> start_;
void start() {
   start_ = std::chrono::high_resolution_clock::now();
}

auto stop() {
  auto end = std::chrono::high_resolution_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start_).count();
  return duration;
}

// auto run_3d_3knl_1var_nDecisions(std::string variant, camp::idx_t numDecisions, std::size_t ps) {
  
//   std::string dimensionality = "3D1Var";

//   std::stringstream vv;
//   vv << variant << "_" << numDecisions;
//   variant = vv.str();
//   start();

//   using namespace RAJA;
//   using VIEW = View<double, Layout<3>>;
//   double root = std::pow(ps, 1.0/3.0);
//   camp::idx_t  n = (camp::idx_t) root;
//   VIEW A(new double[n*n*n], n,n,n);
//   VIEW B(new double[n*n*n], n,n,n);
//   VIEW C(new double[n*n*n], n,n,n);
//   VIEW D(new double[n*n*n], n,n,n);
//   VIEW E(new double[n*n*n], n,n,n);
//   VIEW F(new double[n*n*n], n,n,n);
//   VIEW G(new double[n*n*n], n,n,n);

//   auto lam_init_e = [&](auto i, auto j, auto k) {
//     E(i,j,k) = 0;
//   };
//   auto lam_init_f = [&](auto i, auto j, auto k) {
//     F(i,j,k) = 0;
//   };
//   auto lam_init_g = [&](auto i, auto j, auto k) {
//     G(i,j,k) = 0;
//   };

//   auto lam_comp1 = [&](auto l, auto i, auto j, auto k) {
//     E(l,i,j) += A(l,i,k) * B(l,k,j);
//   };
//   auto lam_comp2 = [&](auto l, auto i, auto j, auto k) {
//     F(l,i,j) += C(l,i,k) * D(l,k,j);
//   };
//   auto lam_comp3 = [&](auto l, auto i, auto j, auto k) {
//     G(l,i,j) += E(l,i,k) * F(l,k,j);
//   };

//   using namespace RAJA;
//   using INIT_POL = KernelPolicy<
//     statement::For<0, omp_parallel_for_exec,
//       statement::For<1, loop_exec,
//         statement::For<2, loop_exec,
//           statement::Lambda<0>
//         >
//       >
//     >
//   >;
    
//   using COMP_POL = KernelPolicy<
//     statement::For<0, omp_parallel_for_exec,
//       statement::For<1, loop_exec,
//         statement::For<2, loop_exec,
//           statement::For<3, loop_exec,
//             statement::Lambda<0>
//           >
//         >
//       >
//     >
//   >;
//   auto eseg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));
//   auto fseg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));
//   auto gseg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));

//   auto comp_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));

//   auto init_e = make_kernel<INIT_POL>(eseg, lam_init_e);
//   auto comp1 = make_kernel<COMP_POL>(comp_seg, lam_comp1);
//   auto init_f = make_kernel<INIT_POL>(fseg, lam_init_f);
//   auto comp2 = make_kernel<COMP_POL>(comp_seg, lam_comp2);
//   auto init_g = make_kernel<INIT_POL>(gseg, lam_init_g);
//   auto comp3 = make_kernel<COMP_POL>(comp_seg, lam_comp3);

//   if (variant.find("original") != std::string::npos) {
//     init_e();
//     comp1();
//     init_f();
//     comp2();
//     init_g();
//     comp3();
//     auto run_time = stop();
//     write_datapoint(dimensionality, ps, variant, "Computation", run_time);
//     return;
//   }


//   std::array<idx_t, 3> bsizes{{n,n,n}};
//   auto blayout = make_permuted_layout(bsizes, {{2,1,0}});
//   std::array<idx_t, 3> dsizes{{n,n,n}};
//   auto dlayout = make_permuted_layout(dsizes, {{2,1,0}});
//   std::array<idx_t, 3> fsizes{{n,n,n}};
//   auto flayout01 = make_permuted_layout(fsizes, {{2,0,1}});
//   auto flayout10 = make_permuted_layout(fsizes, {{2,1,0}});

//   if (variant.find("hand") != std::string::npos) {
//     camp::idx_t run_time = 0;
//     camp::idx_t conv_time = 0;
//     start();
//     init_e();

//     comp1();
//     init_f();
//     run_time+= stop();
//     start();
//     permute_view(F, flayout01)();
//     conv_time += stop();
//     start();
//     comp2();
//     run_time += stop();
//     start();
//     permute_view(F, flayout10)();
//     conv_time += stop();
//     start();
//     init_g();
//     comp3();
//     run_time += stop();
//     write_datapoint(dimensionality, ps, variant, "Computation", run_time);
//     write_datapoint(dimensionality, ps, variant, "Conversion", conv_time);
//     return;
//   }

//   auto dec = format_decisions(tie(F), init_e, comp1, init_f, comp2, init_g, comp3);

//   if(numDecisions == -1) {
//     dec.set_format_for(F, flayout01, comp2);
//     dec.set_format_for(F, flayout10, comp3);
//   }
//   if(numDecisions >= 1) {
//     dec.set_format_for(F, flayout01, init_e);
//   }
//   if(numDecisions >= 2) {
//     dec.set_format_for(F, flayout01, comp1);
//   }
//   if(numDecisions >= 3) {
//     dec.set_format_for(F, flayout01, init_f);
//   }
//   if(numDecisions >= 4) {
//     dec.set_format_for(F, flayout01, comp2);
//   }
//   if(numDecisions >= 5) {
//     dec.set_format_for(F, flayout10, init_g);
//   }
//   if(numDecisions >= 6) {
//     dec.set_format_for(F, flayout10, comp3);
//   }
//   if(numDecisions >= 7) {
//     dec.set_output_format(F, flayout10);
//   }
  

//   if(variant.find("lock") != std::string::npos) {
//     dec.lock();
//   }
//   auto initializing_time = stop();
//   start();
//   auto comp = dec.finalize();
//   //auto solve_time = stop();
//   start();
//   //comp();
  

//   auto breakdown = dec.time_execution();
//   auto conversion_time = camp::get<0>(breakdown);
//   auto computation_time = camp::get<1>(breakdown);
//   auto kernel_time = stop();
//   write_datapoint(dimensionality, ps, variant, "Initialization", initializing_time);
//   write_datapoint(dimensionality, ps, variant, "Constraint Setup", dec.setup_time);
//   write_datapoint(dimensionality, ps, variant, "ISL Space Setup", dec.space_time);
//   write_datapoint(dimensionality, ps, variant, "ISL Map Setup", dec.map_time);
//   write_datapoint(dimensionality, ps, variant, "ISL Solve", dec.solve_time);
//   write_datapoint(dimensionality, ps, variant, "Kernels", kernel_time);
//   write_datapoint(dimensionality, ps, variant, "Conversion", conversion_time);
//   write_datapoint(dimensionality, ps, variant, "Computation", computation_time);
// }


// auto run_2d_3knl_1var_nDecisions(std::string variant, camp::idx_t numDecisions, std::size_t ps) {
  
//   std::string dimensionality = "2D1Var";

//   std::stringstream vv;
//   vv << variant << "_" << numDecisions;
//   variant = vv.str();
//   start();

//   using namespace RAJA;
//   using VIEW = View<double, Layout<2>>;
//   double root = std::pow(ps, 1.0/2.0);
//   camp::idx_t  n = (camp::idx_t) root;
//   VIEW A(new double[n*n], n,n);
//   VIEW B(new double[n*n], n,n);
//   VIEW C(new double[n*n], n,n);
//   VIEW D(new double[n*n], n,n);
//   VIEW E(new double[n*n], n,n);
//   VIEW F(new double[n*n], n,n);
//   VIEW G(new double[n*n], n,n);

//   auto lam_init_e = [&](auto i, auto j) {
//     E(i,j) = 0;
//   };
//   auto lam_init_f = [&](auto i, auto j) {
//     F(i,j) = 0;
//   };
//   auto lam_init_g = [&](auto i, auto j) {
//     G(i,j) = 0;
//   };

//   auto lam_comp1 = [&](auto i, auto j, auto k) {
//     E(i,j) += A(i,k) * B(k,j);
//   };
//   auto lam_comp2 = [&](auto i, auto j, auto k) {
//     F(i,j) += C(i,k) * D(k,j);
//   };
//   auto lam_comp3 = [&](auto i, auto j, auto k) {
//     G(i,j) += E(i,k) * F(k,j);
//   };

//   using namespace RAJA;
//   using INIT_POL = KernelPolicy<
//     statement::For<0, omp_parallel_for_exec,
//       statement::For<1, loop_exec,

//           statement::Lambda<0>
//         >

//     >
//   >;
    
//   using COMP_POL = KernelPolicy<
//     statement::For<0, omp_parallel_for_exec,
//       statement::For<1, loop_exec,
//         statement::For<2, loop_exec,

//             statement::Lambda<0>
//         >
//       >
//     >
//   >;
//   auto eseg = make_tuple(RangeSegment(0,n), RangeSegment(0,n));
//   auto fseg = make_tuple(RangeSegment(0,n), RangeSegment(0,n));
//   auto gseg = make_tuple(RangeSegment(0,n), RangeSegment(0,n));

//   auto comp_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));

//   auto init_e = make_kernel<INIT_POL>(eseg, lam_init_e);
//   auto comp1 = make_kernel<COMP_POL>(comp_seg, lam_comp1);
//   auto init_f = make_kernel<INIT_POL>(fseg, lam_init_f);
//   auto comp2 = make_kernel<COMP_POL>(comp_seg, lam_comp2);
//   auto init_g = make_kernel<INIT_POL>(gseg, lam_init_g);
//   auto comp3 = make_kernel<COMP_POL>(comp_seg, lam_comp3);

//   init_e();
//   init_g();
//   init_f();
//   if (variant.find("original") != std::string::npos) {
    
//     comp1();
//     comp2();
//     comp3();
    
   
//     auto run_time = stop();
//     write_datapoint(dimensionality, ps, variant, "Computation", run_time);
//     return;
//   }


//   std::array<idx_t, 2> bsizes{{n,n}};
//   auto blayout = make_permuted_layout(bsizes, {{1,0}});
//   std::array<idx_t, 2> dsizes{{n,n}};
//   auto dlayout = make_permuted_layout(dsizes, {{1,0}});
//   std::array<idx_t, 2> fsizes{{n,n}};
//   auto flayout01 = make_permuted_layout(fsizes, {{0,1}});
//   auto flayout10 = make_permuted_layout(fsizes, {{1,0}});

//   if (variant.find("hand") != std::string::npos) {
//     camp::idx_t run_time = 0;
//     camp::idx_t conv_time = 0;
//     start();
//     init_e();

//     comp1();
//     init_f();
//     run_time+= stop();
//     start();
//     permute_view(F, flayout01)();
//     conv_time += stop();
//     start();
//     comp2();
//     run_time += stop();
//     start();
//     permute_view(F, flayout10)();
//     conv_time += stop();
//     start();
//     init_g();
//     comp3();
//     run_time += stop();
//     write_datapoint(dimensionality, ps, variant, "Computation", run_time);
//     write_datapoint(dimensionality, ps, variant, "Conversion", conv_time);
//     return;
//   }

//   auto dec = format_decisions(tie(F), init_e, comp1, init_f, comp2, init_g, comp3);

//   if(numDecisions == -1) {
//     dec.set_format_for(F, flayout01, comp2);
//     dec.set_format_for(F, flayout10, comp3);
//   }
//   if(numDecisions >= 1) {
//     dec.set_format_for(F, flayout01, init_e);
//   }
//   if(numDecisions >= 2) {
//     dec.set_format_for(F, flayout01, comp1);
//   }
//   if(numDecisions >= 3) {
//     dec.set_format_for(F, flayout01, init_f);
//   }
//   if(numDecisions >= 4) {
//     dec.set_format_for(F, flayout01, comp2);
//   }
//   if(numDecisions >= 5) {
//     dec.set_format_for(F, flayout10, init_g);
//   }
//   if(numDecisions >= 6) {
//     dec.set_format_for(F, flayout10, comp3);
//   }
//   if(numDecisions >= 7) {
//     dec.set_output_format(F, flayout10);
//   }
  

//   if(variant.find("lock") != std::string::npos) {
//     dec.lock();
//   }
//   auto initializing_time = stop();
//   start();
//   auto comp = dec.finalize();
//   //auto solve_time = stop();
//   start();
//   //comp();
  

//   auto breakdown = dec.time_execution();
//   auto conversion_time = camp::get<0>(breakdown);
//   auto computation_time = camp::get<1>(breakdown);
//   auto kernel_time = stop();
//   write_datapoint(dimensionality, ps, variant, "Initialization", initializing_time);
//   write_datapoint(dimensionality, ps, variant, "Constraint Setup", dec.setup_time);
//   write_datapoint(dimensionality, ps, variant, "ISL Space Setup", dec.space_time);
//   write_datapoint(dimensionality, ps, variant, "ISL Map Setup", dec.map_time);
//   write_datapoint(dimensionality, ps, variant, "ISL Solve", dec.solve_time);
//   write_datapoint(dimensionality, ps, variant, "Kernels", kernel_time);
//   write_datapoint(dimensionality, ps, variant, "Conversion", conversion_time);
//   write_datapoint(dimensionality, ps, variant, "Computation", computation_time);
// }
// // Copyright (c) 2016-19, Lawrence Livermore National Security, LLC
// // and RAJA project contributors. See the RAJA/COPYRIGHT file for details.
// //
// // SPDX-License-Identifier: (BSD-3-Clause)
// //////////////////////////////////////////////////////////////////////////////
// #include "RAJA/RAJA.hpp"
// #include <algorithm>
// #include <chrono>
// using namespace RAJA;
// template <idx_t I0, idx_t I1, idx_t I2>
// struct order_to_kpol3 {
//   using Policy = KernelPolicy<
//     statement::For<I0, omp_parallel_for_exec,
//       statement::For<I1, loop_exec,
//         statement::For<I2, loop_exec,
//           statement::Lambda<0>
//         >
//       >
//     >
//   >;
// };

// void write_datapoint(auto dimensionality, auto problemSize, auto views, auto constraints, auto variant, auto component, auto time) {
//   std::cout << dimensionality << "," << problemSize << "," << views << "," << constraints << "," << variant << "," << component << "," << time << "\n";
// }

// std::chrono::time_point<std::chrono::high_resolution_clock> _start;
// void start() {
//    _start = std::chrono::high_resolution_clock::now();
// }

// auto stop() {
//   auto end = std::chrono::high_resolution_clock::now();

//   auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - _start).count();
//   return duration;
// }

// auto run_3d_3knl_1var_nDecisions(std::string variant, camp::idx_t numDecisions, std::size_t ps) {
  
//   std::string dimensionality = "3D1Var";

//   std::stringstream vv;
//   vv << variant << "_" << numDecisions;
//   variant = vv.str();
//   start();

//   using namespace RAJA;
//   using VIEW = View<double, Layout<3>>;
//   double root = std::pow(ps, 1.0/3.0);
//   camp::idx_t  n = (camp::idx_t) root;
//   VIEW A(new double[n*n*n], n,n,n);
//   VIEW B(new double[n*n*n], n,n,n);
//   VIEW C(new double[n*n*n], n,n,n);
//   VIEW D(new double[n*n*n], n,n,n);
//   VIEW E(new double[n*n*n], n,n,n);
//   VIEW F(new double[n*n*n], n,n,n);
//   VIEW G(new double[n*n*n], n,n,n);

//   auto lam_init_e = [&](auto i, auto j, auto k) {
//     E(i,j,k) = 0;
//   };
//   auto lam_init_f = [&](auto i, auto j, auto k) {
//     F(i,j,k) = 0;
//   };
//   auto lam_init_g = [&](auto i, auto j, auto k) {
//     G(i,j,k) = 0;
//   };

//   auto lam_comp1 = [&](auto l, auto i, auto j, auto k) {
//     E(l,i,j) += A(l,i,k) * B(l,k,j);
//   };
//   auto lam_comp2 = [&](auto l, auto i, auto j, auto k) {
//     F(l,i,j) += C(l,i,k) * D(l,k,j);
//   };
//   auto lam_comp3 = [&](auto l, auto i, auto j, auto k) {
//     G(l,i,j) += E(l,i,k) * F(l,k,j);
//   };

//   using namespace RAJA;
//   using INIT_POL = KernelPolicy<
//     statement::For<0, omp_parallel_for_exec,
//       statement::For<1, loop_exec,
//         statement::For<2, loop_exec,
//           statement::Lambda<0>
//         >
//       >
//     >
//   >;
    
//   using COMP_POL = KernelPolicy<
//     statement::For<0, omp_parallel_for_exec,
//       statement::For<1, loop_exec,
//         statement::For<2, loop_exec,
//           statement::For<3, loop_exec,
//             statement::Lambda<0>
//           >
//         >
//       >
//     >
//   >;
//   auto eseg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));
//   auto fseg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));
//   auto gseg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));

//   auto comp_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));

//   auto init_e = make_kernel<INIT_POL>(eseg, lam_init_e);
//   auto comp1 = make_kernel<COMP_POL>(comp_seg, lam_comp1);
//   auto init_f = make_kernel<INIT_POL>(fseg, lam_init_f);
//   auto comp2 = make_kernel<COMP_POL>(comp_seg, lam_comp2);
//   auto init_g = make_kernel<INIT_POL>(gseg, lam_init_g);
//   auto comp3 = make_kernel<COMP_POL>(comp_seg, lam_comp3);

//   if (variant.find("original") != std::string::npos) {
//     init_e();
//     comp1();
//     init_f();
//     comp2();
//     init_g();
//     comp3();
//     auto run_time = stop();
//     write_datapoint(dimensionality, ps, variant, "Computation", run_time);
//     return;
//   }


//   std::array<idx_t, 3> bsizes{{n,n,n}};
//   auto blayout = make_permuted_layout(bsizes, {{2,1,0}});
//   std::array<idx_t, 3> dsizes{{n,n,n}};
//   auto dlayout = make_permuted_layout(dsizes, {{2,1,0}});
//   std::array<idx_t, 3> fsizes{{n,n,n}};
//   auto flayout01 = make_permuted_layout(fsizes, {{2,0,1}});
//   auto flayout10 = make_permuted_layout(fsizes, {{2,1,0}});

//   if (variant.find("hand") != std::string::npos) {
//     camp::idx_t run_time = 0;
//     camp::idx_t conv_time = 0;
//     start();
//     init_e();

//     comp1();
//     init_f();
//     run_time+= stop();
//     start();
//     permute_view(F, flayout01)();
//     conv_time += stop();
//     start();
//     comp2();
//     run_time += stop();
//     start();
//     permute_view(F, flayout10)();
//     conv_time += stop();
//     start();
//     init_g();
//     comp3();
//     run_time += stop();
//     write_datapoint(dimensionality, ps, variant, "Computation", run_time);
//     write_datapoint(dimensionality, ps, variant, "Conversion", conv_time);
//     return;
//   }

//   auto dec = format_decisions(tie(F), init_e, comp1, init_f, comp2, init_g, comp3);

//   if(numDecisions == -1) {
//     dec.set_format_for(F, flayout01, comp2);
//     dec.set_format_for(F, flayout10, comp3);
//   }
//   if(numDecisions >= 1) {
//     dec.set_format_for(F, flayout01, init_e);
//   }
//   if(numDecisions >= 2) {
//     dec.set_format_for(F, flayout01, comp1);
//   }
//   if(numDecisions >= 3) {
//     dec.set_format_for(F, flayout01, init_f);
//   }
//   if(numDecisions >= 4) {
//     dec.set_format_for(F, flayout01, comp2);
//   }
//   if(numDecisions >= 5) {
//     dec.set_format_for(F, flayout10, init_g);
//   }
//   if(numDecisions >= 6) {
//     dec.set_format_for(F, flayout10, comp3);
//   }
//   if(numDecisions >= 7) {
//     dec.set_output_format(F, flayout10);
//   }
  

//   if(variant.find("lock") != std::string::npos) {
//     dec.lock();
//   }
//   auto initializing_time = stop();
//   start();
//   auto comp = dec.finalize();
//   //auto solve_time = stop();
//   start();
//   //comp();
  

//   auto breakdown = dec.time_execution();
//   auto conversion_time = camp::get<0>(breakdown);
//   auto computation_time = camp::get<1>(breakdown);
//   auto kernel_time = stop();
//   write_datapoint(dimensionality, ps, variant, "Initialization", initializing_time);
//   write_datapoint(dimensionality, ps, variant, "Constraint Setup", dec.setup_time);
//   write_datapoint(dimensionality, ps, variant, "ISL Space Setup", dec.space_time);
//   write_datapoint(dimensionality, ps, variant, "ISL Map Setup", dec.map_time);
//   write_datapoint(dimensionality, ps, variant, "ISL Solve", dec.solve_time);
//   write_datapoint(dimensionality, ps, variant, "Kernels", kernel_time);
//   write_datapoint(dimensionality, ps, variant, "Conversion", conversion_time);
//   write_datapoint(dimensionality, ps, variant, "Computation", computation_time);
// }


// auto run_2d_3knl_1var_nDecisions(std::string variant, camp::idx_t numDecisions, std::size_t ps) {
  
//   std::string dimensionality = "2D1Var";

//   std::stringstream vv;
//   vv << variant << "_" << numDecisions;
//   variant = vv.str();
//   start();

//   using namespace RAJA;
//   j
//   double root = std::pow(ps, 1.0/2.0);
//   camp::idx_t  n = (camp::idx_t) root;
//   VIEW A(new double[n*n], n,n);
//   VIEW B(new double[n*n], n,n);
//   VIEW C(new double[n*n], n,n);
//   VIEW D(new double[n*n], n,n);
//   VIEW E(new double[n*n], n,n);
//   VIEW F(new double[n*n], n,n);
//   VIEW G(new double[n*n], n,n);

//   auto lam_init_e = [&](auto i, auto j) {
//     E(i,j) = 0;
//   };
//   auto lam_init_f = [&](auto i, auto j) {
//     F(i,j) = 0;
//   };
//   auto lam_init_g = [&](auto i, auto j) {
//     G(i,j) = 0;
//   };

//   auto lam_comp1 = [&](auto i, auto j, auto k) {
//     E(i,j) += A(i,k) * B(k,j);
//   };
//   auto lam_comp2 = [&](auto i, auto j, auto k) {
//     F(i,j) += C(i,k) * D(k,j);
//   };
//   auto lam_comp3 = [&](auto i, auto j, auto k) {
//     G(i,j) += E(i,k) * F(k,j);
//   };

//   using namespace RAJA;
//   using INIT_POL = KernelPolicy<
//     statement::For<0, omp_parallel_for_exec,
//       statement::For<1, loop_exec,

//           statement::Lambda<0>
//         >

//     >
//   >;
    
//   using COMP_POL = KernelPolicy<
//     statement::For<0, omp_parallel_for_exec,
//       statement::For<1, loop_exec,
//         statement::For<2, loop_exec,

//             statement::Lambda<0>
//         >
//       >
//     >
//   >;
//   auto eseg = make_tuple(RangeSegment(0,n), RangeSegment(0,n));
//   auto fseg = make_tuple(RangeSegment(0,n), RangeSegment(0,n));
//   auto gseg = make_tuple(RangeSegment(0,n), RangeSegment(0,n));

//   auto comp_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));

//   auto init_e = make_kernel<INIT_POL>(eseg, lam_init_e);
//   auto comp1 = make_kernel<COMP_POL>(comp_seg, lam_comp1);
//   auto init_f = make_kernel<INIT_POL>(fseg, lam_init_f);
//   auto comp2 = make_kernel<COMP_POL>(comp_seg, lam_comp2);
//   auto init_g = make_kernel<INIT_POL>(gseg, lam_init_g);
//   auto comp3 = make_kernel<COMP_POL>(comp_seg, lam_comp3);

//   init_e();
//   init_g();
//   init_f();
//   if (variant.find("original") != std::string::npos) {
    
//     comp1();
//     comp2();
//     comp3();
    
   
//     auto run_time = stop();
//     write_datapoint(dimensionality, ps, variant, "Computation", run_time);
//     return;
//   }


//   std::array<idx_t, 2> bsizes{{n,n}};
//   auto blayout = make_permuted_layout(bsizes, {{1,0}});
//   std::array<idx_t, 2> dsizes{{n,n}};
//   auto dlayout = make_permuted_layout(dsizes, {{1,0}});
//   std::array<idx_t, 2> fsizes{{n,n}};
//   auto flayout01 = make_permuted_layout(fsizes, {{0,1}});
//   auto flayout10 = make_permuted_layout(fsizes, {{1,0}});

//   if (variant.find("hand") != std::string::npos) {
//     camp::idx_t run_time = 0;
//     camp::idx_t conv_time = 0;
//     start();
//     init_e();

//     comp1();
//     init_f();
//     run_time+= stop();
//     start();
//     permute_view(F, flayout01)();
//     conv_time += stop();
//     start();
//     comp2();
//     run_time += stop();
//     start();
//     permute_view(F, flayout10)();
//     conv_time += stop();
//     start();
//     init_g();
//     comp3();
//     run_time += stop();
//     write_datapoint(dimensionality, ps, variant, "Computation", run_time);
//     write_datapoint(dimensionality, ps, variant, "Conversion", conv_time);
//     return;
//   }

//   auto dec = format_decisions(tie(F), init_e, comp1, init_f, comp2, init_g, comp3);

//   if(numDecisions == -1) {
//     dec.set_format_for(F, flayout01, comp2);
//     dec.set_format_for(F, flayout10, comp3);
//   }
//   if(numDecisions >= 1) {
//     dec.set_format_for(F, flayout01, init_e);
//   }
//   if(numDecisions >= 2) {
//     dec.set_format_for(F, flayout01, comp1);
//   }
//   if(numDecisions >= 3) {
//     dec.set_format_for(F, flayout01, init_f);
//   }
//   if(numDecisions >= 4) {
//     dec.set_format_for(F, flayout01, comp2);
//   }
//   if(numDecisions >= 5) {
//     dec.set_format_for(F, flayout10, init_g);
//   }
//   if(numDecisions >= 6) {
//     dec.set_format_for(F, flayout10, comp3);
//   }
//   if(numDecisions >= 7) {
//     dec.set_output_format(F, flayout10);
//   }
  

//   if(variant.find("lock") != std::string::npos) {
//     dec.lock();
//   }
//   auto initializing_time = stop();
//   start();
//   auto comp = dec.finalize();
//   //auto solve_time = stop();
//   start();
//   //comp();
  

//   auto breakdown = dec.time_execution();
//   auto conversion_time = camp::get<0>(breakdown);
//   auto computation_time = camp::get<1>(breakdown);
//   auto kernel_time = stop();
//   write_datapoint(dimensionality, ps, variant, "Initialization", initializing_time);
//   write_datapoint(dimensionality, ps, variant, "Constraint Setup", dec.setup_time);
//   write_datapoint(dimensionality, ps, variant, "ISL Space Setup", dec.space_time);
//   write_datapoint(dimensionality, ps, variant, "ISL Map Setup", dec.map_time);
//   write_datapoint(dimensionality, ps, variant, "ISL Solve", dec.solve_time);
//   write_datapoint(dimensionality, ps, variant, "Kernels", kernel_time);
//   write_datapoint(dimensionality, ps, variant, "Conversion", conversion_time);
//   write_datapoint(dimensionality, ps, variant, "Computation", computation_time);
// }
// // Copyright (c) 2016-19, Lawrence Livermore National Security, LLC
// // and RAJA project contributors. See the RAJA/COPYRIGHT file for details.
// //
// // SPDX-License-Identifier: (BSD-3-Clause)
// //////////////////////////////////////////////////////////////////////////////
// #include "RAJA/RAJA.hpp"
// #include <algorithm>
// #include <chrono>
// using namespace RAJA;
// template <idx_t I0, idx_t I1, idx_t I2>
// struct order_to_kpol3 {
//   using Policy = KernelPolicy<
//     statement::For<I0, omp_parallel_for_exec,
//       statement::For<I1, loop_exec,
//         statement::For<I2, loop_exec,
//           statement::Lambda<0>
//         >
//       >
//     >
//   >;
// };

// void write_datapoint(auto dimensionality, auto problemSize, auto views, auto constraints, auto variant, auto component, auto time) {
//   std::cout << dimensionality << "," << problemSize << "," << views << "," << constraints << "," << variant << "," << component << "," << time << "\n";
// }

// std::chrono::time_point<std::chrono::high_resolution_clock> _start;
// void start() {
//    _start = std::chrono::high_resolution_clock::now();
// }

// auto stop() {
//   auto end = std::chrono::high_resolution_clock::now();

//   auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - _start).count();
//   return duration;
// }

// auto run_3d_3knl_1var_nDecisions(std::string variant, camp::idx_t numDecisions, std::size_t ps) {
  
//   std::string dimensionality = "3D1Var";

//   std::stringstream vv;
//   vv << variant << "_" << numDecisions;
//   variant = vv.str();
//   start();

//   using namespace RAJA;
//   using VIEW = View<double, Layout<3>>;
//   double root = std::pow(ps, 1.0/3.0);
//   camp::idx_t  n = (camp::idx_t) root;
//   VIEW A(new double[n*n*n], n,n,n);
//   VIEW B(new double[n*n*n], n,n,n);
//   VIEW C(new double[n*n*n], n,n,n);
//   VIEW D(new double[n*n*n], n,n,n);
//   VIEW E(new double[n*n*n], n,n,n);
//   VIEW F(new double[n*n*n], n,n,n);
//   VIEW G(new double[n*n*n], n,n,n);

//   auto lam_init_e = [&](auto i, auto j, auto k) {
//     E(i,j,k) = 0;
//   };
//   auto lam_init_f = [&](auto i, auto j, auto k) {
//     F(i,j,k) = 0;
//   };
//   auto lam_init_g = [&](auto i, auto j, auto k) {
//     G(i,j,k) = 0;
//   };

//   auto lam_comp1 = [&](auto l, auto i, auto j, auto k) {
//     E(l,i,j) += A(l,i,k) * B(l,k,j);
//   };
//   auto lam_comp2 = [&](auto l, auto i, auto j, auto k) {
//     F(l,i,j) += C(l,i,k) * D(l,k,j);
//   };
//   auto lam_comp3 = [&](auto l, auto i, auto j, auto k) {
//     G(l,i,j) += E(l,i,k) * F(l,k,j);
//   };

//   using namespace RAJA;
//   using INIT_POL = KernelPolicy<
//     statement::For<0, omp_parallel_for_exec,
//       statement::For<1, loop_exec,
//         statement::For<2, loop_exec,
//           statement::Lambda<0>
//         >
//       >
//     >
//   >;
    
//   using COMP_POL = KernelPolicy<
//     statement::For<0, omp_parallel_for_exec,
//       statement::For<1, loop_exec,
//         statement::For<2, loop_exec,
//           statement::For<3, loop_exec,
//             statement::Lambda<0>
//           >
//         >
//       >
//     >
//   >;
//   auto eseg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));
//   auto fseg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));
//   auto gseg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));

//   auto comp_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));

//   auto init_e = make_kernel<INIT_POL>(eseg, lam_init_e);
//   auto comp1 = make_kernel<COMP_POL>(comp_seg, lam_comp1);
//   auto init_f = make_kernel<INIT_POL>(fseg, lam_init_f);
//   auto comp2 = make_kernel<COMP_POL>(comp_seg, lam_comp2);
//   auto init_g = make_kernel<INIT_POL>(gseg, lam_init_g);
//   auto comp3 = make_kernel<COMP_POL>(comp_seg, lam_comp3);

//   if (variant.find("original") != std::string::npos) {
//     init_e();
//     comp1();
//     init_f();
//     comp2();
//     init_g();
//     comp3();
//     auto run_time = stop();
//     write_datapoint(dimensionality, ps, variant, "Computation", run_time);
//     return;
//   }


//   std::array<idx_t, 3> bsizes{{n,n,n}};
//   auto blayout = make_permuted_layout(bsizes, {{2,1,0}});
//   std::array<idx_t, 3> dsizes{{n,n,n}};
//   auto dlayout = make_permuted_layout(dsizes, {{2,1,0}});
//   std::array<idx_t, 3> fsizes{{n,n,n}};
//   auto flayout01 = make_permuted_layout(fsizes, {{2,0,1}});
//   auto flayout10 = make_permuted_layout(fsizes, {{2,1,0}});

//   if (variant.find("hand") != std::string::npos) {
//     camp::idx_t run_time = 0;
//     camp::idx_t conv_time = 0;
//     start();
//     init_e();

//     comp1();
//     init_f();
//     run_time+= stop();
//     start();
//     permute_view(F, flayout01)();
//     conv_time += stop();
//     start();
//     comp2();
//     run_time += stop();
//     start();
//     permute_view(F, flayout10)();
//     conv_time += stop();
//     start();
//     init_g();
//     comp3();
//     run_time += stop();
//     write_datapoint(dimensionality, ps, variant, "Computation", run_time);
//     write_datapoint(dimensionality, ps, variant, "Conversion", conv_time);
//     return;
//   }

//   auto dec = format_decisions(tie(F), init_e, comp1, init_f, comp2, init_g, comp3);

//   if(numDecisions == -1) {
//     dec.set_format_for(F, flayout01, comp2);
//     dec.set_format_for(F, flayout10, comp3);
//   }
//   if(numDecisions >= 1) {
//     dec.set_format_for(F, flayout01, init_e);
//   }
//   if(numDecisions >= 2) {
//     dec.set_format_for(F, flayout01, comp1);
//   }
//   if(numDecisions >= 3) {
//     dec.set_format_for(F, flayout01, init_f);
//   }
//   if(numDecisions >= 4) {
//     dec.set_format_for(F, flayout01, comp2);
//   }
//   if(numDecisions >= 5) {
//     dec.set_format_for(F, flayout10, init_g);
//   }
//   if(numDecisions >= 6) {
//     dec.set_format_for(F, flayout10, comp3);
//   }
//   if(numDecisions >= 7) {
//     dec.set_output_format(F, flayout10);
//   }
  

//   if(variant.find("lock") != std::string::npos) {
//     dec.lock();
//   }
//   auto initializing_time = stop();
//   start();
//   auto comp = dec.finalize();
//   //auto solve_time = stop();
//   start();
//   //comp();
  

//   auto breakdown = dec.time_execution();
//   auto conversion_time = camp::get<0>(breakdown);
//   auto computation_time = camp::get<1>(breakdown);
//   auto kernel_time = stop();
//   write_datapoint(dimensionality, ps, variant, "Initialization", initializing_time);
//   write_datapoint(dimensionality, ps, variant, "Constraint Setup", dec.setup_time);
//   write_datapoint(dimensionality, ps, variant, "ISL Space Setup", dec.space_time);
//   write_datapoint(dimensionality, ps, variant, "ISL Map Setup", dec.map_time);
//   write_datapoint(dimensionality, ps, variant, "ISL Solve", dec.solve_time);
//   write_datapoint(dimensionality, ps, variant, "Kernels", kernel_time);
//   write_datapoint(dimensionality, ps, variant, "Conversion", conversion_time);
//   write_datapoint(dimensionality, ps, variant, "Computation", computation_time);
// }




void experiment1(idx_t problemSize, bool quiet=false) {
  idx_t dimensionality = 2;
  idx_t computation=3;
  idx_t views=3;
  std::string constraints="BestAnalytical";

  double root = std::pow(problemSize, 1.0/dimensionality);
  camp::idx_t n = (camp::idx_t) root;

  //BaseRAJA definitions
  using VIEW = View<double, Layout<2>>;
  VIEW A(new double[n*n], n,n);
  VIEW B(new double[n*n], n,n);
  VIEW C(new double[n*n], n,n);
  VIEW D(new double[n*n], n,n);
  VIEW E(new double[n*n], n,n);
  VIEW F(new double[n*n], n,n);
  VIEW G(new double[n*n], n,n);

  auto lam_init_e = [&](auto i, auto j) {
    E(i,j) = rand();
  };
  auto lam_init_f = [&](auto i, auto j) {
    F(i,j) = rand();
  };
  auto lam_init_g = [&](auto i, auto j) {
    G(i,j) = rand();
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


  auto init_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n));
  auto comp_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));

  auto init_e = make_kernel<INIT_POL>(init_seg, lam_init_e);
  auto init_f = make_kernel<INIT_POL>(init_seg, lam_init_f);
  auto init_g = make_kernel<INIT_POL>(init_seg, lam_init_g);
  
  auto comp1 = make_kernel<COMP_POL>(comp_seg, lam_comp1);
  auto comp2 = make_kernel<COMP_POL>(comp_seg, lam_comp2);
  auto comp3 = make_kernel<COMP_POL>(comp_seg, lam_comp3);

  std::array<idx_t, 2> sizes{{n,n}};
  auto blayout = make_permuted_layout(sizes, {{1,0}});
  auto dlayout = make_permuted_layout(sizes, {{1,0}});
  auto flayout01 = make_permuted_layout(sizes, {{0,1}});
  auto flayout10 = make_permuted_layout(sizes, {{1,0}});

  init_e();
  init_f();
  init_g();
  //BaseRAJA performance
  start();
  comp1();
  comp2();
  comp3();
  auto baseRAJA = stop();
  if(!quiet) {
    write_datapoint(1, dimensionality, problemSize, views, constraints, "BaseRAJA", "Computation", baseRAJA);
  }

  //HandTransformed
  permute_view(F, flayout01)();
  permute_view(B, flayout01)();
  permute_view(D, flayout01)();
  init_e();
  init_f();
  init_g();


  camp::idx_t run_time = 0;
  camp::idx_t conv_time = 0;
  start();
    permute_view(B, flayout10)();
  conv_time += stop();
  start();
    comp1();
  run_time += stop();
  start();
    permute_view(F, flayout01)();
    permute_view(D, flayout10)();
  conv_time += stop();
  start();
    comp2();
  run_time += stop();
  start();
    permute_view(F, flayout10)();
  conv_time += stop();
  start();
    comp3();
  run_time += stop();
  if(!quiet) {
  write_datapoint(1, dimensionality, problemSize, views, constraints, "HandTransformed", "Computation", run_time);
  write_datapoint(1, dimensionality, problemSize, views, constraints, "HandTransformed", "Conversion", conv_time);
  }

  //Experiment1
  permute_view(F, flayout01)();
  permute_view(B, flayout01)();
  permute_view(D, flayout01)();
  init_e();
  init_f();
  init_g();

  auto dec = format_decisions(tie(B,D,F), comp1, comp2, comp3);
  dec.set_format_for(B,flayout10, comp1);
  dec.set_format_for(D,flayout10, comp2);
  dec.set_format_for(F,flayout01, comp2);
  dec.set_format_for(F,flayout10, comp3);

  auto breakdown = dec.time_execution();
  auto conversion_time = get<0>(breakdown);
  auto computation_time = get<1>(breakdown);

  if(!quiet) {
  write_datapoint(1, dimensionality, problemSize, views, constraints, "Experiment1", "Computation", computation_time);
  write_datapoint(1, dimensionality, problemSize, views, constraints, "Experiment1", "Conversion", conversion_time);
  write_datapoint(1, dimensionality, problemSize, views, constraints, "Experiment1", "Cost Estimation", dec.setup_time);
  write_datapoint(1, dimensionality, problemSize, views, constraints, "Experiment1", "ISL Space", dec.space_time);
  write_datapoint(1, dimensionality, problemSize, views, constraints, "Experiment1", "ISL Func", dec.map_time);
  write_datapoint(1, dimensionality, problemSize, views, constraints, "Experiment1", "ISL Solve", dec.solve_time);
  }
}


void experiment2(idx_t problemSize, bool quiet=false) {
  idx_t dimensionality = 2;
  idx_t computation=3;
  idx_t views=3;
  std::string constraints="BestAnalytical";

  double root = std::pow(problemSize, 1.0/dimensionality);
  camp::idx_t n = (camp::idx_t) root;

  //BaseRAJA definitions
  using VIEW = View<double, Layout<2>>;
  VIEW A(new double[n*n], n,n);
  VIEW B(new double[n*n], n,n);
  VIEW C(new double[n*n], n,n);
  VIEW D(new double[n*n], n,n);
  VIEW E(new double[n*n], n,n);
  VIEW F(new double[n*n], n,n);
  VIEW G(new double[n*n], n,n);

  auto lam_init_e = [&](auto i, auto j) {
    E(i,j) = rand();
  };
  auto lam_init_f = [&](auto i, auto j) {
    F(i,j) = rand();
  };
  auto lam_init_g = [&](auto i, auto j) {
    G(i,j) = rand();
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


  auto init_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n));
  auto comp_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));

  auto init_e = make_kernel<INIT_POL>(init_seg, lam_init_e);
  auto init_f = make_kernel<INIT_POL>(init_seg, lam_init_f);
  auto init_g = make_kernel<INIT_POL>(init_seg, lam_init_g);
  
  auto comp1 = make_kernel<COMP_POL>(comp_seg, lam_comp1);
  auto comp2 = make_kernel<COMP_POL>(comp_seg, lam_comp2);
  auto comp3 = make_kernel<COMP_POL>(comp_seg, lam_comp3);

  std::array<idx_t, 2> sizes{{n,n}};
  auto blayout = make_permuted_layout(sizes, {{1,0}});
  auto dlayout = make_permuted_layout(sizes, {{1,0}});
  auto flayout01 = make_permuted_layout(sizes, {{0,1}});
  auto flayout10 = make_permuted_layout(sizes, {{1,0}});

  init_e();
  init_f();
  init_g();
  //BaseRAJA performance

  comp1();
  comp2();
  comp3();

}
/*
void experiment2(idx_t constraints, bool quiet=false) {
  idx_t dimensionality = 2;
  idx_t computation=3;
  idx_t views = 1;
  idx_t problemSize = 1000000; //10^6

  double root = std::pow(problemSize, 1.0/dimensionality);
  camp::idx_t n = (camp::idx_t) root;

  //BaseRAJA definitions
  using VIEW = View<double, Layout<2>>;
  VIEW A(new double[n*n], n,n);
  VIEW B(new double[n*n], n,n);
  VIEW C(new double[n*n], n,n);
  VIEW D(new double[n*n], n,n);
  VIEW E(new double[n*n], n,n);
  VIEW F(new double[n*n], n,n);
  VIEW G(new double[n*n], n,n);

  auto lam_init_e = [&](auto i, auto j) {
    E(i,j) = rand();
  };
  auto lam_init_f = [&](auto i, auto j) {
    F(i,j) = rand();
  };
  auto lam_init_g = [&](auto i, auto j) {
    G(i,j) = rand();
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


  auto init_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n));
  auto comp_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));

  auto init_e = make_kernel<INIT_POL>(init_seg, lam_init_e);
  auto init_f = make_kernel<INIT_POL>(init_seg, lam_init_f);
  auto init_g = make_kernel<INIT_POL>(init_seg, lam_init_g);
  
  auto comp1 = make_kernel<COMP_POL>(comp_seg, lam_comp1);
  auto comp2 = make_kernel<COMP_POL>(comp_seg, lam_comp2);
  auto comp3 = make_kernel<COMP_POL>(comp_seg, lam_comp3);

  std::array<idx_t, 2> sizes{{n,n}};
  auto blayout = make_permuted_layout(sizes, {{1,0}});
  auto dlayout = make_permuted_layout(sizes, {{1,0}});
  auto flayout01 = make_permuted_layout(sizes, {{0,1}});
  auto flayout10 = make_permuted_layout(sizes, {{1,0}});

  init_e();
  init_f();
  init_g();
  //BaseRAJA performance
  start();
  comp1();
  comp2();
  comp3();
  auto baseRAJA = stop();
  if(!quiet) {
    write_datapoint(2,dimensionality, problemSize, views, constraints, "BaseRAJA", "Computation", baseRAJA);
  }

  //Experiment2
  permute_view(F, flayout01)();
  permute_view(B, flayout01)();
  permute_view(D, flayout01)();
  init_e();
  init_f();
  init_g();

  
  auto dec = format_decisions(tie(F), comp1, comp2, comp3);
  if(constraints >= 1) {
    dec.set_format_for(F,flayout01, comp1);
  }
  if(constraints >= 2) {
    dec.set_format_for(F,flayout01, comp2);
  }
  if(constraints >= 3) {
    dec.set_format_for(F,flayout10, comp3);
  }
  if(constraints >= 4) {
    dec.set_output_format(F,flayout10);
  }

  auto breakdown = dec.time_execution();
  auto conversion_time = get<0>(breakdown);
  auto computation_time = get<1>(breakdown);

  if(!quiet) {
  write_datapoint(2,dimensionality, problemSize, views, constraints, "Experiment2", "Computation", computation_time);
  write_datapoint(2,dimensionality, problemSize, views, constraints, "Experiment2", "Conversion", conversion_time);
  write_datapoint(2,dimensionality, problemSize, views, constraints, "Experiment2", "Cost Estimation", dec.setup_time);
  write_datapoint(2,dimensionality, problemSize, views, constraints, "Experiment2", "ISL Space", dec.space_time);
  write_datapoint(2,dimensionality, problemSize, views, constraints, "Experiment2", "ISL Func", dec.map_time);
  write_datapoint(2,dimensionality, problemSize, views, constraints, "Experiment2", "ISL Solve", dec.solve_time);
  }
}

void experiment4_3(bool quiet=false) {
  idx_t dimensionality = 3;
  idx_t computation=3;
  idx_t views = 3;
  idx_t Views = views;
  idx_t problemSize = std::pow(2,30); //2^24 
std::string constraints = "BestLayoutsAnalytical";
  using VIEW = View<double, Layout<3>>;
  double root = std::pow(problemSize, 1.0/4.0);
  camp::idx_t  n = (camp::idx_t) root;
  VIEW A(new double[n*n*n], n,n,n);
  VIEW B(new double[n*n*n], n,n,n);
  VIEW C(new double[n*n*n], n,n,n);
  VIEW D(new double[n*n*n], n,n,n);
  VIEW E(new double[n*n*n], n,n,n);
  VIEW F(new double[n*n*n], n,n,n);
  VIEW G(new double[n*n*n], n,n,n);

  auto lam_init_e = [&](auto i, auto j, auto k) {
    E(i,j,k) = std::rand();
  };
  auto lam_init_f = [&](auto i, auto j, auto k) {
    F(i,j,k) = std::rand();
  };
  auto lam_init_g = [&](auto i, auto j, auto k) {
    G(i,j,k) = std::rand();
  };

  auto lam_comp1 = [&](auto l, auto i, auto j, auto k) {
    E(l,i,j) += A(l,i,k) * B(l,k,j);
  };
  auto lam_comp2 = [&](auto l, auto i, auto j, auto k) {
    F(l,i,j) += C(l,i,k) * D(l,k,j);
  };
  auto lam_comp3 = [&](auto l, auto i, auto j, auto k) {
    G(l,i,j) += E(l,i,k) * F(l,k,j);
  };

  using namespace RAJA;
  using INIT_POL = KernelPolicy<
    statement::For<0, omp_parallel_for_exec,
      statement::For<1, loop_exec,
        statement::For<2, loop_exec,
          statement::Lambda<0>
        >
      >
    >
  >;
  
  using COMP_POL = KernelPolicy<
    statement::For<0, omp_parallel_for_exec,
      statement::For<1, loop_exec,
        statement::For<2, loop_exec,
          statement::For<3, loop_exec,
            statement::Lambda<0>
          >
        >
      >
    >
  >;
  auto init_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));
  auto comp_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));

  auto init_e = make_kernel<INIT_POL>(init_seg, lam_init_e);
  auto init_f = make_kernel<INIT_POL>(init_seg, lam_init_f);
  auto init_g = make_kernel<INIT_POL>(init_seg, lam_init_g);
  
  auto comp1 = make_kernel<COMP_POL>(comp_seg, lam_comp1);
  auto comp2 = make_kernel<COMP_POL>(comp_seg, lam_comp2);
  auto comp3 = make_kernel<COMP_POL>(comp_seg, lam_comp3);

  std::array<idx_t, 3> sizes{{n,n,n}};
  auto blayout = make_permuted_layout(sizes, {{0,2,1}});
  auto flayout01 = make_permuted_layout(sizes, {{0,1,2}});
  auto flayout10 = make_permuted_layout(sizes, {{0,2,1}});


  
  //BaseRAJA performance
  init_e();
  init_f();
  init_g();
  start();
  comp1();
  comp2();
  comp3();
  auto baseRAJA = stop();
  if(!quiet) {
    write_datapoint(4, dimensionality, problemSize, views, constraints, "BaseRAJA", "Computation", baseRAJA);
  }

  //HandTransformed
  permute_view(F, flayout01)();
  permute_view(B, flayout01)();
  permute_view(D, flayout01)();
  init_e();
  init_f();
  init_g();

  camp::idx_t run_time = 0;
  camp::idx_t conv_time = 0;
  
  start();
    permute_view(B, flayout10)();
  conv_time += stop();
  start();
    comp1();
  run_time += stop();
  start();
    permute_view(D, flayout10)();
    permute_view(F, flayout01)();
  conv_time += stop();
  start();
    comp2();
  run_time += stop();
  start();
    permute_view(F, flayout10)();
  conv_time += stop();
  start();
    comp3();
  run_time += stop();
  if(!quiet) {
  write_datapoint(4, dimensionality, problemSize, views, constraints, "HandTransformed", "Computation", run_time);
  write_datapoint(4, dimensionality, problemSize, views, constraints, "HandTransformed", "Conversion", conv_time);
  }
  
  //Experiment 
  permute_view(F, flayout01)();
  permute_view(B, flayout01)();
  permute_view(D, flayout01)();
  init_e();
  init_f();
  init_g();

  auto dec = format_decisions(tie(B,D,F), comp1, comp2, comp3);
  dec.set_format_for(B,flayout10, comp1);
  dec.set_format_for(D,flayout10, comp2);
  dec.set_format_for(F,flayout01, comp2);
  dec.set_format_for(F,flayout10, comp3);
  auto breakdown = dec.time_execution();
  auto conversion_time = get<0>(breakdown);
  auto computation_time = get<1>(breakdown);
  if(!quiet) {
  write_datapoint(4, dimensionality, problemSize, Views, constraints, "Experiment4", "Computation", computation_time);
  write_datapoint(4, dimensionality, problemSize, Views, constraints, "Experiment4", "Conversion", conversion_time);
  write_datapoint(4, dimensionality, problemSize, Views, constraints, "Experiment4", "Cost Estimation", dec.setup_time);
  write_datapoint(4, dimensionality, problemSize, Views, constraints, "Experiment4", "ISL Space", dec.space_time);
  write_datapoint(4, dimensionality, problemSize, Views, constraints, "Experiment4", "ISL Func", dec.map_time);
  write_datapoint(4, dimensionality, problemSize, Views, constraints, "Experiment4", "ISL Solve", dec.solve_time);
  }
  

  
}

void experiment4_2(bool quiet=false) {
  idx_t dimensionality = 2;
  idx_t computation=3;
  idx_t views = 3;
  idx_t Views = views;
  idx_t problemSize = std::pow(2,30); //2^24
  double root = std::pow(problemSize, 1.0/3.0);
  size_t n = (size_t) root;
std::string constraints = "BestLayoutsAnalytical";
 using VIEW = View<double, Layout<2>>;
  VIEW A(new double[n*n], n,n);
  VIEW B(new double[n*n], n,n);
  VIEW C(new double[n*n], n,n);
  VIEW D(new double[n*n], n,n);
  VIEW E(new double[n*n], n,n);
  VIEW F(new double[n*n], n,n);
  VIEW G(new double[n*n], n,n);

  auto lam_init_e = [&](auto i, auto j) {
    E(i,j) = rand();
  };
  auto lam_init_f = [&](auto i, auto j) {
    F(i,j) = rand();
  };
  auto lam_init_g = [&](auto i, auto j) {
    G(i,j) = rand();
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


  auto init_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n));
  auto comp_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));

  auto init_e = make_kernel<INIT_POL>(init_seg, lam_init_e);
  auto init_f = make_kernel<INIT_POL>(init_seg, lam_init_f);
  auto init_g = make_kernel<INIT_POL>(init_seg, lam_init_g);
  
  auto comp1 = make_kernel<COMP_POL>(comp_seg, lam_comp1);
  auto comp2 = make_kernel<COMP_POL>(comp_seg, lam_comp2);
  auto comp3 = make_kernel<COMP_POL>(comp_seg, lam_comp3);

  std::array<idx_t, 2> sizes;
  sizes[0] = n; sizes[1] = n;//{{n,n}};
  auto blayout = make_permuted_layout(sizes, {{1,0}});
  auto dlayout = make_permuted_layout(sizes, {{1,0}});
  auto flayout01 = make_permuted_layout(sizes, {{0,1}});
  auto flayout10 = make_permuted_layout(sizes, {{1,0}});

  
  //BaseRAJA performance
  init_e();
  init_f();
  init_g();
  start();
  comp1();
  comp2();
  comp3();
  auto baseRAJA = stop();
  if(!quiet) {
    write_datapoint(4, dimensionality, problemSize, views, constraints, "BaseRAJA", "Computation", baseRAJA);
  }

  //HandTransformed
  permute_view(F, flayout01)();
  permute_view(B, flayout01)();
  permute_view(D, flayout01)();
  init_e();
  init_f();
  init_g();

  camp::idx_t run_time = 0;
  camp::idx_t conv_time = 0;
  
  start();
    permute_view(B, flayout10)();
  conv_time += stop();
  start();
    comp1();
  run_time += stop();
  start();
    permute_view(D, flayout10)();
    permute_view(F, flayout01)();
  conv_time += stop();
  start();
    comp2();
  run_time += stop();
  start();
    permute_view(F, flayout10)();
  conv_time += stop();
  start();
    comp3();
  run_time += stop();
  if(!quiet) {
  write_datapoint(4, dimensionality, problemSize, views, constraints, "HandTransformed", "Computation", run_time);
  write_datapoint(4, dimensionality, problemSize, views, constraints, "HandTransformed", "Conversion", conv_time);
  }
  
  //Experiment 
  permute_view(F, flayout01)();
  permute_view(B, flayout01)();
  permute_view(D, flayout01)();
  init_e();
  init_f();
  init_g();

  auto dec = format_decisions(tie(B,D,F), comp1, comp2, comp3);
  dec.set_format_for(B,flayout10, comp1);
  dec.set_format_for(D,flayout10, comp2);
  dec.set_format_for(F,flayout01, comp2);
  dec.set_format_for(F,flayout10, comp3);
  auto breakdown = dec.time_execution();
  auto conversion_time = get<0>(breakdown);
  auto computation_time = get<1>(breakdown);
  if(!quiet) {
  write_datapoint(4, dimensionality, problemSize, Views, constraints, "Experiment4", "Computation", computation_time);
  write_datapoint(4, dimensionality, problemSize, Views, constraints, "Experiment4", "Conversion", conversion_time);
  write_datapoint(4, dimensionality, problemSize, Views, constraints, "Experiment4", "Cost Estimation", dec.setup_time);
  write_datapoint(4, dimensionality, problemSize, Views, constraints, "Experiment4", "ISL Space", dec.space_time);
  write_datapoint(4, dimensionality, problemSize, Views, constraints, "Experiment4", "ISL Func", dec.map_time);
  write_datapoint(4, dimensionality, problemSize, Views, constraints, "Experiment4", "ISL Solve", dec.solve_time);
  }
  

  
}
*/


template <idx_t Views>
void view_count_experiment(bool quiet=false) {

  idx_t dimensionality = 2;
  idx_t computation=3;
  idx_t views = Views;
  std::string constraints="0";
  idx_t problemSize = 10000;

  double root = std::pow(problemSize, 1.0/dimensionality);
  camp::idx_t n = (camp::idx_t) root;

  //BaseRAJA definitions
  using VIEW = View<double, Layout<2>>;
  VIEW A(new double[n*n], n,n);
  VIEW B(new double[n*n], n,n);
  VIEW C(new double[n*n], n,n);
  VIEW D(new double[n*n], n,n);
  VIEW E(new double[n*n], n,n);
  VIEW F(new double[n*n], n,n);
  VIEW G(new double[n*n], n,n);

  auto lam_init_e = [&](auto i, auto j) {
    E(i,j) = rand();
  };
  auto lam_init_f = [&](auto i, auto j) {
    F(i,j) = rand();
  };
  auto lam_init_g = [&](auto i, auto j) {
    G(i,j) = rand();
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


  auto init_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n));
  auto comp_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));

  auto init_e = make_kernel<INIT_POL>(init_seg, lam_init_e);
  auto init_f = make_kernel<INIT_POL>(init_seg, lam_init_f);
  auto init_g = make_kernel<INIT_POL>(init_seg, lam_init_g);

  auto comp1 = make_kernel<COMP_POL>(comp_seg, lam_comp1);
  auto comp2 = make_kernel<COMP_POL>(comp_seg, lam_comp2);
  auto comp3 = make_kernel<COMP_POL>(comp_seg, lam_comp3);

 
  if constexpr (Views == 1) {
    auto dec = format_decisions(tie(B), comp1, comp2, comp3);

    auto breakdown = dec.time_execution();
    auto conversion_time = get<0>(breakdown);
    auto computation_time = get<1>(breakdown);
    if(!quiet) {
    write_datapoint(3, dimensionality, problemSize, Views, constraints, "View Count", "Computation", computation_time);
    write_datapoint(3, dimensionality, problemSize, Views, constraints, "View Count", "Conversion", conversion_time);
    write_datapoint(3, dimensionality, problemSize, Views, constraints, "View Count", "Cost Estimation", dec.setup_time);
    write_datapoint(3, dimensionality, problemSize, Views, constraints, "View Count", "ISL Space", dec.space_time);
    write_datapoint(3, dimensionality, problemSize, Views, constraints, "View Count", "ISL Func", dec.map_time);
    write_datapoint(3, dimensionality, problemSize, Views, constraints, "View Count", "ISL Solve", dec.solve_time);
    }
  } else if constexpr (Views == 2) {
    auto dec = format_decisions(tie(B,D), comp1, comp2, comp3);

    auto breakdown = dec.time_execution();
    auto conversion_time = get<0>(breakdown);
    auto computation_time = get<1>(breakdown);
    if(!quiet) {
    write_datapoint(3, dimensionality, problemSize, Views, constraints, "View Count", "Computation", computation_time);
    write_datapoint(3, dimensionality, problemSize, Views, constraints, "View Count", "Conversion", conversion_time);
    write_datapoint(3, dimensionality, problemSize, Views, constraints, "View Count", "Cost Estimation", dec.setup_time);
    write_datapoint(3, dimensionality, problemSize, Views, constraints, "View Count", "ISL Space", dec.space_time);
    write_datapoint(3, dimensionality, problemSize, Views, constraints, "View Count", "ISL Func", dec.map_time);
    write_datapoint(3, dimensionality, problemSize, Views, constraints, "View Count", "ISL Solve", dec.solve_time);
    }
  }
  else if constexpr (Views == 3) {
    auto dec = format_decisions(tie(B,D,F), comp1, comp2, comp3);

    auto breakdown = dec.time_execution();
    auto conversion_time = get<0>(breakdown);
    auto computation_time = get<1>(breakdown);
    if(!quiet) {
    write_datapoint(3, dimensionality, problemSize, Views, constraints, "View Count", "Computation", computation_time);
    write_datapoint(3, dimensionality, problemSize, Views, constraints, "View Count", "Conversion", conversion_time);
    write_datapoint(3, dimensionality, problemSize, Views, constraints, "View Count", "Cost Estimation", dec.setup_time);
    write_datapoint(3, dimensionality, problemSize, Views, constraints, "View Count", "ISL Space", dec.space_time);
    write_datapoint(3, dimensionality, problemSize, Views, constraints, "View Count", "ISL Func", dec.map_time);
    write_datapoint(3, dimensionality, problemSize, Views, constraints, "View Count", "ISL Solve", dec.solve_time);
    }
  } 

}

int main(int RAJA_UNUSED_ARG(argc), char** RAJA_UNUSED_ARG(argv[]))
{


  //problem sizes are 2^{15..20}
  std::cerr << "Warmup" << "\n";
  experiment1(100,true);
  experiment1(100,true);
  experiment1(100,true);

  std::cout << "Experiment,Dimensionality,Problem Size,Views,Constraints,Variant,Component,Time (microseconds)\n";

  std::cerr << "View Count Experiment\n";
  view_count_experiment<1>();
  view_count_experiment<2>();
  view_count_experiment<3>();
/*
  std::cerr << "Experiment2, warmup" << "\n";
  experiment2(3, true);
  experiment2(3, true);
  experiment2(3, true);
  for(int i = 0; i <= 4; i++) {
    std::cerr << "Experiment2, c=" << i << "\n";
    experiment2(i);
    experiment2(i);
    experiment2(i);
  }

  experiment3<1>(true);
  std::cerr << "Experiment3 v=1\n";
  experiment3<1>();
  experiment3<1>();
  experiment3<1>();
  std::cerr << "Experiment3 v=2\n";
  experiment3<2>(true);
  experiment3<2>();
  experiment3<2>();
  experiment3<2>();
  std::cerr << "Experiment3 v=3\n";
  experiment3<3>(true);
  experiment3<3>();
  experiment3<3>();
  experiment3<3>();

  experiment4_3();
  experiment4_3();
  experiment4_3();
  
  experiment4_2();
  experiment4_2();
  experiment4_2();
*/
}
