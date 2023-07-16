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

auto write_datapoint = [](auto experiment, auto dimensionality, auto problemSize, auto views, auto constraints, auto loops, auto component, auto time) {
  std::cout << experiment << "," << dimensionality << "," << problemSize << "," << views << "," << constraints << "," << loops << "," << component << "," << time << "\n";
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
// G  }
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
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "Computation", computation_time);
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "Conversion", conversion_time);
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "Cost Estimation", dec.setup_time);
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "ISL Space", dec.space_time);
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "ISL Func", dec.map_time);
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "ISL Solve", dec.solve_time);
    }
  } else if constexpr (Views == 2) {
    auto dec = format_decisions(tie(B,D), comp1, comp2, comp3);

    auto breakdown = dec.time_execution();
    auto conversion_time = get<0>(breakdown);
    auto computation_time = get<1>(breakdown);
    if(!quiet) {
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "Computation", computation_time);
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "Conversion", conversion_time);
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "Cost Estimation", dec.setup_time);
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "ISL Space", dec.space_time);
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "ISL Func", dec.map_time);
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "ISL Solve", dec.solve_time);
    }
  }
  else if constexpr (Views == 3) {
    auto dec = format_decisions(tie(B,D,F), comp1, comp2, comp3);

    auto breakdown = dec.time_execution();
    auto conversion_time = get<0>(breakdown);
    auto computation_time = get<1>(breakdown);
    if(!quiet) {
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "Computation", computation_time);
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "Conversion", conversion_time);
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "Cost Estimation", dec.setup_time);
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "ISL Space", dec.space_time);
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "ISL Func", dec.map_time);
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "ISL Solve", dec.solve_time);
    }
  } 
  else if constexpr (Views == 4) {
    auto dec = format_decisions(tie(B,D,F,C), comp1, comp2, comp3);

    auto breakdown = dec.time_execution();
    auto conversion_time = get<0>(breakdown);
    auto computation_time = get<1>(breakdown);
    if(!quiet) {
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "Computation", computation_time);
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "Conversion", conversion_time);
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "Cost Estimation", dec.setup_time);
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "ISL Space", dec.space_time);
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "ISL Func", dec.map_time);
    write_datapoint("View Count", dimensionality, problemSize, Views, constraints, computation, "ISL Solve", dec.solve_time);
    }
  } 

}

template <idx_t Constraints>
void constraint_count_experiment(bool quiet = false) {

  idx_t dimensionality = 2;
  idx_t computation=3;
  idx_t Views = 1;
  std::string constraints=std::to_string(Constraints);
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
  VIEW H(new double[n*n], n,n);

  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      A(i,j) = std::rand();
      B(i,j) = std::rand();
      C(i,j) = std::rand();
      D(i,j) = std::rand();
      E(i,j) = std::rand();
      F(i,j) = std::rand();
      G(i,j) = std::rand();
      H(i,j) = std::rand();
    }
  }
  auto lam1 = [&](auto i, auto j, auto k) {
    C(i,j) += A(i,k) * B(k,j);
  };
  auto lam2 = [&](auto i, auto j, auto k) {
    D(i,j) += E(i,k) * C(k,j);
  };
  auto lam3 = [&](auto i, auto j, auto k) { 
    F(i,j) += C(i,k) * G(k,j);
  };
  auto lam4 = [&](auto i, auto j, auto k) { 
    H(i,j) += E(i,k) * C(k,j);
  };

  auto segs = tuple_repeat<3>(RangeSegment(0,n));
  auto knl1 = make_kernel<COMP_POL>(segs, lam1);
  auto knl2 = make_kernel<COMP_POL>(segs, lam2);
  auto knl3 = make_kernel<COMP_POL>(segs, lam3);
  auto knl4 = make_kernel<COMP_POL>(segs, lam4);

  auto dec = format_decisions(tie(C), knl1, knl2, knl3, knl4);
  if constexpr (Constraints > 0) {
    dec.set_format_for(C, {0,1}, knl1);
  }
  if constexpr (Constraints > 1) {
    dec.set_format_for(C, {1,0}, knl2);
  }
  if constexpr (Constraints > 2) {
    dec.set_format_for(C, {0,1}, knl3);
  }
  if constexpr (Constraints > 3) {
    dec.set_format_for(C, {1,0}, knl4);
  } 
  auto breakdown = dec.time_execution();
  auto conversion_time = get<0>(breakdown);
  auto computation_time = get<1>(breakdown);
  if(!quiet) {
    write_datapoint("Constraint Count", dimensionality, problemSize, Views, constraints, computation, "Computation", computation_time);
    write_datapoint("Constraint Count", dimensionality, problemSize, Views, constraints, computation, "Conversion", conversion_time);
    write_datapoint("Constraint Count", dimensionality, problemSize, Views, constraints, computation, "Cost Estimation", dec.setup_time);
    write_datapoint("Constraint Count", dimensionality, problemSize, Views, constraints, computation, "ISL Space", dec.space_time);
    write_datapoint("Constraint Count", dimensionality, problemSize, Views, constraints, computation, "ISL Func", dec.map_time);
    write_datapoint("Constraint Count", dimensionality, problemSize, Views, constraints, computation, "ISL Solve", dec.solve_time);
  }
} // constraint_count_experiment

template <idx_t... Loops> 
void loop_count_experiment_impl(camp::idx_seq<Loops...>, bool quiet=false) {

   idx_t dimensionality = 2;
  idx_t computation=sizeof...(Loops);
  idx_t Views = 1;
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
  VIEW H(new double[n*n], n,n);

  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      A(i,j) = std::rand();
      B(i,j) = std::rand();
      C(i,j) = std::rand();
      D(i,j) = std::rand();
      E(i,j) = std::rand();
      F(i,j) = std::rand();
      G(i,j) = std::rand();
      H(i,j) = std::rand();
    }
  }
  auto lam1 = [&](auto i, auto j, auto k) {
    C(i,j) += A(i,k) * B(k,j);
  };
  auto lam2 = [&](auto i, auto j, auto k) {
    D(i,j) += E(i,k) * C(k,j);
  };
  auto lam3 = [&](auto i, auto j, auto k) { 
    F(i,j) += C(i,k) * G(k,j);
  };
  auto lam4 = [&](auto i, auto j, auto k) { 
    H(i,j) += E(i,k) * C(k,j);
  };
  auto lam5 = [&](auto i, auto j, auto k) {
    C(j,i) += A(i,k) * B(k,j);
  };
  auto lam6 = [&](auto i, auto j, auto k) {
    D(i,j) += E(i,k) * C(i,k);
  };
  auto lam7 = [&](auto i, auto j, auto k) { 
    F(i,j) += C(j,i) * G(k,j);
  };
  auto lam8 = [&](auto i, auto j, auto k) { 
    H(i,j) += E(i,k) * C(k,i);
  };

  auto segs = tuple_repeat<3>(RangeSegment(0,n));
  auto knl1 = make_kernel<COMP_POL>(segs, lam1);
  auto knl2 = make_kernel<COMP_POL>(segs, lam2);
  auto knl3 = make_kernel<COMP_POL>(segs, lam3);
  auto knl4 = make_kernel<COMP_POL>(segs, lam4);
  auto knl5 = make_kernel<COMP_POL>(segs, lam5);
  auto knl6 = make_kernel<COMP_POL>(segs, lam6);
  auto knl7 = make_kernel<COMP_POL>(segs, lam7);
  auto knl8 = make_kernel<COMP_POL>(segs, lam8);

  auto allKnls = make_tuple(knl1, knl2, knl3, knl4, knl5, knl6, knl7, knl8);

  auto dec = format_decisions(tie(C), camp::get<Loops>(allKnls)...);

  auto breakdown = dec.time_execution();
  auto conversion_time = get<0>(breakdown);
  auto computation_time = get<1>(breakdown);
  if(!quiet) {
    write_datapoint("Loop Count", dimensionality, problemSize, Views, constraints, computation, "Computation", computation_time);
    write_datapoint("Loop Count", dimensionality, problemSize, Views, constraints, computation, "Conversion", conversion_time);
    write_datapoint("Loop Count", dimensionality, problemSize, Views, constraints, computation, "Cost Estimation", dec.setup_time);
    write_datapoint("Loop Count", dimensionality, problemSize, Views, constraints, computation, "ISL Space", dec.space_time);
    write_datapoint("Loop Count", dimensionality, problemSize, Views, constraints, computation, "ISL Func", dec.map_time);
    write_datapoint("Loop Count", dimensionality, problemSize, Views, constraints, computation, "ISL Solve", dec.solve_time);
  } 

} // loop_count_experiment_impl
template <idx_t NumLoops>
void loop_count_experiment(bool quiet=false) {
  auto seq = idx_seq_from_to<0,NumLoops>();
  loop_count_experiment_impl(seq, quiet);
}

void dim_count_experiment_2(bool quiet=false) {
 idx_t dimensionality = 2;
  idx_t computation=3;
  idx_t views = 3;
  idx_t Views = views;
  idx_t problemSize = std::pow(65,2); //2^24
  idx_t constraints = 0;
  using VIEW = View<double, Layout<2>>;
  double root = std::pow(problemSize, 1.0/2.0);
  camp::idx_t  n = (camp::idx_t) root;
  VIEW A(new double[n*n], n,n);
  VIEW B(new double[n*n], n,n);
  VIEW C(new double[n*n], n,n);
  VIEW D(new double[n*n], n,n);
  VIEW E(new double[n*n], n,n);
  VIEW F(new double[n*n], n,n);
  VIEW G(new double[n*n], n,n);
  auto lam_comp1 = [&]( auto i, auto j, auto k) {
    E(i,j) += A(i,k) * B(k,j);
  };
  auto lam_comp2 = [&](auto i, auto j, auto k) {
    F(i,j) += C(i,k) * D(k,j);
  };
  auto lam_comp3 = [&]( auto i, auto j, auto k) {
    G(i,j) += E(i,k) * F(k,j);
  };

  using COMP_POL = KernelPolicy<
    statement::For<0, omp_parallel_for_exec,
      statement::For<1, loop_exec,
        statement::For<2, loop_exec,
            statement::Lambda<0>
        >
      >
    >
  >;
  auto init_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));
  auto comp_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));

  auto knl1 = make_kernel<COMP_POL>(comp_seg, lam_comp1);
  auto knl2 = make_kernel<COMP_POL>(comp_seg, lam_comp2);
  auto knl3 = make_kernel<COMP_POL>(comp_seg, lam_comp3);

  auto dec = format_decisions(tie(B,D,F), knl1, knl2, knl3);

  auto breakdown = dec.time_execution();
  auto conversion_time = get<0>(breakdown);
  auto computation_time = get<1>(breakdown);
  if(!quiet) {
    write_datapoint("Dimension Count", dimensionality, problemSize, Views, constraints, computation, "Computation", computation_time);
    write_datapoint("Dimension Count", dimensionality, problemSize, Views, constraints, computation, "Conversion", conversion_time);
    write_datapoint("Dimension Count", dimensionality, problemSize, Views, constraints, computation, "Cost Estimation", dec.setup_time);
    write_datapoint("Dimension Count", dimensionality, problemSize, Views, constraints, computation, "ISL Space", dec.space_time);
    write_datapoint("Dimension Count", dimensionality, problemSize, Views, constraints, computation, "ISL Func", dec.map_time);
    write_datapoint("Dimension Count", dimensionality, problemSize, Views, constraints, computation, "ISL Solve", dec.solve_time);
  } 

}



void dim_count_experiment_3(bool quiet=false) {
 idx_t dimensionality = 3;
  idx_t computation=3;
  idx_t views = 3;
  idx_t Views = views;
  idx_t problemSize = std::pow(65,3); //2^24
  idx_t constraints = 0;
  using VIEW = View<double, Layout<3>>;
  double root = std::pow(problemSize, 1.0/3.0);
  camp::idx_t  n = (camp::idx_t) root;
  VIEW A(new double[n*n*n], n,n,n);
  VIEW B(new double[n*n*n], n,n,n);
  VIEW C(new double[n*n*n], n,n,n);
  VIEW D(new double[n*n*n], n,n,n);
  VIEW E(new double[n*n*n], n,n,n);
  VIEW F(new double[n*n*n], n,n,n);
  VIEW G(new double[n*n*n], n,n,n);
  auto lam_comp1 = [&](auto l, auto i, auto j, auto k) {
    E(l,i,j) += A(l,i,k) * B(l,k,j);
  };
  auto lam_comp2 = [&](auto l, auto i, auto j, auto k) {
    F(l,i,j) += C(l,i,k) * D(l,k,j);
  };
  auto lam_comp3 = [&](auto l, auto i, auto j, auto k) {
    G(l,i,j) += E(l,i,k) * F(l,k,j);
  };

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
  n /= 2;
  auto init_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));
  auto comp_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));

  auto knl1 = make_kernel<COMP_POL>(comp_seg, lam_comp1);
  auto knl2 = make_kernel<COMP_POL>(comp_seg, lam_comp2);
  auto knl3 = make_kernel<COMP_POL>(comp_seg, lam_comp3);

  auto dec = format_decisions(tie(B,D,F), knl1, knl2, knl3);

  auto breakdown = dec.time_execution();
  auto conversion_time = get<0>(breakdown);
  auto computation_time = get<1>(breakdown);
  if(!quiet) {
    write_datapoint("Dimension Count", dimensionality, problemSize, Views, constraints, computation, "Computation", computation_time);
    write_datapoint("Dimension Count", dimensionality, problemSize, Views, constraints, computation, "Conversion", conversion_time);
    write_datapoint("Dimension Count", dimensionality, problemSize, Views, constraints, computation, "Cost Estimation", dec.setup_time);
    write_datapoint("Dimension Count", dimensionality, problemSize, Views, constraints, computation, "ISL Space", dec.space_time);
    write_datapoint("Dimension Count", dimensionality, problemSize, Views, constraints, computation, "ISL Func", dec.map_time);
    write_datapoint("Dimension Count", dimensionality, problemSize, Views, constraints, computation, "ISL Solve", dec.solve_time);
  } 
}

void dim_count_experiment_4(bool quiet=false) {
 idx_t dimensionality = 4;
  idx_t computation=3;
  idx_t views = 3;
  idx_t Views = views;
  idx_t problemSize = std::pow(65,4); //2^24
  idx_t constraints = 0;
  using VIEW = View<double, Layout<4>>;
  double root = std::pow(problemSize, 1.0/4.0);
  camp::idx_t  n = (camp::idx_t) root;
  VIEW A(new double[n*n*n*n], n,n,n,n);
  VIEW B(new double[n*n*n*n], n,n,n,n);
  VIEW C(new double[n*n*n*n], n,n,n,n);
  VIEW D(new double[n*n*n*n], n,n,n,n);
  VIEW E(new double[n*n*n*n], n,n,n,n);
  VIEW F(new double[n*n*n*n], n,n,n,n);
  VIEW G(new double[n*n*n*n], n,n,n,n);
  auto lam_comp1 = [&](auto m, auto l, auto i, auto j, auto k) {
    E(m,l,i,j) += A(m,l,i,k) * B(m,l,k,j);
  };
  auto lam_comp2 = [&](auto m, auto l, auto i, auto j, auto k) {
    F(l,i,m,j) += C(l,m,i,k) * D(l,k,j,m);
  };
  auto lam_comp3 = [&](auto m, auto l, auto i, auto j, auto k) {
    G(l,m,i,j) += E(m,l,i,k) * F(l,k,j,m);
  };

  using COMP_POL = KernelPolicy<
    statement::For<0, omp_parallel_for_exec,
      statement::For<1, loop_exec,
        statement::For<2, loop_exec,
          statement::For<3, loop_exec,
            statement::For<4, loop_exec,
              statement::Lambda<0>
            >
          >
        >
      >
    >
  >;
  n /= 4;
  auto comp_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));

  auto knl1 = make_kernel<COMP_POL>(comp_seg, lam_comp1);
  auto knl2 = make_kernel<COMP_POL>(comp_seg, lam_comp2);
  auto knl3 = make_kernel<COMP_POL>(comp_seg, lam_comp3);

  auto dec = format_decisions(tie(B,D,F), knl1, knl2, knl3);

  auto breakdown = dec.time_execution();
  auto conversion_time = get<0>(breakdown);
  auto computation_time = get<1>(breakdown);
  if(!quiet) {
    write_datapoint("Dimension Count", dimensionality, problemSize, Views, constraints, computation, "Computation", computation_time);
    write_datapoint("Dimension Count", dimensionality, problemSize, Views, constraints, computation, "Conversion", conversion_time);
    write_datapoint("Dimension Count", dimensionality, problemSize, Views, constraints, computation, "Cost Estimation", dec.setup_time);
    write_datapoint("Dimension Count", dimensionality, problemSize, Views, constraints, computation, "ISL Space", dec.space_time);
    write_datapoint("Dimension Count", dimensionality, problemSize, Views, constraints, computation, "ISL Func", dec.map_time);
    write_datapoint("Dimension Count", dimensionality, problemSize, Views, constraints, computation, "ISL Solve", dec.solve_time);
  } 
}


int main(int RAJA_UNUSED_ARG(argc), char** RAJA_UNUSED_ARG(argv[]))
{


  //problem sizes are 2^{15..20}
  std::cerr << "Warmup" << "\n";
  view_count_experiment<1>(true);

  std::cout << "Experiment,Dimensionality,Problem Size,Views,Constraints,Loops,Component,Time (microseconds)\n";

  std::cerr << "View Count Experiment\n";
  view_count_experiment<1>();
  view_count_experiment<2>();
  view_count_experiment<3>();
  view_count_experiment<4>();

  std::cerr << "Constraint Count Experiment\n";
  constraint_count_experiment<0>();
  constraint_count_experiment<1>();
  constraint_count_experiment<2>();
  constraint_count_experiment<3>();
  constraint_count_experiment<4>();

  std::cerr << "Loop Count Experiment \n";
  loop_count_experiment<2>();
  loop_count_experiment<3>();
  loop_count_experiment<4>();
  loop_count_experiment<5>();
  std::cerr << "lc 6\n";
  loop_count_experiment<6>();
  std::cerr << "lc 7\n";
  loop_count_experiment<7>();
  std::cerr << "lc 8\n";
  loop_count_experiment<8>();

  std::cout << "Dim Count Experiment\n";
  dim_count_experiment_2();
  dim_count_experiment_3();
  //dim_count_experiment_4();
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
