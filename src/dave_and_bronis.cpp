#include "RAJA/RAJA.hpp"
#include <algorithm>
#include <chrono>

using namespace RAJA;

template <idx_t I, idx_t... Is>
struct order_to_kpol_helper {
  
  using SubPolicy = typename order_to_kpol_helper<Is...>::Policy;
  using Policy = statement::For<I, loop_exec, SubPolicy>;
};
template <idx_t I>
struct order_to_kpol_helper<I>{
  using Policy = statement::For<I, loop_exec, statement::Lambda<0>>;
};

template <idx_t I, idx_t...Is>
struct order_to_kpol {
  using SubPolicy = typename order_to_kpol_helper<Is...>::Policy;
  using Policy = KernelPolicy<
    statement::For<I,omp_parallel_for_exec,
      SubPolicy
    >
  >;
};

void warmup() {

  isl_ctx * ctx = isl_ctx_alloc();

  isl_ctx_free(ctx);
  idx_t n = 100;
  using VIEW = View<double, Layout<2>>;
  VIEW A(new double[n*n], n,n);
  using POL = typename order_to_kpol<0,1>::Policy;
  auto segs = make_tuple(RangeSegment(0,n), RangeSegment(0,n));
  auto knl = make_kernel<POL>(segs, [&](auto i, auto j) { A(i,j) = i + j;});
  auto dec = format_decisions(tie(A), knl);

  auto breakdown = dec.time_execution();
  camp::sink(breakdown);
  std::cerr << "warmup cost estimation: " << dec.setup_time << "\n";
  return;
};

auto write_header = []() {
  std::cout << "Experiment Name, Loop Nest Count,Loop Nest Depth,Memory Footprint,Data Dimensionality,View Count,Constraint Count,Component,Execution Time\n";
};
auto write_datapoint = [](auto experimentName, auto numLoopNests, auto depthLoopNests, auto memoryFootprint, auto dataDimensionality, auto viewCount, auto constraintCount, auto component, auto time) {
  std::cout << experimentName << "," << numLoopNests << "," << depthLoopNests << "," << memoryFootprint << "," << dataDimensionality << "," << viewCount << "," << constraintCount << "," << component << "," << time << "\n";
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

//first experiment to measure the impact of changing the number of loops in the sequence. 
//Constants:
//  Loop Nest Depth: 3
//  Memory Footprint: 3 * 100 * 100
//  Data Dimensionality: 2
//  View Count: 1 modelled, 3 in computation
//  Constraint Count: 0
template <idx_t... Loops>
void loop_count_experiment_impl(camp::idx_seq<Loops...>) {

  idx_t n = 100;
  auto numLoops = sizeof...(Loops);
  auto loopDepth = 3;
  auto footprint = 3 * n * n;
  auto dataDimensionality = 2;
  auto viewCount = 1;
  auto constraintCount = 0;
  using VIEW = View<double, Layout<2>>;
  VIEW A(new double[n*n], n,n);
  VIEW B(new double[n*n], n,n);
  VIEW C(new double[n*n], n,n);

   for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      A(i,j) = std::rand();
      B(i,j) = std::rand();
      C(i,j) = std::rand();
    }
  }
  auto lam1 = [&](auto i, auto j, auto k) {
    A(i,j) += B(i,k) * C(k,j);
  };
  auto lam2 = [&](auto i, auto j, auto k) {
    A(i,j) += C(i,k) * B(k,j);
  };
  auto lam3 = [&](auto i, auto j, auto k) {
    A(i,j) += B(i,k) * C(k,j);
  };
  auto lam4 = [&](auto i, auto j, auto k) {
    A(i,j) += C(i,k) * B(k,j);
  };
  auto lam5 = [&](auto i, auto j, auto k) {
    A(i,j) += B(i,k) * C(k,j);
  };
  auto lam6 = [&](auto i, auto j, auto k) {
    A(i,j) += C(i,k) * B(k,j);
  };
  auto lam7 = [&](auto i, auto j, auto k) {
    A(i,j) += B(i,k) * C(k,j);
  };
  auto lam8 = [&](auto i, auto j, auto k) {
    A(i,j) += C(i,k) * B(k,j);
  };

  using COMP_POL = typename order_to_kpol<0,1,2>::Policy;
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

  camp::sink(breakdown);

  
  write_datapoint("Loop Count",numLoops,loopDepth, footprint, dataDimensionality, viewCount, constraintCount, "Cost Estimation", dec.setup_time);
  write_datapoint("Loop Count",numLoops,loopDepth, footprint, dataDimensionality, viewCount, constraintCount, "ISL Setup", dec.space_time + dec.map_time);
  write_datapoint("Loop Count",numLoops,loopDepth, footprint, dataDimensionality, viewCount, constraintCount, "ISL Solve", dec.solve_time);

}
template<idx_t NumLoops>
void loop_count_experiment() {
  auto seq = idx_seq_from_to<0,NumLoops>();
  loop_count_experiment_impl(seq);
}

//second experiment to measure the impact of changing the number of loops in the sequence. 
//Constants:
//  Loop Count: 3
//  Memory Footprint: 3 * 100 * 100
//  Data Dimensionality: 2
//  View Count: 1 modelled, 3 in computation
//  Constraint Count: 0
template <idx_t ...Is>
void nest_depth_experiment_impl(camp::idx_seq<Is...>) {
  idx_t n = 100;
  auto numLoops = 3;
  idx_t constexpr loopDepth = sizeof...(Is);
  auto footprint = 3 * n * n;
  auto dataDimensionality = 2;
  auto viewCount = 1;
  auto constraintCount = 0;
  using VIEW = View<double, Layout<2>>;
  VIEW A(new double[n*n], n,n);
  VIEW B(new double[n*n], n,n);
  VIEW C(new double[n*n], n,n);

   for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      A(i,j) = std::rand();
      B(i,j) = std::rand();
      C(i,j) = std::rand();
    }
  }

  auto segs = tuple_repeat<loopDepth>(RangeSegment(0,n));

  auto lam1_2 = [&](auto i, auto j) { A(i,j) += B(j,i) * C(i,j); };
  auto lam2_2 = [&](auto i, auto j) { A(i,j) += C(j,i) * B(i,j); };
  auto lam3_2 = [&](auto i, auto j) { A(i,j) += B(j,i) * C(i,j); };

  auto lam1_3 = [&](auto i, auto j, auto k) { A(i,j) += B(i,k) * C(k,j); };
  auto lam2_3 = [&](auto i, auto j, auto k) { A(i,j) += C(i,k) * B(k,j); };
  auto lam3_3 = [&](auto i, auto j, auto k) { A(i,j) += B(i,k) * C(k,j); };

  auto lam1_4 = [&](auto i, auto j, auto k, auto l) { A(i,j) += B(i,l) * C(k,j); };
  auto lam2_4 = [&](auto i, auto j, auto k, auto l) { A(i,j) += C(i,l) * B(k,j); };
  auto lam3_4 = [&](auto i, auto j, auto k, auto l) { A(i,j) += B(i,l) * C(k,j); };
  
  camp::sink(lam1_2, lam2_2, lam3_2, lam1_3, lam2_3, lam3_3, lam1_4, lam2_4, lam3_4);
  using POL = typename order_to_kpol<Is...>::Policy;
  if constexpr (loopDepth == 2) {
    auto knl1 = make_kernel<POL>(segs, lam1_2);
    auto knl2 = make_kernel<POL>(segs, lam2_2);
    auto knl3 = make_kernel<POL>(segs, lam3_2);
    auto dec = format_decisions(tie(C), knl1, knl2, knl3);
    auto breakdown = dec.time_execution();
    camp::sink(breakdown);
    write_datapoint("Nest Depth",numLoops,loopDepth, footprint, dataDimensionality, viewCount, constraintCount, "Cost Estimation", dec.setup_time);
    write_datapoint("Nest Depth",numLoops,loopDepth, footprint, dataDimensionality, viewCount, constraintCount, "ISL Setup", dec.space_time + dec.map_time);
    write_datapoint("Nest Depth",numLoops,loopDepth, footprint, dataDimensionality, viewCount, constraintCount, "ISL Solve", dec.solve_time);

  }
  if constexpr (loopDepth == 3) {
    auto knl1 = make_kernel<POL>(segs, lam1_3);
    auto knl2 = make_kernel<POL>(segs, lam2_3);
    auto knl3 = make_kernel<POL>(segs, lam3_3);
    auto dec = format_decisions(tie(C), knl1, knl2, knl3);
    auto breakdown = dec.time_execution();
    camp::sink(breakdown);
    write_datapoint("Nest Depth",numLoops,loopDepth, footprint, dataDimensionality, viewCount, constraintCount, "Cost Estimation", dec.setup_time);
    write_datapoint("Nest Depth",numLoops,loopDepth, footprint, dataDimensionality, viewCount, constraintCount, "ISL Setup", dec.space_time + dec.map_time);
    write_datapoint("Nest Depth",numLoops,loopDepth, footprint, dataDimensionality, viewCount, constraintCount, "ISL Solve", dec.solve_time);

  }
  if constexpr (loopDepth == 4) {
    auto knl1 = make_kernel<POL>(segs, lam1_4);
    auto knl2 = make_kernel<POL>(segs, lam2_4);
    auto knl3 = make_kernel<POL>(segs, lam3_4);
    auto dec = format_decisions(tie(C), knl1, knl2, knl3);
    auto breakdown = dec.time_execution();
    camp::sink(breakdown);
    write_datapoint("Nest Depth",numLoops,loopDepth, footprint, dataDimensionality, viewCount, constraintCount, "Cost Estimation", dec.setup_time);
    write_datapoint("Nest Depth",numLoops,loopDepth, footprint, dataDimensionality, viewCount, constraintCount, "ISL Setup", dec.space_time + dec.map_time);
    write_datapoint("Nest Depth",numLoops,loopDepth, footprint, dataDimensionality, viewCount, constraintCount, "ISL Solve", dec.solve_time);

  }

}

template<idx_t Depth>
void nest_depth_experiment() {
  auto seq = idx_seq_from_to<0,Depth>();
  nest_depth_experiment_impl(seq);
}

//second experiment to measure the impact of changing the number of loops in the sequence. 
//Constants:
//  Loop Count: 3
//  Loop Nest Depth: 3
//  Data Dimensionality: 2
//  View Count: 1 modelled, 3 in computation
//  Constraint Count: 0
template <idx_t n>
void footprint_experiment() {
  auto numLoops = 3;
  auto loopDepth = 3;
  auto footprint = 3 * n * n;
  auto dataDimensionality = 2;
  auto viewCount = 1;
  auto constraintCount = 0;
  using VIEW = View<double, Layout<2>>;
  VIEW A(new double[n*n], n,n);
  VIEW B(new double[n*n], n,n);
  VIEW C(new double[n*n], n,n);

   for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      A(i,j) = std::rand();
      B(i,j) = std::rand();
      C(i,j) = std::rand();
    }
  }
  auto lam1 = [&](auto i, auto j, auto k) {
    A(i,j) += B(i,k) * C(k,j);
  };
  auto lam2 = [&](auto i, auto j, auto k) {
    A(i,j) += C(i,k) * B(k,j);
  };
  auto lam3 = [&](auto i, auto j, auto k) {
    A(i,j) += B(i,k) * C(k,j);
  };

  using COMP_POL = typename order_to_kpol<0,1,2>::Policy;
  auto segs = tuple_repeat<3>(RangeSegment(0,n));
  auto knl1 = make_kernel<COMP_POL>(segs, lam1);
  auto knl2 = make_kernel<COMP_POL>(segs, lam2);
  auto knl3 = make_kernel<COMP_POL>(segs, lam3);

  auto dec = format_decisions(tie(C), knl1, knl2, knl3);
  auto breakdown = dec.time_execution();
  camp::sink(breakdown);
  write_datapoint("Footprint",numLoops,loopDepth, footprint, dataDimensionality, viewCount, constraintCount, "Cost Estimation", dec.setup_time);
  write_datapoint("Footprint",numLoops,loopDepth, footprint, dataDimensionality, viewCount, constraintCount, "ISL Setup", dec.space_time + dec.map_time);
  write_datapoint("Footprint",numLoops,loopDepth, footprint, dataDimensionality, viewCount, constraintCount, "ISL Solve", dec.solve_time);

}


//creates the decision object for the different dimensionalities.
//the point of this function is to contain boilerplate so 
// the data_dim_experiment function is comprehensible
template <typename ViewType, idx_t... DataDimIdxs>
auto data_dims_experiment_decision_object(ViewType & A, ViewType & B, ViewType & C, idx_t n) {
  constexpr idx_t NumDataDims = sizeof...(DataDimIdxs);
  
  using POL = typename order_to_kpol<DataDimIdxs...>::Policy;
  auto segs = tuple_repeat<NumDataDims>(RangeSegment(0,n));

  if constexpr (NumDataDims == 2) {
    auto knl1 = make_kernel<POL>(segs, [&](auto i, auto j) {
      A(i,j) = B(i,j) * C(j,i);
    });
    auto knl2 = make_kernel<POL>(segs, [&](auto i, auto j) {
      A(j,i) = B(j,i) * C(j,i);
    });
    auto dec = format_decisions(tie(A,B,C), knl1, knl2);
    return dec; 
  } else if constexpr (NumDataDims == 3) {
    auto knl1 = make_kernel<POL>(segs, [&](auto i, auto j, auto k) {
      A(k,i,j) = B(i,k,j) * C(j,i,k);
    });
    auto knl2 = make_kernel<POL>(segs, [&](auto i, auto j, auto k) {
      A(j,i,k) = B(k,j,i) * C(j,i,k);
    });
    auto dec = format_decisions(tie(A,B,C), knl1, knl2);
    return dec; 
  } else if constexpr (NumDataDims == 4) {
    auto knl1 = make_kernel<POL>(segs, [&](auto i, auto j, auto k, auto l) {
      A(k,i,j,l) = B(l,i,k,j) * C(j,i,k,l);
    });
    auto knl2 = make_kernel<POL>(segs, [&](auto i, auto j, auto k, auto l) {
      A(j,l,i,k) = B(k,j,i,l) * C(j,i,l,k);
    });
    auto dec = format_decisions(tie(A,B,C), knl1, knl2);
    return dec; 
  } else if constexpr (NumDataDims == 5) {
    auto knl1 = make_kernel<POL>(segs, [&](auto i, auto j, auto k, auto l, auto m) {
      A(k,i,j,l,m) = B(m,l,i,k,j) * C(j,i,m,k,l);
    });
    auto knl2 = make_kernel<POL>(segs, [&](auto i, auto j, auto k, auto l, auto m) {
      A(j,m,l,i,k) = B(k,m,j,i,l) * C(j,i,l,k,m);
    });
    auto dec = format_decisions(tie(A,B,C), knl1, knl2);
    return dec; 
  } else if constexpr (NumDataDims == 6) {
    auto knl1 = make_kernel<POL>(segs, [&](auto i, auto j, auto k, auto l, auto m, auto n) {
      A(k,i,j,l,n,m) = B(m,l,n,i,k,j) * C(j,i,m,k,l,n);
    });
    auto knl2 = make_kernel<POL>(segs, [&](auto i, auto j, auto k, auto l, auto m, auto n) {
      A(j,m,l,i,k,n) = B(k,m,j,n,i,l) * C(n,j,i,l,k,m);
    });
    auto dec = format_decisions(tie(A,B,C), knl1, knl2);
    return dec; 
  } else {
  return 0.0;
  }
}


template <camp::idx_t...DataDimIdxs>
void data_dims_experiment(idx_t footprint, camp::idx_seq<DataDimIdxs...>) {
  constexpr idx_t NumDataDims = sizeof...(DataDimIdxs);
  auto numLoops = 2;
  auto loopDepth = NumDataDims;
  idx_t n = std::pow((footprint / 3.0),(double) 1.0 / NumDataDims);
  footprint = std::pow(n, NumDataDims);
  auto dataDimensionality = NumDataDims;
  auto viewCount = 1;
  auto constraintCount = 0;
  using VIEW = View<double, Layout<NumDataDims>>;
  
  auto nTuple = tuple_repeat<NumDataDims>(n);
  VIEW A(new double[footprint], camp::get<DataDimIdxs>(nTuple)...);
  VIEW B(new double[footprint], camp::get<DataDimIdxs>(nTuple)...);
  VIEW C(new double[footprint], camp::get<DataDimIdxs>(nTuple)...);
  

  auto dec = data_dims_experiment_decision_object<VIEW, DataDimIdxs...>(A,B,C,n);
  
  auto breakdown = dec.time_execution();
  camp::sink(breakdown);
  //auto conversion_time = get<0>(breakdown);
  //auto computation_time = get<1>(breakdown);
  write_datapoint("Data Dimensionality",numLoops,loopDepth, footprint, dataDimensionality, viewCount, constraintCount, "Cost Estimation", dec.setup_time);
  write_datapoint("Data Dimensionality",numLoops,loopDepth, footprint, dataDimensionality, viewCount, constraintCount, "ISL Setup", dec.space_time + dec.map_time);
  write_datapoint("Data Dimensionality",numLoops,loopDepth, footprint, dataDimensionality, viewCount, constraintCount, "ISL Solve", dec.solve_time);

}
template <camp::idx_t NumDataDims>
void data_dims_experiment(camp::idx_t footprint) {
  data_dims_experiment(footprint, idx_seq_from_to<0,NumDataDims>());
}



template <camp::idx_t constraintCount>
void constraint_count_experiment() {
  idx_t n = 100;
  auto numLoops = 5;
  idx_t constexpr loopDepth = 2;
  auto footprint = 3 * n * n;
  auto dataDimensionality = 2;
  auto viewCount = 1;
  using VIEW = View<double, Layout<2>>;
  VIEW A(new double[n*n], n,n);
  VIEW B(new double[n*n], n,n);
  VIEW C(new double[n*n], n,n);

   for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      A(i,j) = std::rand();
      B(i,j) = std::rand();
      C(i,j) = std::rand();
    }
  }

  auto segs = tuple_repeat<loopDepth>(RangeSegment(0,n));

  auto lam1 = [&](auto i, auto j) { A(i,j) += B(j,i) * C(i,j); };
  auto lam2 = [&](auto i, auto j) { A(i,j) += C(j,i) * B(i,j); };
  auto lam3 = [&](auto i, auto j) { A(i,j) += B(j,i) * C(i,j); };
  auto lam4 = [&](auto i, auto j) { C(i,j) += B(j,i) * C(i,j); };
  auto lam5 = [&](auto i, auto j) { A(i,j) += C(j,i) * B(i,j); };

  using COMP_POL = typename order_to_kpol<0,1>::Policy;
  auto knl1 = make_kernel<COMP_POL>(segs, lam1);
  auto knl2 = make_kernel<COMP_POL>(segs, lam2);
  auto knl3 = make_kernel<COMP_POL>(segs, lam3);
  auto knl4 = make_kernel<COMP_POL>(segs, lam4);
  auto knl5 = make_kernel<COMP_POL>(segs, lam5);
  auto dec = format_decisions(tie(C), knl1, knl2, knl3);

  if constexpr (constraintCount > 0) {
    dec.set_format_for(C,{0,1},knl1);
  }
  if constexpr (constraintCount > 1) {
    dec.set_format_for(C,{1,0},knl2);
  }
  if constexpr (constraintCount > 2) {
    dec.set_format_for(C,{0,1},knl3);
  }
  if constexpr (constraintCount > 3) {
    dec.set_format_for(C,{1,0},knl4);
  }
  if constexpr (constraintCount > 4) {
    dec.set_format_for(C,{0,1},knl5);
  }
  if constexpr (constraintCount > 5) {
    dec.set_output_format(C,{1,0});
  }

  auto breakdown = dec.time_execution();

  camp::sink(breakdown);

  write_datapoint("Constraint Count",numLoops,loopDepth, footprint, dataDimensionality, viewCount, constraintCount, "Cost Estimation", dec.setup_time);
  write_datapoint("Constraint Count",numLoops,loopDepth, footprint, dataDimensionality, viewCount, constraintCount, "ISL Setup", dec.space_time + dec.map_time);
  write_datapoint("Constraint Count",numLoops,loopDepth, footprint, dataDimensionality, viewCount, constraintCount, "ISL Solve", dec.solve_time);


}

int main() {

  warmup();
  write_header();
  
  std::cerr << "Running Loop Count Experiment...\n";
  std::cerr << "Loop Count = 1...\n";
  loop_count_experiment<1>();
  std::cerr << "Loop Count = 2...\n";
  loop_count_experiment<2>();
  std::cerr << "Loop Count = 3...\n";
  loop_count_experiment<3>();
  std::cerr << "Loop Count = 4...\n";
  loop_count_experiment<4>();
  std::cerr << "Loop Count = 5...\n";
  loop_count_experiment<5>();
  std::cerr << "Loop Count = 6...\n";
  loop_count_experiment<6>();
  std::cerr << "Loop Count = 7...\n";
  loop_count_experiment<7>();
  std::cerr << "Loop Count = 8...\n";
  loop_count_experiment<8>();
  std::cerr << "Loop Count Experiment Complete!\n";

  std::cerr << "Running Nest Depth Experiment...\n";
  std::cerr << "Nest Depth = 2...\n";
  nest_depth_experiment<2>();
  std::cerr << "Nest Depth = 3...\n";
  nest_depth_experiment<3>();
  std::cerr << "Nest Depth = 4...\n";
  nest_depth_experiment<4>();
  std::cerr << "Neth Depth Experiment Complete!\n";

  std::cerr << "Running Memory Footprint Experiment...\n";
  std::cerr << "View Side Length = 100...\n";
  footprint_experiment<100>();
  std::cerr << "View Side Length = 141...\n";
  footprint_experiment<141>();
  std::cerr << "View Side Length = 173...\n";
  footprint_experiment<173>();
  std::cerr << "View Side Length = 200...\n";
  footprint_experiment<200>();
  std::cerr << "View Side Length = 1000...\n";
  footprint_experiment<1000>();
  std::cerr << "Memory Footprint Experiment Complete!\n";

  std::cerr << "Running Constraint Count Experiment...\n";
  std::cerr << "Num Constraints = 0...\n";
  constraint_count_experiment<0>();
  std::cerr << "Num Constraints = 1...\n";
  constraint_count_experiment<1>();
  std::cerr << "Num Constraints = 2...\n";
  constraint_count_experiment<2>();
  std::cerr << "Num Constraints = 3...\n";
  constraint_count_experiment<3>();
  std::cerr << "Num Constraints = 4...\n";
  constraint_count_experiment<4>();
  std::cerr << "Num Constraints = 5...\n";
  constraint_count_experiment<5>();
  std::cerr << "Num Constraints = 6...\n";
  constraint_count_experiment<6>();
  std::cerr << "Constraint Count Experiment Complete!\n";

  std::cerr <<"Running Data Dimensionality Experiment...\n";
  auto footprint = 30000000;
  std::cerr <<"Num Data Dims = 2...\n";
  data_dims_experiment<2>(footprint);
  std::cerr <<"Num Data Dims = 3...\n";
  data_dims_experiment<3>(footprint);
  std::cerr <<"Num Data Dims = 4...\n";
  data_dims_experiment<4>(footprint);
  std::cerr <<"Num Data Dims = 5...\n";
  data_dims_experiment<5>(footprint);
  std::cerr <<"Num Data Dims = 6...\n";
  data_dims_experiment<6>(footprint);
  std::cerr << "Data Dimensionality Experiment Complete!\n";
}
