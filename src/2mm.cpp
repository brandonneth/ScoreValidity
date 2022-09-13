#include "RAJA/RAJA.hpp"
#include <algorithm>
#include <chrono>

std::chrono::time_point<std::chrono::high_resolution_clock> start_;
void start() {
   start_ = std::chrono::high_resolution_clock::now();
}

auto stop() {
  auto end = std::chrono::high_resolution_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start_).count();
  return duration;
}


void original() {
  using namespace RAJA;
  camp::idx_t n = (camp::idx_t) 1000;

  //BaseRAJA definitions
  using VIEW = View<double, Layout<2>>;
  VIEW A(new double[n*n], n,n);
  VIEW B(new double[n*n], n,n);
  VIEW C(new double[n*n], n,n);
  VIEW T(new double[n*n], n,n);
  VIEW D(new double[n*n], n,n);
  
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      T(i,j) = 0;
      D(i,j) = 0;
      A(i,j) = rand();
      B(i,j) = rand();
      C(i,j) = rand();
    }
  }

  using COMP_POL = KernelPolicy<
    statement::For<0, omp_parallel_for_exec,
      statement::For<1, loop_exec,
        statement::For<2, loop_exec,
            statement::Lambda<0>
        >
      >
    >
  >;

  auto lam1 = [&](auto i, auto j, auto k) {
    T(i,j) += A(i,k) * B(k,j);
  };

  auto lam2 = [&](auto i, auto j, auto k) {
    D(i,j) += T(i,k) * C(k,j);
  };
  auto comp_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));

  auto comp1 = make_kernel<COMP_POL>(comp_seg, lam1);
  auto comp2 = make_kernel<COMP_POL>(comp_seg, lam2);

  start();
  comp1();
  comp2();
  auto elapsed = stop();

  std::cout << "Original,Computation," << elapsed << "\n";
}

void hand() {
  using namespace RAJA;
  camp::idx_t n = (camp::idx_t) 1000;

  //BaseRAJA definitions
  using VIEW = View<double, Layout<2>>;
  VIEW A(new double[n*n], n,n);
  VIEW B(new double[n*n], n,n);
  VIEW C(new double[n*n], n,n);
  VIEW T(new double[n*n], n,n);
  VIEW D(new double[n*n], n,n);
  
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      T(i,j) = 0;
      D(i,j) = 0;
      A(i,j) = rand();
      B(i,j) = rand();
      C(i,j) = rand();
    }
  }

  using COMP_POL = KernelPolicy<
    statement::For<0, omp_parallel_for_exec,
      statement::For<1, loop_exec,
        statement::For<2, loop_exec,
            statement::Lambda<0>
        >
      >
    >
  >;

  auto lam1 = [&](auto i, auto j, auto k) {
    T(i,j) += A(i,k) * B(k,j);
  };

  auto lam2 = [&](auto i, auto j, auto k) {
    D(i,j) += T(i,k) * C(k,j);
  };
  auto comp_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));

  auto comp1 = make_kernel<COMP_POL>(comp_seg, lam1);
  auto comp2 = make_kernel<COMP_POL>(comp_seg, lam2);

  std::array<camp::idx_t, 2> sizes{{n,n}};
  auto l10 = make_permuted_layout(sizes, {{1,0}});
  auto l01 = make_permuted_layout(sizes, {{0,1}});

  camp::idx_t run_time = 0;
  camp::idx_t conv_time = 0;
  start();
    permute_view(B, l10)();
  conv_time += stop();
  start();
    comp1();
  run_time += stop();
  start();
    permute_view(C, l10)();
  conv_time += stop();
  start();
    comp2();
  run_time += stop();
  std::cout << "Hand,Computation," << run_time << "\n";
  std::cout << "Hand,Conversion," << conv_time << "\n";
  std::cout << "Hand,Total," << conv_time+run_time << "\n";
}

void locked() {
  using namespace RAJA;
  camp::idx_t n = (camp::idx_t) 1000;

  //BaseRAJA definitions
  using VIEW = View<double, Layout<2>>;
  VIEW A(new double[n*n], n,n);
  VIEW B(new double[n*n], n,n);
  VIEW C(new double[n*n], n,n);
  VIEW T(new double[n*n], n,n);
  VIEW D(new double[n*n], n,n);
  
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      T(i,j) = 0;
      D(i,j) = 0;
      A(i,j) = rand();
      B(i,j) = rand();
      C(i,j) = rand();
    }
  }

  using COMP_POL = KernelPolicy<
    statement::For<0, omp_parallel_for_exec,
      statement::For<1, loop_exec,
        statement::For<2, loop_exec,
            statement::Lambda<0>
        >
      >
    >
  >;

  auto lam1 = [&](auto i, auto j, auto k) {
    T(i,j) += A(i,k) * B(k,j);
  };

  auto lam2 = [&](auto i, auto j, auto k) {
    D(i,j) += T(i,k) * C(k,j);
  };
  auto comp_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));

  auto comp1 = make_kernel<COMP_POL>(comp_seg, lam1);
  auto comp2 = make_kernel<COMP_POL>(comp_seg, lam2);

  auto fd = format_decisions(tie(B,C), comp1, comp2);
  fd.set_format_for(B,{1,0},comp1);
  fd.set_format_for(B,{1,0},comp2);
  fd.set_format_for(C,{0,1},comp1);
  fd.set_format_for(C,{1,0},comp2);
  fd.set_output_format(B,{1,0});
  fd.set_output_format(C,{1,0});
  fd.lock();

  auto times = fd.time_execution();

  camp::idx_t run_time = camp::get<1>(times);
  camp::idx_t conv_time = camp::get<0>(times);
  
  std::cout << "Locked,Computation," << run_time << "\n";
  std::cout << "Locked,Conversion," << conv_time << "\n";
  std::cout << "Locked,Total," << conv_time+run_time << "\n";

  auto c = fd.finalize();

  start();
  c();
  auto t = stop();

  std::cout << "Locked,Finalize," << t << "\n";
}

int main() {

  original();
  original();
  original();
  hand();
  hand();
  hand();
  locked();
  locked();
  locked();
  return 0;
}