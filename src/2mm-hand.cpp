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



void hand(camp::idx_t n) {
  using namespace RAJA;

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
  std::cout << "Hand,Computation," << run_time << "," << n << "\n";
  std::cout << "Hand,Conversion," << conv_time << "," << n << "\n";
  std::cout << "Hand,Total," << conv_time+run_time << "," << n << "\n";
}

int main(int argc, char ** argv) {
  
  camp::idx_t n = 1000;
  if(argc > 1) {
    n *= std::atoi(argv[1]);
  }
  hand(n);
  hand(n);
  hand(n);
  return 0;
}
