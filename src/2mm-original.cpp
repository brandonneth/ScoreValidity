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


void original(camp::idx_t n) {
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

  start();
  comp1();
  comp2();
  auto elapsed = stop();

  std::cout << "Original,Total," << elapsed << "," << n << "\n";
}

int main(int argc, char ** argv) {
  
  camp::idx_t n = 1000;
  if (argc > 1) {
    n = n * std::atof(argv[1]);
  }
  original(n);
  original(n);
  original(n);
  return 0;
}
