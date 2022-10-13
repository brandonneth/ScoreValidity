#include "RAJA/RAJA.hpp"
#include <algorithm>
#include <chrono>
#include <papi.h>
 long long counters[3];
  int PAPI_events[] = {
                PAPI_L1_DCM,
                PAPI_L2_DCM,
                PAPI_L2_DCA };

std::chrono::time_point<std::chrono::high_resolution_clock> start_;
void start() {
   start_ = std::chrono::high_resolution_clock::now();
  PAPI_start_counters(PAPI_events, 3);
}

auto stop() {
  auto end = std::chrono::high_resolution_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start_).count();
  PAPI_read_counters(counters, 3);

  return duration;
}

void original(camp::idx_t n) {
  using namespace RAJA;

  //BaseRAJA definitions
  using VIEW = View<double, Layout<2>>;
  std::array<camp::idx_t, 2> sizes{{n,n}};
  auto l10 = make_permuted_layout(sizes, {{1,0}});
  auto l01 = make_permuted_layout(sizes, {{0,1}}); 
  VIEW A(new double[n*n], n,n);
  VIEW B(new double[n*n], l10);
  VIEW C(new double[n*n], l10);
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
    statement::For<0, loop_exec, //omp_parallel_for_exec,
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

  std::cout << "Right,Total," << elapsed << "," << n << ",";
  std::cout << counters[0] << "," << counters[1] << "," << counters[2] << "," << counters[1] / counters[2] << "\n";

  delete[] A.get_data();
  delete[] B.get_data();
  delete[] C.get_data();
  delete[] D.get_data();
  delete[] T.get_data();
}

int main(int argc, char ** argv) {
  PAPI_library_init(PAPI_VER_CURRENT);
  camp::idx_t n = 1000;
  if (argc > 1) {
    n = n * std::atof(argv[1]);
  }
  original(n);
  original(n);
  original(n);
  return 0;
}
