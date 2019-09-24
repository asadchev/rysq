#include "rysq/rysq.h"

#if __has_include(<papi.h>)
# define RYSQ_HAS_PAPI
# include <papi.h>
#endif
#include <iostream>
#include <cassert>
#include <functional>
#include <stdexcept>
#include <vector>
#include <chrono>

#ifdef RYSQ_HAS_PAPI
struct PAPI {

  static void Init() {
    PAPI_library_init(PAPI_VER_CURRENT);
  }

  static auto Error(int r) {
    return std::runtime_error(
      std::string("PAPI: ") + PAPI_strerror(r)
    );
  }

  struct Flops {

    float real_time, proc_time, mflops;
    long long flpins;

    Flops() {
      int retval = PAPI_flops(&this->real_time, &this->proc_time, &this->flpins, &this->mflops);
      if (retval != PAPI_OK) {
        throw Error(retval);
      }
    }
  };

};
#endif  // RYSQ_HAS_PAPI

inline double timeit(size_t N, std::function<void(size_t)> f) {
  typedef std::chrono::high_resolution_clock clock;
  auto t = clock::now();
  for (size_t i = 0; i < N; ++i) {
    f(i);
  }
  //return std::chrono::duration_cast< std::chrono::duration<double> >(clock::now() - t).count();
  return std::chrono::duration<double>(clock::now() - t).count();
}

template<class Kernel, class ... Rs>
double time_kernel(Kernel &&kernel, const Rs& ... rs) {
  size_t K = (1<<21)/flops(*kernel);
  K = std::max<size_t>(K,1);
  double t = timeit(
    K,
    [&](size_t k){
      kernel->compute(rs...);
    }
  );
  return t/K;
}

int main() {

  using namespace rysq;

  Vector3 r0 = { 0, 0, 0 };

#ifdef RYSQ_HAS_PAPI
  PAPI::Init();
#endif

  double average = 0;
  for (int p = 0; p <= RYSQ_MAX_AM; ++p) {
    for (int q = 0; q <= RYSQ_MAX_AM; ++q) {
      Shell P = { p, 1.0, 1.0 };
      Shell R = { q, 1.0, 1.0 };
      average += time_kernel(rysq::kernel({P}, {R}), r0, r0);
    }
  }
  std::cout << "Kernel [X|X] av. time: " << average << std::endl;


  average = 0;
  for (int p = 0; p <= RYSQ_MAX_AM; ++p) {
    for (int q = 0; q <= RYSQ_MAX_AM; ++q) {
      for (int r = 0; r <= RYSQ_MAX_AM; ++r) {
        Shell P = { p, 1.0, 1.0 };
        Shell Q = { q, 1.0, 1.0 };
        Shell R = { r, 1.0, 1.0 };
        average += time_kernel(rysq::kernel({P,Q}, {R}), r0, r0, r0);
      }
    }
  }
  std::cout << "Kernel [XX|X] av. time: " << average << std::endl;

  average = 0;
  for (int p = 0; p <= RYSQ_MAX_AM; ++p) {
    for (int q = 0; q <= RYSQ_MAX_AM; ++q) {
      for (int r = 0; r <= RYSQ_MAX_AM; ++r) {
        for (int s = 0; s <= RYSQ_MAX_AM; ++s) {
          Shell P = { p, 1.0, 1.0 };
          Shell Q = { q, 1.0, 1.0 };
          Shell R = { r, 1.0, 1.0 };
          Shell S = { s, 1.0, 1.0 };
          average += time_kernel(rysq::kernel({P,Q}, {R,S}), r0, r0, r0, r0);
        }
      }
    }
  }
  std::cout << "Kernel [XX|XX] av. time: " << average << std::endl;

}
