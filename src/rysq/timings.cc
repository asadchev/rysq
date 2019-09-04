#include "rysq/rysq.h"

#include <papi.h>
#include <iostream>
#include <assert.h>
#include <stdexcept>
#include <vector>
#include <chrono>

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


inline double timeit(size_t N, std::function<void(size_t)> f) {
  typedef std::chrono::high_resolution_clock clock;
  auto t = clock::now();
  for (size_t i = 0; i < N; ++i) {
    f(i);
  }
  //return std::chrono::duration_cast< std::chrono::duration<double> >(clock::now() - t).count();
  return std::chrono::duration<double>(clock::now() - t).count();
}

int main() {

  using namespace rysq;

  Shell a{ 1, 1.0, 1.0 };
  Shell b{ 1, 0.5, 2.0 };
  Shell c{ 1, 0.1, 3.0 };

  Vector3 r0 = { 0, 0, 0 };

  std::vector<double> buffer(RYSQ_MAX_CART*RYSQ_MAX_CART*RYSQ_MAX_CART);

  //PAPI::Init();

  double average = 0;

  for (int p = 0; p <= RYSQ_MAX_AM; ++p) {
    for (int q = 0; q <= RYSQ_MAX_AM; ++q) {
      for (int r = 0; r <= RYSQ_MAX_AM; ++r) {

        BraKet<Shell,Shell,Shell> braket = {
          { p, 1.0, 1.0 },
          { q, 1.0, 1.0 },
          { r, 1.0, 1.0 }
        };

        std::fill(buffer.begin(), buffer.end(), 0);

        auto eri = kernel(braket);

        std::cout << "Kernel " << braket.str() << ": "
                  << "~flops=" << eri->flops();

        size_t K = (1<<21)/eri->flops();

        double t = timeit(
          K,
          [&](size_t k){
            eri->compute(r0, r0, r0, buffer.data());
          }
        );

        average += t/K;

        std::cout << ", time=" << t << std::endl;

        // std::cout << buffer[0] << std::endl;

      }
    }
  }

  std::cout << "Kerenl av. time: " << average << std::endl;

}
