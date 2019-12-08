#include "rysq/stieltjes.h"
#include <cstdlib>

namespace rysq {

  template<int N>
  struct stieltjes_grid_size {
    static const int value = 20+N*5;
  };

#define RYSQ_STIELTJES_GRID_SIZE(N,K)           \
  template<> struct stieltjes_grid_size<N> {    \
    static const int value = K;                 \
  };

  RYSQ_STIELTJES_GRID_SIZE(0,20);
  RYSQ_STIELTJES_GRID_SIZE(1,20);
  RYSQ_STIELTJES_GRID_SIZE(2,25);

  template<int N>
  struct Roots {
    Roots() : stieltjes_(Stieltjes<K>::instance())
    {
    }
    RYSQ_GPU_ENABLED
    bool compute(double x, double *X, double *W) const {
      return stieltjes_.template compute<N ? N : 1>(x, X, W);
    }
  private:
    static const int K = stieltjes_grid_size<N>::value;
    const Stieltjes<K> stieltjes_;
  };

}
