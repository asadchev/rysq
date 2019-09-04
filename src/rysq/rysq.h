#ifndef RYSQ_RYSQ_H
#define RYSQ_RYSQ_H

#include <stdint.h>
#include <memory>

#ifndef RYSQ_MAX_AM
#define RYSQ_MAX_AM 5
#endif

#if (RYSQ_MAX_AM < 0)
#error "RYSQ_MAX_AM is invalid"
#endif

#define RYSQ_MAX_CART (((RYSQ_MAX_AM+2)*(RYSQ_MAX_AM+1))/2)

#ifdef __CUDACC__
#define RYSQ_GPU_ENABLED __host__ __device__
#else
#define RYSQ_GPU_ENABLED
#endif

namespace rysq {

  struct Vector3 {
    double x,y,z;
  };

  inline int nbf(int L) {
    int n = L+1;
    return ((n*n+n)/2);
  }

  struct Shell {

    struct Orbital {
      //uint16_t x:4, y:4, z:4;
      uint16_t x, y, z;
    }; //  __attribute__((aligned(16)));

    const int L;
    const double a, C;

    Shell(int L, double a, double C)
      : L(L), a(a), C(C)
    {
      auto it = this->orbitals_;
      for (int k = 0; k <= L; ++k) {
        for (int z = 0; z <= k; ++z) {
          it->x = L-k;
          it->y = k-z;
          it->z = z;
          ++it;
        }
      }
    }

    auto begin() const {
      return this->orbitals_;
    }

    auto end() const {
      return this->orbitals_ + nbf(this->L);
    }

  private:
    Orbital orbitals_[RYSQ_MAX_CART];

  };

  inline int nbf(const Shell &s) {
    return nbf(s.L);
  }

  template<class ... Args>
  int nbf(const Shell &s, Args&& ... args) {
    return nbf(s)*nbf(args...);
  }


  template<class P_, class Q_, class R_>
  struct BraKet {
    P_ P;
    Q_ Q;
    R_ R;
    int L() const {
      return P.L + Q.L + R.L;
    }
    std::string str() const {
      return (
        "(" + std::to_string(P.L) + std::to_string(Q.L) +
        "|" +
        std::to_string(R.L) + ")"
        );
    }
  };

  struct Kernel {

    const Shell P,Q,R;

    virtual bool compute(const Vector3&, const Vector3&, const Vector3&, double *buffer) = 0;

    Kernel(const BraKet<Shell,Shell,Shell> &braket)
      : P(braket.P), Q(braket.Q), R(braket.R)
    {
    }

    virtual ~Kernel() {}

    size_t flops() const {
      int N = (P.L + Q.L + R.L)/2 + 1;
      return N*3*nbf(P,Q,R);
    }

  };

  std::unique_ptr<Kernel> kernel(const BraKet<Shell,Shell,Shell> &braket);

}

#endif /* RYSQ_RYSQ_H */
