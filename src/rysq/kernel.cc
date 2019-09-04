#include "rysq/rysq.h"
#include "rysq/roots.h"
#include <vector>

namespace rysq {

  template<int Alignment = 1>
  constexpr int stride(int N) {
    return N+N%Alignment;
  }

  template<int N>
  static void recurrence(
    int m,
    const double *B0, const double *B1, const double *C,
    const double *Gm1, const double *Gm2,
    double *Gm)
  {
  }

  template<int N>
  static inline void recurrence(
    int m, int n,
    double A, double B,
    double rAi, double rAB, double rBK,
    const double *t2, const double *W,
    double *G)
  {
    // bra recurrence G(a,m,0)
#define G(a,i) G[(a)+(i)*stride(N)]
    for(int a = 0; a < N; ++a) {
      G(a,0) = 1.0;
    }
    if (m > 0) {
      for(int a = 0; a < N; ++a) {
        G(a,1) = (rAi - rAB*B*t2[a]);
      }
      double A2 = 1.0/(2*A);
      for (int i = 1; i < m; ++i) {
        for(int a = 0; a < N; ++a) {
          double B1 = (1.0 - B*t2[a])*A2;
          double C = G(a,1);
          G(a,i+1) = C*G(a,i) + double(i)*B1*G(a,i-1);
        }
      }
    }
#undef G

    if (n == 0) return;

    // ket recurrence G(a,m,n)
    double B0[N], B1[N], C[N];
#define G(a,i,j) G[a+(i)*stride(N)+(j)*stride(N)*(m+1)]
    for (int a = 0; a < N; ++a) {
      G(a,0,1) = C[a]*G(a,0,0);
    }
    for (int j = 2; j <= n; ++j) {
      // Gm[a,i] = i*B0[a]*Gm-1[a,i-1] + (m-1)*B1[a]*Gm-2[a,i] + C[a]*Gm-1[a,i]
      double m1 = j - 1;
      for (int a = 0; a < N; ++a) {
        G(a,0,j) = m1*B1[a]*G(a,0,j-2) + C[a]*G(a,0,j-1);
      }
      for (int i = 1; i <= m; ++i) {
        for (int a = 0; a < N; ++a) {
          G(a,i,j) = double(i)*B0[a]*G(a,i-1,j-1) + m1*B1[a]*G(a,i,j-2) + C[a]*G(a,i,j-1);
        }
      }
    }
#undef G

  }

  template<int N>
  static inline void transfer(
    int mi, int mj, int n,
    double *G, double *I)
  {
    double *Ik = I;
    for (int k = 0; k < n; ++k) {
      for (int j = 0; j <= mj; ++j) {
        double *Gk = nullptr;
        Gk = G  + k*stride(N)*(mi+mj+1);
        for (int i = 0; i <= mi; ++i) {
          for (int a = 0; a < N; ++a) {
            I[a] = Gk[a];
          }
          Ik += stride(N);
          Gk += stride(N);
        }
        Gk = G  + k*stride(N)*(mi+mj+1);
        for (int i = 0; i <= mi+mj-j; ++i) {
          for (int a = 0; a < N; ++a) {
            Gk[a] = Gk[a] + Gk[a+stride(N)];
          }
          Gk += stride(N);
        }
      }
    }
#undef G
#undef I
  }

  template<int N>
  struct KernelImpl : Kernel {

    Roots<N> roots_;
    double *buffer_ = NULL;

    using Kernel::Kernel;

    ~KernelImpl() {
      delete[] buffer_;
    }

    struct Index {
      int stride[3];
      int operator()(int p, int q, int r) const {
        return p*stride[0] + q*stride[1] + r*stride[2];
      }
    };

    struct Buffer {
      double *Ix, *Iy, *Iz;
      double *G;
    };

    auto buffer() {
      size_t I = (stride(N)*(P.L+1)*(Q.L+1)*(R.L+1));
      size_t G = (stride(N)*(P.L+Q.L+1)*(R.L+1));
      if (!buffer_) {
        buffer_ = new double[G+I*3];
      }
      return Buffer{
        buffer_ + I*0,
        buffer_ + I*1,
        buffer_ + I*2,
        buffer_ + I*3
      };
    }

    bool compute(const Vector3 &r0, const Vector3 &r1, const Vector3 &r2, double *V) override {

      double t2[N], W[N];

      double X = 1;
      roots_.compute(X, t2, W);

      auto buffer = this->buffer();

      double *I[3] = { buffer.Ix, buffer.Iy, buffer.Iz };

      for (int k = 0; k < 3; ++k) {
        double A = 1, B = 1, rAi = 1, rAB = 1, rBk = 1;
        if (Q.L == 0) {
          // no bra transfer
          auto *G = I[k];
          recurrence<N>((P.L), (R.L),  A, B, rAi, rAB, rBk, t2, W, G);
        }
        else {
          auto *G = buffer.G;
          recurrence<N>((P.L + Q.L), (R.L),  A, B, rAi, rAB, rBk, t2, W, G);
          transfer<N>(P.L, Q.L, R.L, G, I[k]);
        }
      }

      Index index = { stride(N), stride(N)*(P.L+1), stride(N)*(P.L+1)*(Q.L+1) };

      double *ptr = V;
      for (auto r = R.begin(); r != R.end(); ++r) {
        for (auto q = Q.begin(); q != Q.end(); ++q) {
          for (auto p = P.begin(); p != P.end(); ++p) {
            auto *Ix = I[0] + index(p->x, q->x, r->x);
            auto *Iy = I[1] + index(p->y, q->y, r->y);
            auto *Iz = I[2] + index(p->z, q->z, r->z);
            { // for (size_t k = 0; k < 36; ++k) {
              for (int a = 0; a < N; ++a) {
                *ptr += Ix[a]*Iy[a]*Iz[a];
              }
            }
            ++ptr;
          }
        }
      }

      return true;

    }

  };

  template<size_t N, class ... Args>
  std::unique_ptr<Kernel> kernel(const Args& ... args) {
    return std::make_unique< KernelImpl<N> >(args...);
  }

  template<class BraKet, size_t ... Ns>
  std::unique_ptr<Kernel> kernel(const BraKet& braket, std::index_sequence<Ns...>) {
    typedef std::function<std::unique_ptr<Kernel>(const BraKet&)> F;
    static std::vector<F> kernel_table = { F(kernel<Ns,BraKet>)... };
    size_t N = (braket.L())/2 + 1;
    if (N < kernel_table.size()) {
      return kernel_table.at(N)(braket);
    }
    return nullptr;
  }


}

std::unique_ptr<rysq::Kernel> rysq::kernel(const BraKet<Shell,Shell,Shell> &braket) {
  static const int MAX_ROOTS = (RYSQ_MAX_AM*3)/2 + 1;
  return kernel(braket, std::make_index_sequence<MAX_ROOTS+1>{});
}
