#include "rysq/rysq.h"
#include "rysq/roots.h"
#include "rysq/vector.h"

#include <cassert>
#include <functional>
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
    double rAi, double rAB, double rBk,
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
    // Gn[a,i] = C[a]*Gn-1[a,i] + i*B0[a]*Gn-1[a,i-1] + (n-1)*B1[a]*Gn-2[a,i]

    double B2 = 1.0/(2*B);
    double B0[N], C[N], B1[N];
#define G(a,i,j) G[a+(i)*stride(N)+(j)*stride(N)*(m+1)]
    for (int a = 0; a < N; ++a) {
      C[a] = rBk + A*rAB*t2[a];
      B0[a] = 0.5*t2[a];
      B1[a] = (1.0 - A*t2[a])*B2;
      G(a,0,1) = C[a]*G(a,0,0);
    }

    {
      int j = 1;
      for (int i = 1; i <= m; ++i) {
        for (int a = 0; a < N; ++a) {
          G(a,i,j) = C[a]*G(a,i,j-1) + double(i)*B0[a]*G(a,i-1,j-1);
        }
      }
    }

    for (int j = 2; j <= n; ++j) {
      double n1 = j - 1;
      for (int a = 0; a < N; ++a) {
        G(a,0,j) = n1*B1[a]*G(a,0,j-2) + C[a]*G(a,0,j-1);
      }
      for (int i = 1; i <= m; ++i) {
        for (int a = 0; a < N; ++a) {
          G(a,i,j) = C[a]*G(a,i,j-1) + double(i)*B0[a]*G(a,i-1,j-1) + n1*B1[a]*G(a,i,j-2);
        }
      }
    }

#undef G

  }

  static inline void transfer(
    int N,
    int mi, int mj,
    double rij,
    double *G, double *I)
  {
#define G(a,i) G[a+(i)*N]
    for (int i = 0; i <= mi; ++i) {
      for (int a = 0; a < N; ++a) {
        I[a] = G(a,i);
      }
      I += N;
    }
    for (int j = 1; j <= mj; ++j) {
      for (int i = 0; i <= mi+mj-j; ++i) {
        for (int a = 0; a < N; ++a) {
          G(a,i) = rij*G(a,i) + G(a,i+1);
        }
      }
      for (int i = 0; i <= mi; ++i) {
        for (int a = 0; a < N; ++a) {
          I[a] = G(a,i);
        }
        I += N;
      }
    }
#undef G
#undef I
  }

  template<int N>
  static inline void transfer_bra(
    int mi, int mj,
    int nk, int nl,
    double rij,
    double *G, double *I)
  {
    for (int k = 0; k <= (nk+nl); ++k) {
      transfer(stride(N), mi, mj, rij, G, I);
      I += stride(N)*(mi+1)*(mj+1);
    }
  }

  template<int N, class Bra, class Ket>
  struct KernelImpl : Kernel2e<Bra,Ket> {

    Roots<N> roots_;
    double *buffer_ = NULL;

    using Kernel2e<Bra,Ket>::Kernel2e;
    using Kernel2e<Bra,Ket>::bra;
    using Kernel2e<Bra,Ket>::ket;

    ~KernelImpl() {
      delete[] buffer_;
    }

    struct Index {
      int stride[4];
      int operator()(int p, int q, int r, int s) const {
        return p*stride[0] + q*stride[1] + r*stride[2] + s*stride[3];
      }
    };

    struct Buffer {
      double *results;
      double *Ix, *Iy, *Iz;
      double *G;
    };

    auto buffer() {
      size_t nbf = shell::nbf(bra)*shell::nbf(ket);
      const auto &P = shell::get<0>(bra);
      const auto &Q = shell::get<1>(bra);
      const auto &R = shell::get<0>(ket);
      const auto &S = shell::get<1>(ket);
      size_t I = (stride(N)*(P.L+1)*(Q.L+1)*(R.L+1)*(S.L+1));
      size_t G = (stride(N)*(P.L+Q.L+1)*(R.L+S.L+1));
      if (!buffer_) {
        buffer_ = new double[nbf+G+I*3];
      }
      return Buffer{
        buffer_,
        buffer_ + nbf + I*0,
        buffer_ + nbf + I*1,
        buffer_ + nbf + I*2,
        buffer_ + nbf + I*3
      };
    }

    const double* compute(const shell::centers<Bra,Ket> &centers) override {

      const auto &P = shell::get<0>(bra);
      const auto &Q = shell::get<1>(bra);
      const auto &R = shell::get<0>(ket);
      const auto &S = shell::get<1>(ket);

      const auto &ri = centers.ri;
      const auto &rj = centers.rj;
      const auto &rk = centers.rk;
      const auto &rl = centers.rl;

      auto ai = P.a;
      auto aj = Q.a;
      auto ak = R.a;
      auto al = S.a;

      double A = ai + aj;
      double B = ak + al;

      double C = (P.C*Q.C*R.C*S.C);
      if (aj) C *= exp(-ai*aj*dot(ri-rj)/A);
      if (al) C *= exp(-ak*al*dot(rk-rl)/B);
      C /= (A*B*sqrt(A+B));

      Vector3 rA = center_of_charge(ai, ri, aj, rj);
      Vector3 rB = center_of_charge(ak, rk, al, rl);
      Vector3 rAB = rA - rB;

      double rho = (A*B)/(A + B);
      double X =  rho*dot(rAB);

      Vector<double,N> t2;
      Vector<double,N> W;

      if (!roots_.compute(X, t2, W)) {
        assert(false);
        return nullptr;
      }

      // if (N == 0) {
      //   double *ptr = this->buffer().results;
      //   ptr[0] = C*W[0];
      //   return ptr;
      // }

      t2 /= (A+B);

      auto buffer = this->buffer();
      double *I[3] = { buffer.Ix, buffer.Iy, buffer.Iz };

      for (int k = 0; k < 3; ++k) {
        double rAi = rA[k] - ri[k];
        double rBk = rB[k] - rk[k];
        double rAB = rA[k] - rB[k];
        const bool xsxs = !(Q.L || S.L);
        double *G = xsxs ? I[k] : buffer.G;
        recurrence<N>((P.L+Q.L), (R.L+S.L),  A, B, rAi, rAB, rBk, t2, W, G);
        if (xsxs) continue;
        transfer_bra<N>(P.L, Q.L, R.L, S.L, (ri[k]-rj[k]), G, I[k]);
        //transfer_ket<N>(P.L, Q.L, R.L, S.L, (rk[k]-rl[k]), buffer.transfer, I[k]);
        //transfer<N>(N*(P.L+1)*(Q.L+1), R.L, S.L, (rk[k]-rl[k]), G, I[k]);
      }

      Index index = {
        stride(N),
        stride(N)*(P.L+1),
        stride(N)*(P.L+1)*(Q.L+1),
        stride(N)*(P.L+1)*(Q.L+1)*(R.L+1)
      };

      double *ptr = buffer.results;
      for (auto s = S.begin(); s != S.end(); ++s) {
        for (auto r = R.begin(); r != R.end(); ++r) {
          for (auto q = Q.begin(); q != Q.end(); ++q) {
            for (auto p = P.begin(); p != P.end(); ++p) {
              auto *Ix = I[0] + index(p->x, q->x, r->x, s->x);
              auto *Iy = I[1] + index(p->y, q->y, r->y, s->y);
              auto *Iz = I[2] + index(p->z, q->z, r->z, s->z);
              double v = 0;
              for (int a = 0; a < N; ++a) {
                v += Ix[a]*Iy[a]*Iz[a]*W[a];
              }
              *ptr++ = SQRT_4PI5*C*v;
            }
          }
        }
      }

      return buffer.results;

    }

  };

  template<size_t N, class ... Args>
  std::unique_ptr< Kernel2e<Args...> > kernel(const Args& ... args) {
    return std::make_unique< KernelImpl<N,Args...> >(args...);
  }

  template<class Bra, class Ket, size_t ... Ns>
  std::unique_ptr< Kernel2e<Bra,Ket> >
  kernel(const Bra& bra, const Ket &ket, std::index_sequence<Ns...>) {
    typedef std::function<std::unique_ptr< Kernel2e<Bra,Ket> >(const Bra&, const Ket&)> F;
    static std::vector<F> kernel_table = { F(kernel<Ns,Bra,Ket>)... };
    size_t N = (L(bra) + L(ket))/2 + 1;
    if (N < kernel_table.size()) {
      return kernel_table.at(N)(bra,ket);
    }
    return nullptr;
  }


}

std::unique_ptr<rysq::Kernel2e2> rysq::kernel(const Bra<1> &bra, const Ket<1> &ket) {
  static const int MAX_ROOTS = (RYSQ_MAX_AM*2)/2 + 1;
  return kernel(bra, ket, std::make_index_sequence<MAX_ROOTS+1>{});
}

std::unique_ptr<rysq::Kernel2e3> rysq::kernel(const Bra<2> &bra, const Ket<1> &ket) {
  static const int MAX_ROOTS = (RYSQ_MAX_AM*3)/2 + 1;
  return kernel(bra, ket, std::make_index_sequence<MAX_ROOTS+1>{});
}

std::unique_ptr<rysq::Kernel2e4> rysq::kernel(const Bra<2> &bra, const Ket<2> &ket) {
  static const int MAX_ROOTS = (RYSQ_MAX_AM*4)/2 + 1;
  return kernel(bra, ket, std::make_index_sequence<MAX_ROOTS+1>{});
}
