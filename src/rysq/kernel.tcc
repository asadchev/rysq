#include "rysq/rysq.h"
#include "rysq/roots/roots.h"
#include "rysq/vector.h"
#include "rysq/memory.h"

#include <vector>
#include <functional>
#include <cassert>

namespace rysq {

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
    const double *t2,
    double *G)
  {
    // bra recurrence G(a,m,0)
#define G(a,i) G[(a)+(i)*(N)]
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
#define G(a,i,j) G[a+(i)*(N)+(j)*(N)*(m+1)]
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
      transfer(N, mi, mj, rij, G, I);
      G += N*(mi+mj+1);
      I += N*(mi+1)*(mj+1);
    }
  }

  template<int N>
  static inline void transfer_ket(
    int mi, int mj,
    int nk, int nl,
    double r,
    double *G, double *I)
  {
    transfer(N*(mi+1)*(mj+1), nk, nl, r, G, I);
  }

  template<int N, class Bra, class Ket>
  struct KernelImpl : Kernel2e<Bra,Ket> {

    using Kernel2e<Bra,Ket>::bra;
    using Kernel2e<Bra,Ket>::ket;

    KernelImpl(const Bra &bra, const Ket &ket, const Parameters &params)
      : Kernel2e<Bra,Ket>(bra,ket,params),
        params_(params)
    {
      if (!this->params_.K) {
        const auto &P = shell::get<0>(bra);
        const auto &Q = shell::get<1>(bra);
        // fit working set in cache
        size_t K1 = 1*(N*(P.L+1)*3 + nbf(P));
        size_t K2 = 1*(N*(P.L+1)*(Q.L+1)*3 + nbf(bra));
        size_t K = std::min(
          l1_cache_size()/(sizeof(double)*K1),
          l2_cache_size()/(sizeof(double)*K2)
        );
        K = std::min<size_t>(nprims(bra)*nprims(ket), K);
        this->params_.K = std::max<size_t>(1,K);
      }
      assert(this->params_.K);
      //printf("K1=%i, K2=%i, K=%i\n", int(8*K*K1), int(8*K*K2), int(K));
    }

    Parameters params_;

    Roots<N> roots_;
    std::unique_ptr<double[]> buffer_;

    struct Index {
      int stride[4];
      int operator()(int p, int q, int r, int s) const {
        return p*stride[0] + q*stride[1] + r*stride[2] + s*stride[3];
      }
    };

    struct Buffer {
      double *results;
      double* I[3];
      double *Ix, *Gx, *Tx;
    };

    auto buffer() {
      size_t nbf = shell::nbf(bra)*shell::nbf(ket);
      const auto &P = shell::get<0>(bra);
      const auto &Q = shell::get<1>(bra);
      const auto &R = shell::get<0>(ket);
      const auto &S = shell::get<1>(ket);

      size_t K = std::max<size_t>(1, this->params_.K);
      size_t I = (N*(P.L+1)*(Q.L+1)*(R.L+1)*(S.L+1));
      size_t T = (N*(P.L+1)*(Q.L+1)*(R.L+S.L+1));
      size_t G = (N*(P.L+Q.L+1)*(R.L+S.L+1));

      if (!buffer_) {
        buffer_.reset(new double[nbf + K*I*3 + G+T+I]);
      }
      double *ptr = buffer_.get();

      Buffer buffer;
      buffer.results = ptr;
      ptr += nbf;

      for (size_t x = 0; x < 3; ++x) {
        buffer.I[x] = ptr;
        ptr += K*I;
      }

      double *Gx = ptr;
      if (!(Q.L || S.L)) {
        Gx = nullptr;
      }
      buffer.Gx = Gx;
      ptr += G;

      double *Tx = ptr;
      if (!Q.L) {
        Tx = Gx;
      }
      if (!S.L) {
        Tx = nullptr;
      }
      buffer.Tx = Tx;
      ptr += T;

      buffer.Ix = ptr;

      return buffer;

    }

    const double* compute(const shell::centers<Bra,Ket> &centers) override {

      const auto &ri = centers.ri;
      const auto &rj = centers.rj;
      const auto &rk = centers.rk;
      const auto &rl = centers.rl;

      const auto &bra_primitives = shell::primitives(bra, ri, rj);
      const auto &ket_primitives = shell::primitives(ket, rk, rl);

      std::vector<shell::Primitives2> primitives;
      primitives.reserve(this->params_.K);

      Buffer buffer = { nullptr };

      for (size_t kl = 0; kl < ket_primitives.size(); ++kl) {
        for (size_t ij = 0; ij < bra_primitives.size(); ++ij) {
          auto prims = shell::primitives(bra_primitives[ij], ket_primitives[kl]);
          if (SQRT_4_POW_PI_5*prims.C < this->params_.cutoff) continue;
          primitives.push_back(prims);
          bool last = ((kl+1 == ket_primitives.size()) && (ij+1 == bra_primitives.size()));
          if (!last && (primitives.size() < this->params_.K)) continue;
          if (!buffer.results) {
            buffer = this->buffer();
            std::fill(buffer.results, buffer.results+nbf(bra)*nbf(ket), double(0));
          }
          intermediates(primitives, centers, buffer);
          contract_intermediates(primitives.size(), buffer);
          primitives.clear();
        }
      }

      return buffer.results;

    }

  private:

    int intermediates(
      const std::vector<shell::Primitives2> &primitives,
      const shell::centers<Bra,Ket> &centers,
      Buffer &buffer)
    {

      const auto &P = shell::get<0>(bra);
      const auto &Q = shell::get<1>(bra);
      const auto &R = shell::get<0>(ket);
      const auto &S = shell::get<1>(ket);

      const auto &ri = centers.ri;
      const auto &rj = centers.rj;
      const auto &rk = centers.rk;
      const auto &rl = centers.rl;

      double *Gx = buffer.Gx;
      double *Tx = buffer.Tx;
      double *Ix = buffer.Ix;

      for (size_t k = 0; k < primitives.size(); ++k) {

        const auto &prims = primitives[k];

        const auto &C = prims.C;
        const auto &A = prims.A;
        const auto &B = prims.B;
        const auto &rA = prims.rA;
        const auto &rB = prims.rB;

        Vector3 rAB = rA - rB;

        Vector<double,N> t2;
        Vector<double,N> W;
        {
          double rho = (A*B)/(A + B);
          double X = rho*dot(rAB);
          if (!roots_.compute(X, t2, W)) {
            assert(false);
            return 0;
          }
          t2 /= (A+B);
        }

        for (int x = 0; x < 3; ++x) {
          double rAi = rA[x] - ri[x];
          double rBk = rB[x] - rk[x];
          double rAB = rA[x] - rB[x];
          assert(Gx || Ix);
          recurrence<N>((P.L+Q.L), (R.L+S.L),  A, B, rAi, rAB, rBk, t2, Gx ? Gx : Ix);
          if (Q.L) {
            assert(Tx || Ix);
            double r = ri[x]-rj[x];
            transfer_bra<N>(P.L, Q.L, R.L, S.L, r, Gx, Tx ? Tx : Ix);
          }
          if (S.L) {
            assert(Tx && Ix);
            double r = rk[x]-rl[x];
            transfer_ket<N>(P.L, Q.L, R.L, S.L, r, Tx, Ix);
          }
          size_t NK = N*primitives.size();
          assert(buffer.I[x]);
          for (int i = 0; i < (P.L+1)*(Q.L+1)*(R.L+1)*(S.L+1); ++i) {
            for (int a = 0; a < N; ++a) {
              buffer.I[x][a+k*N+i*NK] = Ix[a+i*N];
            }
            if (x != 0) continue;
            for (int a = 0; a < N; ++a) {
              double s = C*W[a];
              buffer.I[x][a+k*N+i*NK] *= s;
            }
          }
        }

      } // primitives

      return 1;

    }

    void contract_intermediates(int K, Buffer &buffer) {

      const auto &P = shell::get<0>(bra);
      const auto &Q = shell::get<1>(bra);
      const auto &R = shell::get<0>(ket);
      const auto &S = shell::get<1>(ket);

      auto I = buffer.I;

      Index index = {
        N*K,
        N*K*(P.L+1),
        N*K*(P.L+1)*(Q.L+1),
        N*K*(P.L+1)*(Q.L+1)*(R.L+1)
      };

      double * __restrict__ ptr = buffer.results;
      for (auto s = S.begin(); s != S.end(); ++s) {
        for (auto r = R.begin(); r != R.end(); ++r) {
          for (auto q = Q.begin(); q != Q.end(); ++q) {
            for (auto p = P.begin(); p != P.end(); ++p) {
              auto *Ix = I[0] + index(p->x, q->x, r->x, s->x);
              auto *Iy = I[1] + index(p->y, q->y, r->y, s->y);
              auto *Iz = I[2] + index(p->z, q->z, r->z, s->z);
              double v = inner_loop(K, Ix, Iy, Iz);
              *ptr++ += SQRT_4_POW_PI_5*v;
            }
          }
        }
      }

    }

    static inline double inner_loop(int K, const double *Ix, const double *Iy, const double *Iz) {
      int NK = N*K;
      double v = 0;
      for (int a = 0; a < NK; ++a) {
        v += Ix[a]*Iy[a]*Iz[a];
      }
      return v;
    }

  };

  template<size_t N, class Bra, class Ket>
  std::unique_ptr< Kernel2e<Bra,Ket> > kernel(const Bra &bra, const Ket& ket, const Parameters &params) {
    return std::make_unique< KernelImpl<N,Bra,Ket> >(bra, ket, params);
  }

  template<class Bra, class Ket, size_t ... Ns>
  std::unique_ptr< Kernel2e<Bra,Ket> >
  kernel(const Bra& bra, const Ket &ket, const Parameters &params, std::index_sequence<Ns...>) {
    typedef std::function<std::unique_ptr< Kernel2e<Bra,Ket> >(const Bra&, const Ket&, const Parameters &)> F;
    static std::vector<F> kernel_table = { F(kernel<Ns,Bra,Ket>)... };
    size_t N = (L(bra) + L(ket))/2 + 1;
    if (N < kernel_table.size()) {
      return kernel_table.at(N)(bra,ket,params);
    }
    return nullptr;
  }


}

std::unique_ptr<rysq::Kernel2e2> rysq::kernel(const Bra<1> &bra, const Ket<1> &ket, const Parameters &params) {
  static const int MAX_ROOTS = (RYSQ_MAX_AM*2)/2 + 1;
  return kernel(bra, ket, params,  std::make_index_sequence<MAX_ROOTS+1>{});
}

std::unique_ptr<rysq::Kernel2e3> rysq::kernel(const Bra<2> &bra, const Ket<1> &ket, const Parameters &params) {
  static const int MAX_ROOTS = (RYSQ_MAX_AM*3)/2 + 1;
  return kernel(bra, ket, params,  std::make_index_sequence<MAX_ROOTS+1>{});
}

std::unique_ptr<rysq::Kernel2e4> rysq::kernel(const Bra<2> &bra, const Ket<2> &ket, const Parameters &params) {
  static const int MAX_ROOTS = (RYSQ_MAX_AM*4)/2 + 1;
  return kernel(bra, ket, params,  std::make_index_sequence<MAX_ROOTS+1>{});
}
