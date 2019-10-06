#include "rysq/rysq.h"
#include "rysq/roots/roots.h"
#include "rysq/roots/generated.h"

namespace rysq {
namespace cuda {

  struct Shell {
    int L;
    int K;
    RYSQ_GPU_ENABLED
    const auto& operator[](int idx) const {
      return orbitals_[idx];
    }
    struct Primitive {
      double a,C;
    };
    Primitive prims[20];
  private:
    rysq::Shell::Orbital orbitals_[RYSQ_MAX_CART];
  };


  RYSQ_GPU_ENABLED
  inline int nbf(const Shell &s) {
    return shell::nbf(s.L);
  }

  RYSQ_GPU_ENABLED
  inline int nprims(const Shell &s) {
    return s.K;
  }

  template<int N, int K, int P, class Orbitals>
  __device__
  void inner_kernel(const Orbitals &p, const double *Ix, const double *Iy, const double *Iz, double (&G)[P]) {

#define Ix(k,a,p) Ix[a + k*N + K*N*(p.x)]
#define Iy(k,a,p) Iy[a + k*N + K*N*(p.y)]
#define Iz(k,a,p) Iz[a + k*N + K*N*(p.z)]

    int k = threadIdx.x%K;

#pragma unroll 1
    for (int i = 0; i < P; ++i) {
      double gi = 0;
#pragma unroll 1
      for (int a = 0; a < N; ++a) {
        gi += Ix(a,k,p[i])*Iy(a,k,p[i])*Iz(a,k,p[i]);
      }
      G[i] += gi;
      //g[i] += Ix(1,k,p)*Iy(1,k,p)*Iz(1,k,p);
    }

  }

  template<int N, int PL>
  __global__
  void kernel(const Shell &P, const Shell &Q, double *G) {

    constexpr int NP = ((PL+1)*(PL+2))/2;
    constexpr int K = (32/NP);

    __shared__ double I[3][N*K*(PL+1)];
    __shared__ shell::Primitives primitives[32/NP];

    __shared__ double gp[NP];

    if (threadIdx.x < 30) {
      for (int q = threadIdx.y; q < nbf(Q); q += blockDim.y) {
        inner_kernel<N,K,NP>(P, I[0], I[1], I[2], gp);
      }
    }

    __syncthreads();

    for (int i = threadIdx.x; i < NP; i += blockDim.x) {
      G[i] = gp[i];
    }

  }

  // template
  // __global__
  // void kernel<2,4>(const Shell &P, const Shell &Q, double *G);

  template<int N>
  RYSQ_GPU_ENABLED
  double contract(const double *Ix, const double *Iy, const double *Iz) {
    double g = 0;
    #pragma unroll 1
    for (int i = 0; i < N; ++i) {
      g += Ix[i*3]*Iy[i*3]*Iz[i*3];
    }
    return g;
  }

  template<int N>
  RYSQ_GPU_ENABLED
  inline void recurrence(
    int m, int n,
    const double &A, const double &B,
    double rAB, double rAi, double rBk,
    double t2,
    double *G)
  {
    // bra recurrence G(a,m,0)
#define G(i) G[(i)*N]
    double C = (rAi - rAB*B*t2);
    G(0) = 1.0;
    G(1) = C;
    double A2 = 1.0/(2*A);
    double B1 = (1.0 - B*t2)*A2;
    #pragma unroll 1
    for (int i = 1; i < m; ++i) {
      G(i+1) = C*G(i) + double(i)*B1*G(i-1);
    }
#undef G
  }


  template<int N, int B>
  __global__
  void kernel_2c(const Shell &P, const Shell &R, double *G) {

    __shared__ double ri[3], rk[3];
    __shared__ double X[N], W[N];
    extern __shared__ double I[];

    const int thread = threadIdx.x+threadIdx.y*blockDim.x;
    double g[B] = { 0 };

#define I(i,j,x) (I + x + 3*N*(i + j*(P.L+1)))

    for (int k = 0; k < nprims(R); ++k) {
      for (int i = 0; i < nprims(P); ++i) {
        // compute I
        {
          if (thread < N) {
            double x = P.prims[i].a;
            roots3(x, X, W, thread);
          }
          __syncthreads();
          if (thread < 3*N) {
            int x = thread%3;
            int a = thread/3;
            double ai = P.prims[i].a;
            double ak = R.prims[k].a;
            double rA = ri[x];
            double rB = rk[x];
            double rAB = rA - rB;
            double rAi = 0;//rA - ri[x]; // 0?
            double rBk = 0;//rB - rk[x]; // 0?
            recurrence<N*3>(P.L, R.L, ai, ak, rAB, rAi, rBk, X[a], &I[x+a*3]);
          }
        }
        __syncthreads();
        double C = (P.prims[i].C*R.prims[k].C);
        auto p = P[threadIdx.x];
        //#pragma unroll
        for (int b = 0; b < B; ++b) {
          auto r = R[threadIdx.y+b*blockDim.y];
          double gb = C*contract<N>(I(p.x,r.x,0), I(p.y,r.y,1), I(p.z,r.z,2));
          //int idx = threadIdx.x + threadIdx.y*blockDim.x + b*(blockDim.x*blockDim.y);
          g[b] += gb;
        }
      }
    }

    __syncthreads();

    //#pragma unroll 1
    for (int b = 0; b < B; ++b) {
      int idx = threadIdx.x + threadIdx.y*blockDim.x + b*(blockDim.x*blockDim.y);
      G[idx] = g[b];
    }

  }

  template
  __global__
  void kernel_2c<3,3>(const Shell &P, const Shell &Q, double *G);


}
}
