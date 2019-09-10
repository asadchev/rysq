#ifndef RYSQ_CONFIG_H
#define RYSQ_CONFIG_H

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

#endif /* RYSQ_CONFIG_H */
