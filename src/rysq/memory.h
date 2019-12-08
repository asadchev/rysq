#include <unistd.h>
#include <cassert>

#if defined(_SC_LEVEL1_DCACHE_SIZE) && defined(_SC_LEVEL2_DCACHE_SIZE)
# define RYSQ_HAS_SC_LEVELX_DCACHE_SIZE
#else
# if __has_include(<sys/types.h>) && __has_include(<sys/sysctl.h>)
#  include <sys/types.h>
#  include <sys/sysctl.h>
#  define RYSQ_HAS_SYSCTL
# else
#  error "do not know how to query L{1,2} (data) cache size"
# endif
#endif

namespace rysq {

  inline size_t l1_cache_size() {
#ifdef RYSQ_HAS_SC_LEVELX_DCACHE_SIZE
    return sysconf(_SC_LEVEL1_DCACHE_SIZE);
#endif
#ifdef RYSQ_HAS_SYSCTL
    {
      int64_t cacheSizeFromSysctl = 0;
      size_t sz = sizeof(cacheSizeFromSysctl);
      const auto retcode = sysctlbyname("hw.l1dcachesize", &cacheSizeFromSysctl,
                                        &sz, nullptr, 0);
      assert(retcode == 0);
      return cacheSizeFromSysctl;
    }
#endif
  }

  inline size_t l2_cache_size() {
#ifdef RYSQ_HAS_SC_LEVELX_DCACHE_SIZE
    return sysconf(_SC_LEVEL2_DCACHE_SIZE);
#endif
#ifdef RYSQ_HAS_SYSCTL
    {
      int64_t cacheSizeFromSysctl = 0;
      size_t sz = sizeof(cacheSizeFromSysctl);
      const auto retcode = sysctlbyname("hw.l2cachesize", &cacheSizeFromSysctl,
                                        &sz, nullptr, 0);
      assert(retcode == 0);
      return cacheSizeFromSysctl;
    }
#endif
  }

  template<typename T>
  T* align(T *ptr, size_t alignment) {
    ptr += alignment-1;
    size_t r = size_t(ptr)%(sizeof(T)*alignment);
    ptr -= r/sizeof(T);
    assert(size_t(ptr)%(sizeof(T)*alignment) == 0);
    return ptr;
  }

}
