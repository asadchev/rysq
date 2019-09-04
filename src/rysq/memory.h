#include <unistd.h>

namespace rysq {

  inline size_t l1_cache_size() {
    return sysconf(_SC_LEVEL1_DCACHE_SIZE);
  }

  inline size_t l2_cache_size() {
    return sysconf(_SC_LEVEL2_CACHE_SIZE);
  }

}
