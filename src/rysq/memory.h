#include <unistd.h>
#include <cassert>

namespace rysq {

  inline size_t l1_cache_size() {
    return sysconf(_SC_LEVEL1_DCACHE_SIZE);
  }

  inline size_t l2_cache_size() {
    return sysconf(_SC_LEVEL2_CACHE_SIZE);
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
