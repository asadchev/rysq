#ifndef RYSQ_SHELL_H
#define RYSQ_SHELL_H

#include "rysq/constants.h"
#include "rysq/vector.h"

#include <stdint.h>
#include <type_traits>
#include <string>

namespace rysq {
namespace shell {

  inline int nbf(int L) {
    int n = L+1;
    return ((n*n+n)/2);
  }

  struct Shell {

    struct Unit;

    struct Orbital {
      //uint16_t x:4, y:4, z:4;
      uint16_t x, y, z;
    } __attribute__((aligned(16)));

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


  struct Shell::Unit {

    static const int L = 0;

    constexpr static const Zero a = {};
    constexpr static const double C = 1;

    struct Iterator {
      static const int x = 0, y = 0, z = 0;
      Iterator& operator++() {
        this->index_ = 1;
        return *this;
      }
      bool operator!=(const Iterator &it) const {
        return (it.index_ != this->index_);
      }
      const Iterator* operator->() const {
        return this;
      }
    private:
      friend class Unit;
      explicit Iterator(int index)
        : index_(index)
      {
      }
      int index_;
    };

    auto begin() const { return Iterator{0}; }
    auto end() const { return Iterator{1}; }

  };


  inline int nbf(const Shell &s) {
    return nbf(s.L);
  }

  inline int nbf(const Shell::Unit &s) {
    return 1;
  }

  template<class ... Args>
  int nbf(const Shell &s, Args&& ... args) {
    return nbf(s)*nbf(args...);
  }


  template<int>
  struct Tuple;

  template<>
  struct Tuple<1> {
    static const size_t size = 1;
    Shell first;
    const auto& get(std::integral_constant<int,0>) const { return first; }
    auto get(std::integral_constant<int,1>) const { return Shell::Unit{}; }
  };

  template<>
  struct Tuple<2> {
    static const size_t size = 2;
    Shell first;
    Shell second;
    const auto& get(std::integral_constant<int,0>) const { return first; }
    const auto& get(std::integral_constant<int,1>) const { return second; }
  };

  template<int K, int N>
  auto get(const Tuple<N> &t)
    -> decltype(t.get(std::integral_constant<int,K>{}))
  {
    return t.get(std::integral_constant<int,K>{});
  }

  template<int N>
  inline int L(const Tuple<N> &t) {
    return (get<0>(t).L + get<1>(t).L);
  }

  template<int N>
  inline int nbf(const Tuple<N> &t) {
    return nbf(get<0>(t))*nbf(get<1>(t));
  }

  inline std::string str(const Tuple<2> &p) {
    return std::to_string(p.first.L) + std::to_string(p.second.L);
  }

  inline std::string str(const Tuple<1> &p) {
    return std::to_string(p.first.L);
  }


  template<int K>
  using Bra = Tuple<K>;

  template<int K>
  using Ket = Tuple<K>;


  template<class Bra, class Ket>
  struct centers;

  template<>
  struct centers< Bra<2>, Ket<2> > {
    const Vector<double,3> &ri, &rj, &rk, &rl;
  };

  template<>
  struct centers< Bra<2>, Ket<1> > {
    const Vector<double,3> &ri, &rj, &rk;
    Vector<Zero,3> rl;
  };

  template<>
  struct centers< Bra<1>, Ket<1> > {
    const Vector<double,3> &ri, &rk;
    Vector<Zero,3> rj, rl;
  };

}
}

#endif /* RYSQ_SHELL_H */
