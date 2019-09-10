#ifndef RYSQ_RYSQ_H
#define RYSQ_RYSQ_H

#include "rysq/config.h"
#include "rysq/constants.h"
#include "rysq/shell.h"
#include "rysq/vector.h"

#include <memory>
#include <string>

namespace rysq {

  typedef Vector<double,3> Vector3;

  using ::rysq::shell::Shell;
  using ::rysq::shell::Bra;
  using ::rysq::shell::Ket;

  template<class Bra, class Ket>
  struct Kernel2e {

    const Bra bra;
    const Ket ket;

    Kernel2e(const Bra &bra, const Ket &ket)
      : bra(bra), ket(ket)
    {
    }

    virtual ~Kernel2e() {}

    template<class ... Rs>
    const double* compute(const Rs& ... rs) {
      static_assert(
        sizeof...(Rs) == (Bra::size + Ket::size),
        "wrong number of shell centers"
      );
      return this->compute(shell::centers<Bra,Ket>{rs...});
    }

    virtual const double* compute(const shell::centers<Bra,Ket>&) = 0;

  };



  template<class Bra, class Ket>
  size_t flops(const Kernel2e<Bra,Ket> &kernel) {
    int N = (L(kernel.bra) + L(kernel.ket))/2 + 1;
    return N*3*nbf(kernel.bra)*nbf(kernel.ket);
  }

  template<class Bra, class Ket>
  std::string str(const Kernel2e<Bra,Ket> &kernel) {
    return "[" + shell::str(kernel.bra) + "|" + shell::str(kernel.ket) + "]";
  }

  typedef Kernel2e< Bra<1>, Ket<1> > Kernel2e2;
  typedef Kernel2e< Bra<2>, Ket<1> > Kernel2e3;
  typedef Kernel2e< Bra<2>, Ket<2> > Kernel2e4;

  std::unique_ptr<Kernel2e2> kernel(const Bra<1> &bra, const Ket<1> &ket);

  std::unique_ptr<Kernel2e3> kernel(const Bra<2> &bra, const Ket<1> &ket);

  std::unique_ptr<Kernel2e4> kernel(const Bra<2> &bra, const Ket<2> &ket);

}

#endif /* RYSQ_RYSQ_H */
