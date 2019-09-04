#ifndef RYSQ_CONSTANTS_H
#define RYSQ_CONSTANTS_H

namespace rysq {

  static const double SQRT_4_POW_PI_5 = 34.986836655249725;

  struct Zero {
    template<typename T>
    operator T() const { return T(0); }
    operator double() const { return double(0); }
  };

}

#endif /* RYSQ_CONSTANTS_H */
