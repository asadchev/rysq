#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "rysq/rysq.h"
#include <vector>
#include <cmath>

using namespace rysq;

Vector3 r0 = { 0, 0, 0 };
Vector3 r1 = { 1, 1, 1 };
Vector3 r2 = { 2, 2, 2 };
Vector3 r3 = { 3, 3, 3 };

Shell S(0, 0.1, 1);
Shell P(1, 0.1, 1);
Shell D(2, 0.1, 1);
Shell F(3, 0.1, 1);

template<class Bra, class Ket>
void validate(
  rysq::Kernel2e<Bra,Ket> &kernel,
  const rysq::shell::centers<Bra,Ket> &rs,
  const std::initializer_list<double> &expected)
{
  auto *result = kernel.compute(rs);
  for (auto expect : expected) {
    double value = *result;
    if (fabs(value) < 1e-100) value = 0;
    BOOST_CHECK_CLOSE(value, expect, 0.5);
    ++result;
  }
}

template<class Bra, class Ket>
void validate(
  rysq::Kernel2e<Bra,Ket> &kernel,
  const rysq::shell::centers<Bra,Ket> &rs,
  double expected)
{
  auto *result = kernel.compute(rs);
  size_t n = nbf(kernel.bra)*nbf(kernel.ket);
  double value = 0;
  for (size_t i = 0; i < n; ++i) {
    value += result[i]/n;
  }
  BOOST_CHECK_CLOSE(value, expected, 0.5);
}

BOOST_AUTO_TEST_CASE(ssss) {
  validate(
    *rysq::kernel({S,S},{S,S}),
    {r0,r0,r0,r0},
    { 1382.98 }
  );
}

BOOST_AUTO_TEST_CASE(psss) {
  validate(
    *rysq::kernel({P,S},{S,S}),
    { {1,0,0}, r0, r0, r0 },
    { -760.321, 0, 0 }
  );
}

BOOST_AUTO_TEST_CASE(dsss) {
  validate(
    *rysq::kernel({D,S},{S,S}),
    { {1,0,0}, r0, r0, r0 },
    { 3171.94, 0, 0, 2721.63, 0, 2721.63}
  );
}

BOOST_AUTO_TEST_CASE(fsss) {
  validate(
    *rysq::kernel({F,S},{S,S}),
    { {1,0,0}, r0, r0, r0 },
    { -4921.31, 0, 0, -1550.04, 0, -1550.04, 0, 0, 0, 0 }
  );
}

BOOST_AUTO_TEST_CASE(psps) {
  validate(
    *rysq::kernel({P,S},{P,S}),
    { {1,0,0}, r0, {0,0,1}, r0 },
    { 440.339, 0, 15.0942, 0, 506.037, 0, 424.007, 0, 440.339 }
  );
}

BOOST_AUTO_TEST_CASE(dsps) {
  validate(
    *rysq::kernel({D,S},{P,S}),
    { {1,0,0}, r0, {0,0,1}, r0 },
    {
      -436.348, 0, 65.237, 177.548, 0, 180.235,
        0, -328.489, 0, 0, 75.4709, 0,
       -1718.14, 0, -285.407, -1463.06, 0, -1322.35
    }
  );
}

BOOST_AUTO_TEST_CASE(dsds) {
  validate(
    *rysq::kernel({D,S},{D,S}),
    { {1,0,0}, r0, {0,0,1}, r0 },
    {
      6965.72, 0, -88.0721, 5564.1, 0, 5588.93,
      0, 326.185, 0, 0, 13.434, 0,
      289.23, 0, 281.997, -113.075, 0, -88.0721,
      6384.41, 0, -113.075, 6294.51, 0, 5564.1,
      0, 215.414, 0, 0, 326.185, 0,
      7385.08, 0, 289.23, 6384.41, 0, 6965.72
    }
  );
}

BOOST_AUTO_TEST_CASE(ppss) {
  validate(
    *rysq::kernel({P,P},{S,S}),
    { {1,0,0}, {0,0,1}, r0, r0 },
    { 2278.41, 0, 221.592,
      0, 2571.02, 0,
      424.007, 0, 2278.41 }
  );
}

BOOST_AUTO_TEST_CASE(ddss) {
  validate(
    *rysq::kernel({D,D},{S,S}),
    { r0, r1, r2, r3 },
    2063.28
  );
}

BOOST_AUTO_TEST_CASE(ddds) {
  validate(
    *rysq::kernel({D,D},{D,S}),
    { r0, r1, r2, r3 },
    2652.25
  );
}

BOOST_AUTO_TEST_CASE(dddd) {
  validate(
    *rysq::kernel({D,D},{D,D}),
    { r0, r1, r2, r3 },
    5164.2
  );
}
