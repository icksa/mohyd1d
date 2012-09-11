#ifndef _0af5b546_db08_42c0_89fb_1e9b83f1b336
#define _0af5b546_db08_42c0_89fb_1e9b83f1b336

const double NAVO = 6.022e+23;
const double KBOL = 1.602e-12;
const double PI = 3.141592653589793;

namespace GeomType {
  enum E {
    CARTESIAN,
    SPHERICAL
  };
}

namespace HydroBc {
  enum E {
    FREE,
    FIXED
  };
}

#endif
