#ifndef _mtbin
#define _mtbin

#include <array>
#include "salamander.hpp"

class MtBin {
  public:
    MtBin();
    MtBin(double binx);

    std::array<Salamander, 1000> bin;
    double binx;
    ///Kill random salamanders with probability prob
    void mortaliate();

    ///Return the temperature of the bin
    double temp() const;
};

#endif
