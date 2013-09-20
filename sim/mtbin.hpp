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

    ///Location of the first dead salamander on the list
    ///if startofdead==bin.size() then there are no more
    ///openings in the population
    int startofdead;

    ///Kill random salamanders with probability prob
    void mortaliate();

    ///Return the temperature of the bin at a given time
    double temp(double t) const;

    ///Return the height of the bin at a given time
    double height(double t) const;

    ///Kill all the salamanders in the bin and reset startofdead
    void killAll();

    ///Add a salamander to the bin. Fail silently if there's no room.
    void addSalamander(const Salamander &s);
};

#endif
