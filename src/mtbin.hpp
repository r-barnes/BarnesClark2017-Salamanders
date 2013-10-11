#ifndef _mtbin
#define _mtbin

#define binmax 1000

#include <array>
#include "salamander.hpp"

class MtBin {
  public:
    MtBin();
    MtBin(double binx);

    std::array<Salamander, binmax> bin;

    ///Height of this bin above sealevel
    double height;

    ///Location of the first dead salamander on the list
    ///if startofdead==bin.size() then there are no more
    ///openings in the population
    unsigned int startofdead;

    ///Kill random salamanders with probability prob
    void mortaliate(double t);

    ///Return the temperature of the bin at a given time
    double temp(double t) const;

    ///Return the area of the bin at a given time
    double area(double t) const;
    
    ///Karrying Kapacity of the bin given its area at a given time
    ///Returns a number [0, binsize].
    unsigned int kkap(double t) const;    

    ///Kill all the salamanders in the bin and reset startofdead
    void killAll();

    ///Add a salamander to the bin. Fail silently if there's no room.
    void addSalamander(const Salamander &s);

    ///Return number of living salamanders in this bin
    unsigned int alive() const;

    ///Breed salamanders within this bin if there is carrying capacity available
    void breed(double t);

    ///Swap some salamanders with those in another bin
    void diffuse(double t, MtBin &a);

  private:
    void killSalamander(int i);
};

#endif