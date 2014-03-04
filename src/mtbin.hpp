#ifndef _mtbin
//Spatial bins in which salamanders reside
//Each bin represents a particular elevational band
#define _mtbin

#define binmax 1000 //Maximum salamander populations per mountain bin

#include <array>
#include "salamander.hpp"

class MtBin {
  public:
    MtBin();
    MtBin(double binx);

    std::array<Salamander, binmax> bin;

    ///Height of this bin above sealevel
    double height;

    ///Salamanders fill the bin vector from 0 to startofdead. Bins beyond this
    ///are empty pre-allocations. startofdead therefore represents the location
    ///of the first dead salamander on the list.
    ///if startofdead==bin.size() then there are no more openings in the population
    unsigned int startofdead;

    ///Kill random salamanders with probability prob, based on conditions at time t
    void mortaliate(double t);

    ///Return the temperature of the bin at a given time, based on conditions at
    //time timeMyrs in millions of year
    double temp(double timeMyrs) const;

    ///Return the area of the bin at a given time timeMyrs, in millions of years.
    double area(double timeMyrs) const;
    
    ///Karrying Kapacity of the bin given its area at a given time t
    ///Returns a number [0, binsize].
    unsigned int kkap(double t) const;    

    ///Kill all the salamanders in the bin and reset startofdead
    void killAll();

    ///Add a salamander to the bin. Fail silently if there's no room.
    ///Salamander is an object type that contains information about the
    ///genetics, lineage, etc. of a salamander.
    void addSalamander(const Salamander &s);

    ///Return number of living salamanders in this bin
    unsigned int alive() const;

    ///Breed salamanders within this bin if there is carrying capacity available.
    ///Carrying capacity depends on area that exists at the elevation band
    ///described by a particular mountain bin.
    ///Child is a new species if it differs from its parents by more than
    ///species_sim_thresh.
    void breed(double t, double species_sim_thresh);

    ///Swap some salamanders with those in another bin
    void diffuse(double t, MtBin &a);

  private:
    void killSalamander(int i);
};

#endif
