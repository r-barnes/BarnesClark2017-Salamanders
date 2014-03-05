/**
Spatial bins in which salamanders reside. Each bin represents a particular
elevational band of the mountain.

*/
#ifndef _mtbin
#define _mtbin

#include <vector>
#include "salamander.hpp"

class MtBin {
  public:
    double const binmax = 1000; //Maximum salamander populations per mountain bin

    MtBin(double heightkm0);

    ///Returns the height of this bin IN KILOMETERS
    double height() const;

    typedef std::vector<Salamander> container;

    container bin;

    ///Salamanders fill the bin vector from 0 to startofdead-1. Bins beyond this
    ///are empty pre-allocations. startofdead therefore represents the location
    ///of the first dead salamander on the list.
    ///if startofdead==bin.size() then there are no more openings in the population
    unsigned int startofdead;

    ///Apply mortality to salamander within this bin based on how far they differ
    ///from optimal temperature and also on the the carrying capacity of the bin.
    void mortaliate(double tMyrs);

    ///Return the temperature of the bin at a given time, based on conditions at
    ///time tMyrs, in millions of year
    double temp(double tMyrs) const;

    ///Return the area of the bin at a elevationkm kilometers at time tMyrs, in millions of years.
    double area(double elevationkm, double tMyrs) const;
    
    ///Karrying Kapacity of the bin given its area at a given time tMyrs
    ///Returns a number [1, binmax]. //TODO: Is this range okay? Maybe not!
    unsigned int kkap(double tMyrs) const;    

    ///Kill all the salamanders in the bin and reset startofdead
    void killAll();

    ///Add a salamander to the bin. Fail silently if there's no room.
    ///Salamander is an object type that contains information about the
    ///genetics, lineage, &c. of a salamander.
    void addSalamander(const Salamander &s);

    ///Return number of living salamanders in this bin
    unsigned int alive() const;

    ///Breed salamanders within this bin if there is carrying capacity available.
    ///Carrying capacity depends on area that exists at the elevation band
    ///described by a particular mountain bin.
    ///Child is a new species if it differs from its similarity to its parents
    ///is less than species_sim_thresh, which takes values [0,1].
    void breed(double tMyrs, double species_sim_thresh);

    ///Swap some salamanders with those in another bin
    void diffuse(double tMyrs, MtBin &a);

    container::iterator randomSalamander();
    //container::const_iterator randomSalamander();

  private:
    ///Kills the indicated salamander and performs maintenance to keep living
    ///salamanders at the start of the list
    void killSalamander(container::iterator s);

    ///Height of this bin above sealevel across all times IN KILOMETERS
    double heightkm;
};

#endif
