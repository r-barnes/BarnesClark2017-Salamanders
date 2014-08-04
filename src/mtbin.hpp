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
    ///Maximum salamander populations per mountain bin
    double const binmax = 1000;

    ///If this is set to true, the mountains erode over time.
    bool vary_height;

    ///Alias for the type of container we are using to store the salamanders
    ///used in this bin
    typedef std::vector<Salamander> container;

    ///Define the bin used to store the salamanders
    container bin;

    ///Initializes this bin with elevation specified by heightkm0
    MtBin(double heightkm0);

    ///Returns the height of this bin IN KILOMETERS
    double height() const;

    ///Apply mortality to salamander within this bin based on how far they
    ///differ from optimal temperature and also on the the carrying capacity of
    ///the bin.
    void mortaliate(double tMyrs);

    ///Return the temperature of the bin at a given time, based on conditions at
    ///time tMyrs, in millions of year
    double temp(double tMyrs) const;

    ///Return the area of the bin at a elevationkm kilometers at time tMyrs, in
    ///millions of years.
    double area(double elevationkm, double tMyrs) const;
    
    ///Karrying Kapacity of the bin given its area at a given time tMyrs.
    ///Returns a number [0, binmax].
    unsigned int kkap(double tMyrs) const;    

    ///Add a salamander to the bin. Fail silently if there's no room.
    void addSalamander(const Salamander &s);

    ///Return number of living salamanders in this bin
    unsigned int alive() const;

    ///Breed salamanders within this bin if there is carrying capacity
    ///available. Carrying capacity depends on area that exists at the elevation
    ///band described by a particular mountain bin. Child is a new species if it
    ///differs from its similarity to its parents is less than
    ///species_sim_thresh, which takes values [0,1].
    void breed(double tMyrs, double species_sim_thresh);

    ///Salamanders have the opportunity to move up or down the mountain
    void diffuse(double t, MtBin *lower, MtBin *upper);

    container::iterator randomSalamander();

  private:
    ///Kills the indicated salamander by swapping it to the end of bin and then
    ///popping the back of bin. When used with an iterator, the iterator MUST
    ///decrement itself and then advance so that the swapped salamander is
    ///considered.
    void killSalamander(container::iterator s);

    ///Safely transfers salamander s from here to b
    void moveSalamanderTo(const MtBin::container::iterator &s, MtBin &b);

    ///Height of this bin above sealevel across all times IN KILOMETERS
    double heightkm;

    friend class MtBinUnitTest;
};


class MtBinUnitTest{
 public:
  void run() const;
};


#endif
