//A simulation is comprised of a number of spatial bins representing elevation
//bands of the mountain range. Each bin contains many salamanders which interact
//within the bin and migrate between bins. This class specifies the spatial bin.
#ifndef _mtbin
#define _mtbin

#include <vector>
#include "salamander.hpp"
#include "params.hpp"

class MtBin {
 public:
	///Alias for the type of container we are using to store the salamanders
	///used in this bin
	typedef std::vector<Salamander> container;

	///Define the bin used to store the salamanders
	container bin;

	MtBin();

	///Initializes this bin with elevation specified by heightkm0
	MtBin(double heightkm);

	///Returns the height of this bin IN KILOMETERS
	double heightkm() const;

	///Returns the maximum height of the mountain range at the given time
	static double heightMaxKm(double tMyrs);

	///Apply mortality to salamander within this bin based on how far they
	///differ from optimal temperature and also on the the carrying capacity of
	///the bin.
	void mortaliate(double tMyrs, int max_species, int species_sim_thresh);

	///Return the temperature of the bin at a given time, based on conditions at
	///time tMyrs, in millions of year
	double temp(double tMyrs) const;

	///Return the area of the bin at a elevationkm kilometers at time tMyrs, in
	///millions of years.
	double area(double elevationkm, double tMyrs) const;

	//TODO: Cut?
	///Karrying Kapacity of the bin given its area at a given time tMyrs.
	///Returns a number [0, binmax].
	//unsigned int kkap(double tMyrs) const;

	///Add a salamander to the bin. Fail silently if there's no room.
	void addSalamander(const Salamander &s);

	///Return number of living salamanders in this bin
	unsigned int alive() const;

	///Breed salamanders within this bin if there is carrying capacity
	///available. Carrying capacity depends on area that exists at the elevation
	///band described by a particular mountain bin. Child is a new species if it
	///differs from its similarity to its parents is less than
	///species_sim_thresh, which takes values [0,1].
	void breed(double tMyrs, int species_sim_thresh);

	///Salamanders have the opportunity to move up or down the mountain if
	///advantageous
	void diffuseToBetter(double tMyrs, MtBin *lower, MtBin *upper);

	///Salamanders have the opportunity to move up or down the mountain
	void diffuseLocal(double tMyrs, MtBin *lower, MtBin *upper);

	///Salamanders have the opportunity to move all over the mountain
	void diffuseGlobal(double tMyrs, std::vector<MtBin> &mts);

	///Kills all of the salamanders in the bin
	void killAll();

	//Method for moving salamanders into a special separate bin representing the
	//surrounding lowlands.
	void diffuseToLowlands(MtBin &lowlands);

	//Method to be used by the surrounding lowlands to move salamanders back into
	//the active simulation.
	void diffuseFromLowlands(MtBin &frontrange);

	///Fetch an iterator to a random salamander from this bin
	container::iterator randomSalamander(int maxsal);

 private:
	///Kills the indicated salamander by swapping it to the end of bin and then
	///popping the back of bin. When used with an iterator, the iterator MUST
	///decrement itself and then advance so that the swapped salamander is
	///considered.
	void killSalamander(container::iterator s);

	///Safely transfers salamander s from here to b
	void moveSalamanderTo(const MtBin::container::iterator &s, MtBin &b);

	///Height of this bin above sealevel across all times IN KILOMETERS
	double heightkm_val;
};

#endif
