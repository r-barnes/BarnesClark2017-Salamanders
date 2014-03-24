#ifndef _salamander
#define _salamander

class Salamander {
  public:
    ///Gene Type - used for storing genetic information that is used to
    ///determine whether individuals are of the same species, based on a bitwise
    ///comparison of the binary expression of this number.
    typedef unsigned long long genetype;

    ///Initialize a new salamander. Is initialized as "alive", but with no
    ///parent (=-1), genome = 0, optimum temperature = 0, and
    ///probability of mutation = 1e-4.
    Salamander();

    ///Breed this salamander with another to make a baby! Returns a child
    ///salamander. Child's optimum temperature is the average of its parents,
    ///plus a mutation, drawn from a standard normal distribution. Child genome
    ///is based on a merge of the bit fields of the parents' genomes.
    Salamander breed(const Salamander &b) const;


    ///Determine whether two salamander genomes are similar. If the genomes are
    ///considered similary, then the individuals are members of the same
    ///species. species_sim_thresh is a number in [0,1] describing how similar
    ///two genomes must be before salamanders are deemed to be of the same
    ///species. Returns TRUE if salamanders are of the same species.
    bool pSimilar(const Salamander &b, double species_sim_thresh) const;
    
    ///Comparison of the similarity between two genomes. Used to compare genetic
    ///similarity and determine whether genes come from members of the same
    ///species. species_sim_thresh is a number in [0,1] describing how similar
    ///two genomes must be before salamanders are deemed to be of the same
    ///species. Returns TRUE if genes are of the same species.
    bool pSimilarGenome(const Salamander::genetype &b, double species_sim_thresh) const;

    ///Mutate this salamander's genome. Flips each element of the bit field with
    ///probability mutation_probability.
    void mutate();

    ///Determines whether a salamander dies given an input temperature and its
    ///optimum temperature. Based on a logit curve, parameterized such that a
    ///salamander dies with ~50% probability if it is more than 8 degrees C from
    ///its optimum temperature, and with ~90% probability if it is more than 12
    ///degrees from its optimum temperature. Returns TRUE if the salamander
    ///dies.
    bool pDie(double temp) const;

    ///Neutral genes. Determined by the parents of the salamander and used to
    ///determine if the salamander is of the same species as another salamander.
    genetype genes;

    ///Optimal temperature in degC for this salamander. For the first
    ///salamander, this is initialized as 33.5618604122814, which is the
    ///temperature (in degrees C) that we have projected for sea level in the
    ///Appalachians 65MYA. This value is determined for children by the average
    ///of their parents' otemp, and changes between generations based on a
    ///mutation rate described in mutate().
    double otempdegC;

    ///Mutation probability - used by mutate() to determine the probability of
    ///mutation in the genome. Note - this changes the genome that determines
    ///relatedness, speciation, and ability to breed. It does not directly alter
    ///optimum temperature. TODO: Is this mutations per million years?
    double mutation_probability;

    ///Drift rate for temperature optimum - used by breed() to determine the
    ///change in optimum temperature between children and parents. Note - this
    ///changes the temperature optimum, but does not directly influence
    ///relatedness, speciation, and ability to breed.
    ///TODO: Is this mutations per million years?
    double temperature_drift_sd;

    ///The phylogeny this salamander is part of (i.e., from what lineage did
    ///this salamander arise?) Is set to -1 for initialized salamanders, but
    ///should always be set to a positive number indicating the salamander's
    ///parent  species. The only exception to this is the first salamander
    ///("Eve"), which is its own parent, set at 0. TODO: This is more of a
    ///"speciesid"
    int parent;

    ///This variable adjusts how the the salamander's probability of death is
    ///affected by temperature. For further details, please look at the pDie()
    ///method
    double tempdeathfactor;    
};

#endif
