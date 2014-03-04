#ifndef _salamander
#define _salamander

class Salamander {
  public:
    ///Gene Type - used for storing genetic information that is used to determine
    ///whether individuals are of the same species, based on a bitwise comparison
    ///of the binary expression of this number.
    typedef unsigned long long genetype;

    ///Initialize a new salamander
    Salamander();

    ///Breed this salamander with another to make a baby
    Salamander breed(const Salamander &b) const;

    ///Print the salamanders genome as a bit field
    void printGenome() const;

    ///Number of bits shared between two salamanders
    bool pSimilar(const Salamander &b, double species_sim_thresh) const;
    
    ///Number of bits shared between two genomes
    bool pSimilarGenome(const Salamander::genetype &b, double species_sim_thresh) const;

    ///Mutate this salamander's genome
    void mutate();

    ///Randomize this salamander's genome
    void randomizeGeneome();

    ///Given an input temperature, should the salamander die?
    bool pDie(double temp) const;

    ///Neutral genes
    genetype genes;

    ///Optimal temperature for this salamander
    double otemp;

    bool dead;
    
    double mutation_probability;

    ///The phylogeny this is part of
    int parent;
};

#endif
