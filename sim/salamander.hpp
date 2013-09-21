#ifndef _salamander
#define _salamander

class Salamander {
  public:
    ///Gene Type
    typedef unsigned long long genetype;

    ///Initialize a new salamander
    Salamander();

    ///Breed this salamander with another to make a baby
    Salamander breed(const Salamander &b) const;

    ///Number of bits shared between two salamanders
    bool pSimilar(const Salamander &b) const;
    
    ///Number of bits shared between two genomes
    bool pSimilarGenome(const Salamander::genetype &b) const;

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

    ///The phylogeny this is part of
    int parent;
};

#endif
