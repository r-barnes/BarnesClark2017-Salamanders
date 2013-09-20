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

    ///Number of bits shared between two genomes
    unsigned int similarity(const Salamander &b) const;

    ///Mutate this salamander's genome
    void mutate();

    ///Randomize this salamander's genome
    void randomizeGeneome();

    ///Neutral genes
    genetype genes;

    ///Optimal temperature for this salamander
    double otemp;

    bool dead;
};

#endif
