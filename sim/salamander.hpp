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
  private:
    genetype genes; ///< Neutral genes
    double otemp;      ///< Optimal temperature for this salamander
};

#endif
