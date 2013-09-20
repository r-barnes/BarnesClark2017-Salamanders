#ifndef _salamander
#define _salamander

class Salamander {
  private:
    long double genes; ///< Neutral genes
    double otemp;      ///< Optimal temperature for this salamander
  public:
    ///Initialize a new salamander
    Salamander();

    ///Breed this salamander with another to make a baby
    Salamander breed(const Salamander &b) const;
};

#endif
