#include "phylo.hpp"
#include "salamander.hpp"
#include <cstdlib>

Phylo::Phylo(const Salamander &s, double t);
  genes=s.genes;
  emergence=t;
  parent=s.phylo_strain;
}
