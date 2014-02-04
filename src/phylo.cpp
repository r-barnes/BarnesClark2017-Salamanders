#include "phylo.hpp"
#include "salamander.hpp"
#include <cstdlib>

PhyloNode::PhyloNode(const Salamander &s, double t){
  genes=s.genes;
  emergence=t;
  parent=s.parent;
  lastchild=t;
  otemp=s.otemp;
}




Phylogeny::Phylogeny(const Salamander &s, double t){
  nodes.push_back(PhyloNode(s,t));
}

void Phylogeny::UpdatePhylogeny(double t, std::vector<MtBin> &mts){
  for(auto &m: mts)
  for(auto &s: m.bin){
    //If I am dead or have no parent, skip me
    if(s.dead || s.parent==-1) continue;

    //If I am not similar to my parent
    if(!s.pSimilarGenome(nodes.at(s.parent).genes)) {
      bool trigger=false;

      //See if I am similar to any other salamanders, starting with the most recent
      for(int p=nodes.size()-1;p>=0;--p) {
        if( (t-nodes.at(s.parent).lastchild)>=1 ) continue; //TODO: This needs to be adjusted based on the length of timesteps. This is crappy.
        //If my parent is the same as this salamander's parent and my genes are similar to this salamander's
        //Then this salamander and I are both part of the first generation of a new species of salamander
        //Therefore, I will consider myself this salamander's child, since its genome is already stored in the phylogeny
        if(s.parent==nodes.at(p).parent && s.pSimilarGenome(nodes.at(p).genes) ) {
          s.parent=p;
          trigger=true;
          break;
        }
      }
      if(!trigger)                           //No salamander in the phylogeny was similar to me!
        nodes.push_back(PhyloNode(s, t));    //Therefore, I add myself to the phylogeny as a new species
    } else {  //I am similar to my parent, so mark my parent (species) as having survived this long
      nodes.at(s.parent).lastchild=t;
    }
  }
}

int Phylogeny::numAlive(double t) const {
  int sum=0;
  for(const auto &p: nodes){
    if(p.emergence<=t && t<=p.lastchild)
      ++sum;
  }
  return sum;
}