#include "phylo.hpp"
#include "salamander.hpp"
#include <cstdlib>
#include <queue>
#include <set>
#include <cassert>
#include <algorithm>

PhyloNode::PhyloNode(const Salamander &s, double t){
  genes=s.genes;
  emergence=t;
  parent=s.parent;
  lastchild=t;
  otemp=s.otemp;
}

void PhyloNode::addChild(int childNode){
  children.push_back(childNode);
}



Phylogeny::Phylogeny(const Salamander &s, double t){
  addNode(s,t);
}

void Phylogeny::addNode(const Salamander &s, double t){
  nodes.push_back(PhyloNode(s,t));
  //If this salamander species has a parent species, add this species as a child
  //species
  if(s.parent>=0)
    nodes[s.parent].addChild(nodes.size()-1);
}

void Phylogeny::UpdatePhylogeny(double t, std::vector<MtBin> &mts){
  for(auto &m: mts)     //Loop through parts of the mountain
  for(auto &s: m.bin){  //Loop through the salamanders on that part of the mountain
    //If I am dead or have no parent, skip me
    if(s.dead || s.parent==-1) continue;

    //If I am not similar to my parent
    if(!s.pSimilarGenome(nodes.at(s.parent).genes)) {
      bool trigger=false;

      //See if I am similar to any other salamanders, starting with the most recent
      //TODO: Don't search past my parent
      for(int p=nodes.size()-1;p>=0;--p) {
        //TODO: This needs to be adjusted based on the length of timesteps. This is crappy.
        //TODO: What is the below line even for?
        if( (t-nodes.at(s.parent).lastchild)>=1 ) continue;
        
        //If my parent is the same as this salamander's parent and my genes are
        //similar to this salamander's then this salamander and I are both part
        //of the first generation of a new species of salamander. Therefore, I
        //will consider myself this salamander's child, since its genome is
        //already stored in the phylogeny
        if(s.parent==nodes.at(p).parent && s.pSimilarGenome(nodes.at(p).genes) ) {
          s.parent=p;
          trigger=true;
          break;
        }
      }
      if(!trigger)     //No salamander in the phylogeny was similar to me!
        addNode(s,t);  //Therefore, I add myself to the phylogeny as a new species
    } else {  //I am similar to my parent, so mark my parent (species) as having survived this long
      nodes.at(s.parent).lastchild=t;
    }
  }
}

int Phylogeny::numAlive(double t) const {
  //One might think this function can be sped up by recognising that phylogenic
  //nodes are added in monotonically increasing order of time, but that forgets
  //that old nodes may survive all the way to the present. Thus, an exhaustive
  //search is necessary

  int sum=0;
  //Loop through the nodes of the phylogeny
  for(const auto &p: nodes){
    //If the phylogenic node came into being before the time in question and
    //persists past or up to the given time, then make a note that this species
    //is alive at the given time
    if(p.emergence<=t && t<=p.lastchild)
      ++sum;
  }
  return sum;
}

Phylogeny::mbdStruct Phylogeny::meanBranchDistance(double t) const {
  //Mean branch distance structure containing (avg branch dist, species) pairs
  mbdStruct mbd;

  //We start by finding those species which are alive at the given time
  std::vector<int> alive;
  for(unsigned int i=0;i<nodes.size();++i)
    if(nodes[i].emergence<=t && t<=nodes[i].lastchild)
      alive.push_back(i);
      
  //Enlarge to match size of alive. Initialize everything to 0.
  mbd.resize(alive.size(),std::pair<double,unsigned int>(0,0));

  //For each species that is alive
  for(unsigned int a=0;a<alive.size();++a){
    //Set the appropriate species ID
    mbd[a].second=alive[a];

    //Walk up the tree finding ancestors of A until we get to Eve
    std::set<int> parentsOfA;
    int p=alive[a];
    do {
      parentsOfA.insert(p);
    } while((p=nodes[p].parent)!=-1);

    //For each other node B, find its Lowest Common Ancestor with A
    for(unsigned int b=a+1;b<alive.size();++b){
      int p=alive[b];
      do {
        if(parentsOfA.count(p)) break;
      } while((p=nodes[p].parent)!=-1);
      assert(p!=-1); //No ancestor in common -> impossible!
      mbd[a].first+=t-nodes[p].emergence;
      mbd[b].first+=t-nodes[p].emergence;
    }
  }
  
  for(auto &i: mbd)
    i.first/=alive.size()-1;

  return mbd;
}



void Phylogeny::getECDF(const Phylogeny &p, double t) const {
  mbdStruct mbd=p.meanBranchDistance(t);
  //Build a cumulative distribution function from the existing phylogeny
  const unsigned int number_of_species=mbd.size();
  const unsigned int number_of_bins=100;
  std::vector<double> ecdf(number_of_bins,0);
  
  std::sort(mbd.begin(),mbd.end());
  
  double min_mbd      = mbd.front().first;
  double max_mbd      = mbd.back().first;
  double mbd_interval = (max_mbd-min_mbd)/(number_of_bins-1);

  int nbin = 0; // What ecdf bin are we in?
  for(unsigned int i=0;i<number_of_species; i++) {
      if(mbd[i].first<=min_mbd+i*mbd_interval) {
          ecdf[nbin]++;
      } else {
          ecdf[nbin]/=number_of_species; //standardize ecdf to 1;
          nbin++;
          ecdf[nbin]++;
      }
  }
}