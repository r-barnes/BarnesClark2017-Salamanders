#include "phylo.hpp"
#include "salamander.hpp"
#include "mtbin.hpp"
#include <cstdlib>
#include <set>
#include <cassert>
#include <algorithm>
#include <fstream>
#include <string>
#include <iostream>
#include <cmath>

PhyloNode::PhyloNode(const Salamander &s, double t){
  //Copy relevant parameters from the Salamander that originates this strain
  genes     = s.genes;
  parent    = s.parent;
  otempdegC = s.otempdegC;
  //Set emergence and last known child to the current time
  emergence = t;
  lastchild = t;
}


//Add a link to a descendent species in the phylogeny
void PhyloNode::addChild(int n){
  children.push_back(n);
}


//Is this phylogenetic branch alive at the indicated time?
bool PhyloNode::aliveAt(double t) const {
  //If the phylogenetic node came into being before the time in question and
  //persists past or up to the given time, then this species is alive at the
  //given time
  return emergence<=t && t<=lastchild;
}


//Each time the phylogeny is updated, which happens every timestep, this is
//called for each salamander that is still similar enough to its parent to be
//part of the same species. The first time it is called, the lastchild is
//updated; thereafter, the statistics of the species for this particular
//timestep are updated.
void PhyloNode::updateWithSal(const MtBin &mt, const Salamander &s, double t){
  if(lastchild!=t || stats.size()==0){
    lastchild = t;
    stats.emplace_back(SpeciesStats(t));
  }
  stats.back().update(mt.heightkm(),s.otempdegC);
}



///////////////////
//Phylogeny Class
///////////////////


Phylogeny::Phylogeny(){}


//Create a new phylogeny with the indicated salamander as its progenitor. There
//is no preceding organism; therefore, calling this is something of a Biblical
//event (or some other story of Creation).
Phylogeny::Phylogeny(const Salamander &s, double t){
  addNode(s,t);
}


//Make the indicated salamander the progenitor of a new species which is a
//descendent species of the salamander's parent species
int Phylogeny::addNode(const Salamander &s, double t){
  nodes.push_back(PhyloNode(s,t));
  nodes[s.parent].addChild(nodes.size()-1);
  return nodes.size()-1; //Return the new node id
}


//This method loops through all existent salamanders and updates the
//phylogenetic tree to reflect which species have gone extinct, been born, or
//survived.
void Phylogeny::UpdatePhylogeny(double t, double dt, std::vector<MtBin> &mts){
  for(auto &m: mts)     //Loop through parts of the mountain
  for(auto &s: m.bin){  //Loop through the salamanders in this mountain bin
    //If I have no parent, skip me
    if(s.parent==-1)
      throw "Salamander with bad parent discovered!";

    //I am similar to my parent, so mark my parent (species) as having survived
    //this long
    if(s.pSimilarGenome(nodes.at(s.parent).genes, TheParams.speciesSimthresh())) {
      nodes.at(s.parent).updateWithSal(m,s,t);
      continue;
    }

    //If I am not similar to my parent then I may still be similar to one of my
    //parent's other children, which may already have an entry in the phylogeny.
    //Therefore, I'll search through recent entries in the phylogeny to see if
    //this was the case.
    bool has_parent=false;

    //See if I am similar to any other salamanders, starting with the most
    //recent. If I reach my parent, then I know I can't be related to any
    //species that was added to the phylogeny before my parent except by random
    //chance; therefore, I stop my search at that point.
    for(int p=nodes.size()-1;p>=s.parent;--p) {
      //If the last child of this potential parent was born more than one time
      //step ago, then this parent's lineage is dead and I cannot be a part of
      //it. This works because we are stepping by dt-Myr, so 2*dt-Myr is two
      //time steps.
      if( (t-nodes.at(p).lastchild)>=dt*2 ) continue;

      //If my parent is the same as this salamander's parent and my genes are
      //similar to this salamander's then this salamander and I are both part
      //of the first generation of a new species of salamander. Therefore, I
      //will consider myself this salamander's child, since its genome is
      //already stored in the phylogeny
      if( s.parent==nodes.at(p).parent && 
          s.pSimilarGenome(nodes.at(p).genes, TheParams.speciesSimthresh())
      ){
        s.parent   = p;
        has_parent = true;
        break;
      }
    }
    if(!has_parent)             //No salamander in the phylogeny was similar to me!
      s.parent = addNode(s,t);  //Therefore, I add myself to the phylogeny as a new species and make sure this salamander's parent is the new species
  }
}


//Determine the number of living species present in the phylogeny
int Phylogeny::livingSpecies(double t) const {
  //One might think this function can be sped up by recognising that
  //phylogenetic nodes are added in monotonically increasing order of time, but
  //that forgets that old nodes may survive all the way to the present. Thus, an
  //exhaustive search is necessary

  int sum=0;
  //Loop through the nodes of the phylogeny
  for(const auto &p: nodes)
    sum+=p.aliveAt(t);
  return sum;
}


Phylogeny::mbdStruct Phylogeny::meanBranchDistance(double t) const {
  //Mean branch distance structure containing (avg branch dist, species) pairs
  mbdStruct mbd;

  //We start by finding those species which are alive at the given time
  std::vector<int> alive;
  for(unsigned int i=0;i<nodes.size();++i)
    if(nodes[i].aliveAt(t))
      alive.push_back(i);

  //Enlarge to match size of alive. Initialise everything to 0.
  mbd.resize(alive.size(),std::pair<double,int>(0,0));

  //For each species that is alive
  for(unsigned int a=0;a<alive.size();++a){
    //Set the appropriate species ID
    mbd[a].second=alive[a];

    //Walk up the tree finding ancestors of A until we get to Eve. Eve is not
    //added to the list of ancestors, but that's okay, as we'll explain below.
    std::set<int> parentsOfA;
    int p=alive[a];
    do {
      parentsOfA.insert(p);
      p=nodes[p].parent;
    } while( p != nodes[p].parent);

    //For each other node B, find its Lowest Common Ancestor with A
    for(unsigned int b=a+1;b<alive.size();++b){
      int p=alive[b];
      //Walk up the list until we get to Eve. If we find a common ancestor
      //before Eve, then break early. Otherwise, Eve is the last ancestor we
      //would check and must be a common ancestor, therefore if we don't break
      //we terminate pointing at Eve, who is the LCA.
      do {
        if(parentsOfA.count(p)) break;
        p = nodes[p].parent;
      } while(p!=nodes[p].parent);
      mbd[a].first += t-nodes[p].emergence;
      mbd[b].first += t-nodes[p].emergence;
    }
  }

  for(auto &i: mbd)
    i.first/=alive.size()-1;

  return mbd;
}


//This method compares this phylogenetic tree to the one found by Kozak and
//Wiens using ECDF (empirical cumulative distribution function)
double Phylogeny::compareECDF(double t) const {
  mbdStruct mbd = meanBranchDistance(t);

  //Build a cumulative distribution function from the existing phylogeny
  const unsigned int number_of_species = mbd.size();
  const unsigned int number_of_bins    = 100;
  std::vector<double> ecdf(number_of_bins,0);

  std::sort(mbd.begin(),mbd.end());

  //The following values are taken from
  //data/Kozak_Plethodontid_Data/phylodist_cdf_ecdf.csv
  const double observed_ecdf[number_of_bins]={ 0.0208333333333333, 0.15625,
    0.208333333333333, 0.208333333333333, 0.239583333333333, 0.260416666666667,
    0.270833333333333, 0.270833333333333, 0.270833333333333, 0.270833333333333,
    0.270833333333333, 0.270833333333333, 0.270833333333333, 0.302083333333333,
    0.333333333333333, 0.427083333333333, 0.427083333333333, 0.427083333333333,
    0.489583333333333, 0.604166666666667, 0.625, 0.645833333333333, 0.65625,
    0.65625, 0.6875, 0.6875, 0.697916666666667, 0.697916666666667,
    0.697916666666667, 0.697916666666667, 0.697916666666667, 0.697916666666667,
    0.697916666666667, 0.697916666666667, 0.71875, 0.71875, 0.71875, 0.71875,
    0.71875, 0.71875, 0.71875, 0.71875, 0.71875, 0.71875, 0.71875, 0.71875,
    0.729166666666667, 0.729166666666667, 0.729166666666667, 0.739583333333333,
    0.739583333333333, 0.760416666666667, 0.760416666666667, 0.760416666666667,
    0.770833333333333, 0.770833333333333, 0.78125, 0.78125, 0.78125, 0.78125,
    0.78125, 0.78125, 0.78125, 0.78125, 0.78125, 0.78125, 0.78125, 0.78125,
    0.78125, 0.833333333333333, 0.875, 0.90625, 0.90625, 0.927083333333333,
    0.9375, 0.9375, 0.9375, 0.9375, 0.9375, 0.9375, 0.9375, 0.9375,
    0.979166666666667, 0.979166666666667, 0.979166666666667, 0.979166666666667,
    0.979166666666667, 0.979166666666667, 0.979166666666667, 0.979166666666667,
    0.979166666666667, 0.979166666666667, 0.979166666666667, 0.979166666666667,
    0.989583333333333, 0.989583333333333, 0.989583333333333, 0.989583333333333,
    0.989583333333333, 1};

  //The following values are taken from
  //data/Kozak_Plethodontid_Data/phylodist_cdf_bins.csv
  const double min_mbd      = 73.5209701979;
  const double max_mbd      = 118.6240912604;
  const double mbd_interval = (max_mbd-min_mbd)/(number_of_bins-1);

  int nbin = 0; //What ECDF bin are we in?
  for(unsigned int i=0;i<number_of_species; i++) {
    if(mbd[i].first<=min_mbd+i*mbd_interval) {
      ecdf[nbin]++;
    } else {
      ecdf[nbin]/=number_of_species; //standardize ECDf to 1;
      assert(ecdf[nbin]<=1);
      assert(ecdf[nbin]>=0);
      nbin++;
      ecdf[nbin]++;
    }
  }

  //This divides the last bin, which is not otherwise caught by the above loop
  if(number_of_species>0)
    ecdf[nbin]/=number_of_species;

  //Compare simulated ecdf with actual ecdf
  double sum_squared_difference=0;
  for(unsigned int i=0; i<number_of_bins;i++) {
    double squared_diff=std::pow(ecdf[i]-observed_ecdf[i],2);
    //std::cerr<<ecdf[i]<<" "<<observed_ecdf[i]<<" "<<squared_diff<<std::endl;
    assert(0<=squared_diff && squared_diff<=1);
    assert(i<100);
    sum_squared_difference+=squared_diff;
  }
  //std::cerr<<"SS diff: "<<sum_squared_difference<<std::endl;

  return sum_squared_difference;
}


//Print out a CSV of when each species emerged and died, and the optimal
//temperature of that species when it emerged.
void Phylogeny::persistGraph(int run_num, std::ofstream &out) const {
  //NOTE: Printing the species ID twice is done so that the species can be
  //plotted as a line with 'i', the species ID, on the vertical axis and the
  //species emergence and lastchild on the horizontal axis.
  if(run_num==0)
    out<<"RunNum, Emergence, Species, Last Child, Species, otempdegC"<<std::endl;
  for(unsigned int i=0;i<nodes.size();++i)
    out<<run_num           <<","
       <<nodes[i].emergence<<","
       <<i                 <<","
       <<nodes[i].lastchild<<","
       <<i                 <<","
       <<nodes[i].otempdegC<<std::endl;
}


//Print the phylogeny as a Newick tree
std::string Phylogeny::printNewick(int n, int depth) const {
  if(nodes[n].children.empty()){
    double length=nodes[n].lastchild-nodes[n].emergence;
    return "S"+std::to_string(n) + ":" + std::to_string(length);
  }

  //Copy this node's list of children so that we can guarantee that this method
  //does not alter the phylogeny data
  std::vector<int> my_children(nodes[n].children.begin(), nodes[n].children.end());

  //Reverse the order of the children so that the list of children is order from
  //most to least recent
  std::reverse(my_children.begin(),my_children.end());

  //This string holds the phylogeny of the node as we build it. Since the node
  //may have several children which emerged at different times we assume that
  //the parent also becomes a new species at that time and inject "virtual
  //nodes" into the output of the phylogeny. This string holds the virtual nodes
  //as we build them.
  std::string my_phylo="";

  for(unsigned int c=0;c<my_children.size();c++){
    //Alias for the current child
    int child=my_children[c];

    //Avoid an infinite loop with Eve
    if(n==child) continue;

    //Holds the returned Newick representation of the phylogeny rooted at the
    //current child
    std::string childphylo=printNewick(child,depth+1);

    //If this is the most recent child the parent has had, then the parent's
    //lineage ends when the parent has its last child. Thus, the phylogeny is of
    //the form: (ParentName:TimeToParentsLastChild,ChildName:TimeToChildsLastChild)
    if(c==0){
      my_phylo="(S";
      my_phylo+=std::to_string(n)+":";
      my_phylo+=std::to_string(nodes[n].lastchild-nodes[child].emergence);
      my_phylo+=",";
      my_phylo+=childphylo;
      my_phylo+=")";
    } else {
    //This isn't the most recent child the parent has had. The current child has
    //a phylogeny and the more recent child has a phylogeny as well. Thus, we
    //need to weld these two phylogenies together
      std::string temp="(";
      temp+=my_phylo+":";
      //The previous child returned a phylogeny including a virtual node. This
      //virtual node indicates a branching that took place when the previous
      //child split from the parent line. This branching took place at the point
      //in time at which the previous child emerged. So the length to that node
      //is from the emergence time of this child (the time of this splitting) to
      //the emergence time of the previous child. Note: If the previous child
      //split at the same time this child did, then there is a trifurcation. We
      //do not explicitly handle this case.
      int prevchild=my_children[c-1];
      temp+=std::to_string(nodes[prevchild].emergence-nodes[child].emergence);
      temp+=",";
      temp+=childphylo;
      //Since "childphylo" is a valid phylogeny and includes a distance to
      //root, it is not necessary to specify a distance here
      temp+=")";
      my_phylo=temp;
    }
  }

  //This is the distance between the least recent child's emergence and the
  //emergence of the parent. Having this here ensures that every time this
  //function returns the result is a valid phylogeny including a distance to
  //root
  my_phylo+=":"+std::to_string(nodes[my_children.back()].emergence-nodes[n].emergence);

  //All done: let's wrap things up nicely
  if(depth==0)
    my_phylo+=";";

  return my_phylo;
}


//For each node in the phylogenetic tree, print the summary statistics of that
//species throughout the duration of its existence
void Phylogeny::speciesSummaries(int run_num, std::ofstream &out) const {
  if(run_num==0)
    out<<"RunNum, Species, Time, NumAlive, ElevMin, ElevMax, ElevAvg, TempMin, TempMax, TempAvg\n";
  for(unsigned int i=0;i<nodes.size();++i)
  for(auto &ss: nodes[i].stats){  //There is a stats record of each time the species was alive
     out<<run_num                        <<","
        <<i                              <<","
        <<ss.t                           <<","
        <<ss.num_alive                   <<","
        <<ss.elev_min                    <<","
        <<ss.elev_max                    <<","
        <<(ss.elev_avg/ss.num_alive)     <<","
        <<ss.opt_temp_min                <<","
        <<ss.opt_temp_max                <<","
        <<(ss.opt_temp_avg/ss.num_alive) << std::endl;
  }
}