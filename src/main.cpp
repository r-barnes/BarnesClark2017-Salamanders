#include "salamander.hpp"
#include "mtbin.hpp"
#include "phylo.hpp"
#include <array>
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

//Profile of the mountain
//We allocate this in the global namespace to prevent a stack overflow
static array< MtBin, 1000> mts;

typedef std::vector<Phylo> phylolist;

void UpdatePhylogeny(double t, phylolist &plist){
  for(auto &m: mts)
  for(auto &s: m.bin){
    //If I am dead or have no parent, skip me
    if(s.dead || s.parent==-1) continue;

    //If I am not similar to my parent
    if(!s.pSimilarGenome(plist.at(s.parent).genes)) {
      bool trigger=false;

      //See if I am similar to any other salamanders, starting with the most recent
      for(int p=plist.size()-1;p>=0;--p) {
        if( (t-plist.at(s.parent).lastchild)>=1 ) continue; //TODO: This needs to be adjusted based on the length of timesteps. This is crappy.
        //If my parent is the same as this salamander's parent and my genes are similar to this salamander's
        //Then this salamander and I are both part of the first generation of a new species of salamander
        //Therefore, I will consider myself this salamander's child, since its genome is already stored in the phylogeny
        if(s.parent==plist.at(p).parent && s.pSimilarGenome(plist.at(p).genes) ) {
          s.parent=p;
          trigger=true;
          break;
        }
      }
      if(!trigger)                       //No salamander in the phylogeny was similar to me!
        plist.push_back(Phylo(s, t));    //Therefore, I add myself to the phylogeny as a new species
    } else {  //I am similar to my parent, so mark my parent (species) as having survived this long
      plist.at(s.parent).lastchild=t;
    }
  }
}

int main(){
  phylolist phylos;

  ////////////////////////////////////
  //INITIALIZE
  ////////////////////////////////////

  //Start off by killing everything
  for(auto &m: mts)
    m.killAll();

  Salamander Eve;
  Eve.parent=0;
  Eve.otemp =33.5618604122814;//degC //This is the global average temperature at sea level 65 million years ago.
                                     //Today, the optimal temperature for salamanders is 12.804279 degC
  Eve.genes =(Salamander::genetype)169113934842922489;
  Eve.dead  =false;

  for(unsigned int m=0;m<1;++m)
  for(unsigned int s=0;s<10;++s){
    mts[m].addSalamander(Eve);
  }
  phylos.push_back(Phylo(Eve, 0));





  ////////////////////////////////////
  //MAIN LOOP
  ////////////////////////////////////


  //Find cell corresponding to optimal happy temperature
  //Make a salamander in that cell alive
  //Go forth

  //Salamanders are at their very happiest at 12.804279 degC
  //At sea level 65MYA the temp was 33.5618604122814 degC
  //Temperature drops at 9.8 degC/1000m (adiabatic lapse rate of dry air)
  //So (33.5618604122814-12.804279)/9.8=2.1181*1000m=2.11km

  for(double t=0;t<65.001;t+=0.5){
  	cerr<<"#"<<t<<endl;
    unsigned int population_size=0;
    for(unsigned int m=0;m<mts.size();++m){
      mts[m].mortaliate(t);
      mts[m].breed(t);
      if(m>0)            mts[m].diffuse(t,mts[m-1]);
      if(m<mts.size()-1) mts[m].diffuse(t,mts[m+1]);
      population_size+=mts[m].startofdead;
    }
    UpdatePhylogeny(t, phylos);

    cerr<<"#Species count="<<phylos.size()<<endl;
    cerr<<"#Population size="<<population_size<<endl;
  }



  ////////////////////////////////////
  //OUTPUT
  ////////////////////////////////////
  ofstream phylograph("output/phylograph.dot");
  phylograph<<"digraph graphname {"<<endl;
  for(unsigned int i=0;i<phylos.size();++i)
    phylograph<<i<<"[label=\""<<phylos[i].emergence<<" "<<phylos[i].otemp<<"\"];"<<endl;
  for(unsigned int i=0;i<phylos.size();++i)
    phylograph<<phylos[i].parent<<"->"<<i<<";"<<endl;
  phylograph<<"}"<<endl;
  phylograph.close();

  ofstream persistgraph("output/persistgraph.csv");
  for(unsigned int i=0;i<phylos.size();++i)
    persistgraph<<phylos[i].emergence<<","<<i<<","<<phylos[i].lastchild<<","<<i<<endl;
  persistgraph.close();

}