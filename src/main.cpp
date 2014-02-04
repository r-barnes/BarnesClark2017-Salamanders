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

int main(){
  Phylogeny phylos;
  vector<MtBin> mts;
  mts.resize(1000);

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
  phylos.nodes.push_back(PhyloNode(Eve, 0));





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
    phylos.UpdatePhylogeny(t, mts);

    cerr<<"#Species count="<<phylos.nodes.size()<<endl;
    cerr<<"#Population size="<<population_size<<endl;
  }



  ////////////////////////////////////
  //OUTPUT
  ////////////////////////////////////
  ofstream phylograph("output/phylograph.dot");
  phylograph<<"digraph graphname {"<<endl;
  for(unsigned int i=0;i<phylos.nodes.size();++i)
    phylograph<<i<<"[label=\""<<phylos.nodes[i].emergence<<" "<<phylos.nodes[i].otemp<<"\"];"<<endl;
  for(unsigned int i=0;i<phylos.nodes.size();++i)
    phylograph<<phylos.nodes[i].parent<<"->"<<i<<";"<<endl;
  phylograph<<"}"<<endl;
  phylograph.close();

  ofstream persistgraph("output/persistgraph.csv");
  for(unsigned int i=0;i<phylos.nodes.size();++i)
    persistgraph<<phylos.nodes[i].emergence<<","<<i<<","<<phylos.nodes[i].lastchild<<","<<i<<endl;
  persistgraph.close();

}
