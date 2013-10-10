#include "salamander.hpp"
#include "mtbin.hpp"
#include "phylo.hpp"
#include <array>
#include <vector>
#include <iostream>
using namespace std;

//Profile of the mountain
//We allocate this in the global namespace to prevent a stack overflow
static array< MtBin, 1000> mts;

typedef std::vector<Phylo> phylolist;

void UpdatePhylogeny(double t, phylolist &plist){
  for(auto &m: mts)
  for(auto &s: m.bin){
    if(s.dead || s.parent==-1) continue;

    if(!s.pSimilarGenome(plist.at(s.parent).genes)) {
      bool trigger=false;
      for(int p=plist.size()-1;p>=0;--p) {
        if(s.parent==plist.at(p).parent && s.pSimilarGenome(plist.at(p).genes) ) {
          s.parent=p;
          trigger=true;
          break;
        }
      }
      if(!trigger)
        plist.push_back(Phylo(s, t));
    } else {
      plist.at(s.parent).lastchild=t;
    }
  }
}

void OutputPhylo(phylolist &plist){
  cout<<"digraph graphname {"<<endl;

  for(unsigned int i=0;i<plist.size();++i)
    cout<<i<<"[label=\""<<plist[i].emergence<<" "<<plist[i].otemp<<"\"];"<<endl;

  for(unsigned int i=0;i<plist.size();++i)
    cout<<plist[i].parent<<"->"<<i<<";"<<endl;

  cout<<"}"<<endl;
}


int main(){
  phylolist phylos;

  //Start off by killing everything
  for(auto &m: mts)
    m.killAll();

  Salamander Eve;
  Eve.parent=0;
  Eve.otemp =12.804279;
  Eve.genes =(Salamander::genetype)169113934842922489;
  Eve.dead  =false;

  for(unsigned int m=0;m<mts.size();++m)
  for(unsigned int s=0;s<10/*mts.bin.size()*/;++s){
    mts[m].addSalamander(Eve);
  }
  phylos.push_back(Phylo(Eve, 0));

  //Find cell corresponding to optimal happy temperature
  //Make a salamander in that cell alive
  //Go forth

  //Salamanders are at their very happiest at 12.804279 degC
  //At sea level 65MYA the temp was 33.5618604122814 degC
  //Temperature drops at 9.8 degC/1000m (adiabatic lapse rate of dry air)
  //So (33.5618604122814-12.804279)/9.8=2.1181*1000m=2.11km

  for(double t=0;t<28;t+=0.5){
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

  OutputPhylo(phylos);
}
