// Usage:
// > root -b doAll.C

// C++
#include <iostream>
#include <vector>
#include <map>

// ROOT
#include "TBenchmark.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH3D.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TTreeCache.h"
#include "TLorentzVector.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <tuple>

#define USE_CMS3_WWW100 
// CMS3
#ifdef USE_CMS3_WWW100
#include "CMS3.h"
#else
#include "CMS3_WWW0118.h"
#endif

#ifndef FUNCTION_H
#define FUNCTION_H

struct myevt{
  unsigned int run;
  unsigned int ls;
  long long evt;
};

#endif /* FUNCTION_H */

#define MZ 91.1876

float  dR(LorentzVector vec1,LorentzVector vec2 );
float  dEta(LorentzVector vec1,LorentzVector vec2 );
float  dPhi(LorentzVector vec1,LorentzVector vec2 );
float  deltaPhi(float phi1,float phi2 );
float  mT(LorentzVector p4, LorentzVector met);
int    numJ(float jetpt=20., float jeteta=2.4, float csv=-1, int jec=0);//loose: 0.5426, medium: 0.8484, if csv < 0 - count Njets, else Nbjets

bool  sortMSFOS(float M1, float M2);
bool  sortPt(int i1, int i2);
bool  sortDecreasing(double i1, double i2);

vector<string> split(string mystring, string c="/");
string leptype(unsigned lep_index);
int    gentype_v2(unsigned lep1_index=0,unsigned lep2_index=1, int lep3_index=-1);
bool   splitVH(string filename);//true is WHtoWWW, false everything else
string whatWW(string filename, string samplename);
string process(string filename, bool SS, vector<int> tightlep, vector<int> looselep={});//not needed for new
bool   vetophotonprocess(string filename, bool process);//true if W,DY,tt,etc. is photonfake; true if Wg,Zg,ttg,etc. is not photonfake. else false.

bool  getleptonindices_v0(vector<int> &iSS, vector<int> &i3l, vector<int> &iaSS, vector<int> &ia3l, vector<int> &vSS, vector<int> &v3l, vector<int> &vaSS, vector<int> &va3l, int version=0, float pTSS=30., float pT3l=20.);//deprecated
#ifdef USE_CMS3_WWW100
bool  getleptonindices_v2(vector<int> &iSS, vector<int> &i3l, vector<int> &iaSS, vector<int> &ia3l, vector<int> &vSS, vector<int> &v3l, vector<int> &vaSS, vector<int> &va3l, float pTSS=25., float pT3l=20.);
int   NtightSS(float pTSS=25);
int   NlooseSS(float pTSS=25);
int   Ntight3l(float pT3l=20);
int   Nloose3l(float pT3l=20);
#endif
bool  istightlepton(int i, int version);
bool  islooselepton(int i, int version);
float coneCorrPt(int lepi, int version=0);//apply only to loose but not tight leptons

bool  getalljetnumbers(int &nj, int &nj30, int &nb, int jec=0);
float calcMjj(bool closestDR=true, int jec=0);//Mjj or MjjL
float Detajj(int jec=0);
bool  getMjjAndDeta(float &Mjj, float &MjjL, float &Detajj, int jec=0);

float calcMTmax(vector<int> lepindex, LorentzVector MET, bool compareSSpairs=false);
float calcMTmax(int index1, int index2, LorentzVector MET);

bool  passTriggers(bool MCoffline=true, bool MCtrigger=false, float leadleppt=25., float trailleppt=20.);
bool  passofflineTriggers(vector<int> tightlep, vector<int> looselep={});
bool  passonlineTriggers(vector<int> tightlep, vector<int> looselep={});//I'll duplicate and goodrun selection by hand - I think this way it is easier
bool  passFilters_v0();
bool  passFilters();

bool  passJetSSstate_v0(bool preselect=false, int nj=-1, int nb=-1, float Mjj=-1, float MjjL=-1, float Detajj=-1, bool is3lCR=false, int jec=0, bool btag=false, bool Mjjside=false, int version=0);
bool  passJetSSstate(bool preselect=false, bool is3lCR=false, int jec=0, bool btag=false, bool Mjjside=false, int version=0);
int   isSRSS_v0(vector<int> tightlep, vector<int> vetolep,                        bool preselect=false, float maxMT=-1, int nj=-1, int nb=-1, float Mjj=-1, float MjjL=-1, float Detajj=-1, LorentzVector MET = LorentzVector(0,0,0,0), int jec=0, bool btag=false, bool Mjjside=false, int version=0);
int   isARSS_v0(vector<int> tightlep, vector<int> looselep, vector<int> vetolep,  bool preselect=false, float maxMT=-1, int nj=-1, int nb=-1, float Mjj=-1, float MjjL=-1, float Detajj=-1, LorentzVector MET = LorentzVector(0,0,0,0), int jec=0, bool btag=false, bool Mjjside=false, int version=0);
int   isCRSS_v0(vector<int> tightlep, vector<int> selectlep, vector<int> vetolep, bool preselect=false, float maxMT=-1, int nj=-1, int nb=-1, float Mjj=-1, float MjjL=-1, float Detajj=-1, LorentzVector MET = LorentzVector(0,0,0,0), int jec=0, bool noZ=false, bool btag=false, bool Mjjside=false, int version=0);
int    isSRSS(                              bool preselect=false, int jec=0,                 bool btag=false, bool Mjjside=false, int version=1);
int    isARSS(                              bool preselect=false, int jec=0,                 bool btag=false, bool Mjjside=false, int version=1);
int    isCRSS(                              bool preselect=false, int jec=0, bool noZ=false, bool btag=false, bool Mjjside=false, int version=1);
bool   passAnySS(int &SR, int &AR, int &CR, bool preselect=false, int jec=0, bool noZ=false, bool btag=false, bool Mjjside=false, int version=1);

bool   passJet3lstate(int jec=0);
int    isSR3l_v0(vector<int> tightlep,                                     bool preselect=false, int nj=-1, int nb=-1, LorentzVector MET = LorentzVector(0,0,0,0), int jec=0, bool btag=false, int version=0);
int    isAR3l_v0(vector<int> tightlep, vector<int> looselep,               bool preselect=false, int nj=-1, int nb=-1, LorentzVector MET = LorentzVector(0,0,0,0), int jec=0, bool btag=false, int version=0);
int    isCR3l_v0(vector<int> tightlep,                                     bool preselect=false, int nj=-1, int nb=-1, LorentzVector MET = LorentzVector(0,0,0,0), int jec=0, bool btag=false, int version=0);
bool    checkbothSRCR3l_v0(int &isSR3l, int &isCR3l, vector<int> tightlep, bool preselect=false, int nj=-1, int nb=-1, LorentzVector MET = LorentzVector(0,0,0,0), int jec=0, bool btag=false, int version=0);
vector<float> allMSFOS(vector<int> tightlep, vector<int> looselep={});
int    isSR3l(                              bool preselect=false, int jec=0, bool btag=false, int version=1);
int    isAR3l(                              bool preselect=false, int jec=0, bool btag=false, int version=1);
int    isCR3l(                              bool preselect=false, int jec=0, bool btag=false, int version=1);
bool   passAny3l(int &SR, int &AR, int &CR, bool preselect=false, int jec=0, bool btag=false, int version=1);

map<string, TH1D*> bookhistograms(string samplename, vector<string> histonames, vector<int> hbins, vector<float> hlow, vector<float> hup, TDirectory *rootdir, int splitWW=0);
map<string, TH1D*> bookhistclosure(string samplename, vector<string> histonames, vector<int> hbins, vector<float> hlow, vector<float> hup, TDirectory *rootdir);
bool  deleteHistograms(map<string, TH1D*> histos);
bool  fillSRhisto(map<string, TH1D*> histos, string histoname, string sample, string sn, string sn2, int SRSS, int SR3l, float weight, float weight3l=-2e6, bool fillsample=true);
bool  fillhisto(map<string, TH1D*> histos, string histoname, string sample, string sn, float value, float weight, bool fillsample=true);
bool  SaveHistosToFile(string filename, map<string, TH1D*> histos, bool addunderflow=true, bool addoverflow=true, bool deletefile=false);
bool  fileexists(string filename);

float getlepSFWeightandError(float &error, vector<int> tightlep, vector<int> looselep={});//loose not implemented yet
float getlepFakeRateandError(float &error, int index, bool data=true, bool conecorr=true, int version=0);
float getlepFRWeightandError(float &error, int index, bool data=true, bool conecorr=true, bool addclosureunc=false, int version=0);
float getlepFRClosureError(int index, bool data=true, bool conecorr=true, int version=0);
float getTriggerWeightandError(float &error, vector<int> tightlep, vector<int> looselep={});

float getPUWeight(int variation=0);
float getPUWeightAndError(float &errorup, float &errordown);

bool addeventtocheck(vector<myevt> &eventvector, unsigned int runnumber, unsigned int lumisection, long long eventnumber);
bool checkthisevent(vector<myevt> eventvector, unsigned int runnumber, unsigned int lumisection, long long eventnumber);
void addeventtolist(std::ostringstream *&stream);
void addeventtolist(int SS, std::ostringstream *&streamEE, std::ostringstream *&streamEM, std::ostringstream *&streamMM);
void addeventtolist(int SRSS, int SR3l, std::ostringstream *&streamEE, std::ostringstream *&streamEM, std::ostringstream *&streamMM, std::ostringstream *&stream0SFOS, std::ostringstream *&stream1SFOS, std::ostringstream *&stream2SFOS);
void storeeventlist(string output, string sample, std::ostringstream *&stream);
void storeeventlist(string output, string sample, bool isSS, std::ostringstream *&streamEE, std::ostringstream *&streamEM, std::ostringstream *&streamMM);
void storeeventlist(string output, string sample, std::ostringstream *&streamEE, std::ostringstream *&streamEM, std::ostringstream *&streamMM, std::ostringstream *&stream0SFOS, std::ostringstream *&stream1SFOS, std::ostringstream *&stream2SFOS);

// Calculating SFOS pair number and invariant masses
int   calcNSFOS(std::vector<int> tightlep);
float get0SFOSMll(std::vector<int> lepidx);
float get0SFOSMee(std::vector<int> lepidx);
float get1SFOSMll(std::vector<int> lepidx);
float get2SFOSMll0(std::vector<int> lepidx);
float get2SFOSMll1(std::vector<int> lepidx);

// Getting the truth info on the signal samples. (Verified only on the WWW samples. (not WH))
bool isW(int igen);
bool isParentW(int igen);
bool isParentTau(int igen);
bool isLepton(int igen);
bool isQuark(int igen);
bool isLightLepton(int igen);
bool isTau(int igen);
typedef std::vector<int> vIdx;
std::tuple<vIdx, vIdx> getGenIndices();
