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

// CMS3
#include "CMS3_WWW0116.h"

#ifndef FUNCTION_H
#define FUNCTION_H

struct myevt{
  unsigned int run;
  unsigned int ls;
  long long evt;
};

#endif /* FUNCTION_H */

#define MZ 91.1876

vector<string> split(string mystring, string c="/");
int gentype_v2(unsigned lep1_index=0,unsigned lep2_index=1, int lep3_index=-1);
float dR(LorentzVector vec1,LorentzVector vec2 );
float dEta(LorentzVector vec1,LorentzVector vec2 );
float dPhi(LorentzVector vec1,LorentzVector vec2 );
float deltaPhi(float phi1,float phi2 );
float mT(LorentzVector p4, LorentzVector met);
bool sortMSFOS(float M1, float M2);
int numJ(float jetpt=20., float jeteta=2.4, float csv=-1, int jec=0);//loose: 0.5426, medium: 0.8484, if csv < 0 - count Njets, else Nbjets
bool getalljetnumbers(int &nj, int &nj30, int &nb, int jec=0);
bool getleptonindices(vector<int> &iSS, vector<int> &i3l, vector<int> &iaSS, vector<int> &ia3l, vector<int> &vSS, vector<int> &v3l, vector<int> &vaSS, vector<int> &va3l);
float coneCorrPt(int lepi);//apply only to loose but not tight leptons
float loadFR(float &FRSSerr,int index, TH2D *hMuFR, TH2D *hElFR, float muFRptmin=10.1, float muFRptmax=119.9, float muFRetamin=0.01, float muFRetamax=2.39, float elFRptmin=10.1, float elFRptmax=119.9, float elFRetamin=0.01, float elFRetamax=2.49, bool conecorrected=true);
float calcMjj(bool closestDR=true, int jec=0);//Mjj or MjjL
float Detajj(int jec=0);
bool getMjjAndDeta(float &Mjj, float &MjjL, float &Detajj, int jec=0);
float calcMTmax(vector<int> lepindex, LorentzVector MET, bool compareSSpairs=false);
float calcMTmax(int index1, int index2, LorentzVector MET);
bool passofflineTriggers(vector<int> tightlep, vector<int> looselep={});
bool passonlineTriggers(vector<int> tightlep, vector<int> looselep={});//I'll duplicate and goodrun selection by hand - I think this way it is easier
bool passFilters();
bool splitVH(string filename);//true is WHtoWWW, false everything else
string whatWW(string filename, string samplename);
string process(string filename, bool SS, vector<int> tightlep, vector<int> looselep={});
bool vetophotonprocess(string filename, bool process);//true if W,DY,tt,etc. is photonfake; true if Wg,Zg,ttg,etc. is not photonfake. else false.
int isSRSS(vector<int> tightlep, vector<int> vetolep,                        bool preselect=false, float maxMT=-1, int nj=-1, int nb=-1, float Mjj=-1, float MjjL=-1, float Detajj=-1, LorentzVector MET = LorentzVector(0,0,0,0), int jec=0, bool btag=false, bool Mjjside=false);
int isARSS(vector<int> tightlep, vector<int> looselep, vector<int> vetolep,  bool preselect=false, float maxMT=-1, int nj=-1, int nb=-1, float Mjj=-1, float MjjL=-1, float Detajj=-1, LorentzVector MET = LorentzVector(0,0,0,0), int jec=0, bool btag=false, bool Mjjside=false);
int isCRSS(vector<int> tightlep, vector<int> selectlep, vector<int> vetolep, bool preselect=false, float maxMT=-1, int nj=-1, int nb=-1, float Mjj=-1, float MjjL=-1, float Detajj=-1, LorentzVector MET = LorentzVector(0,0,0,0), int jec=0, bool noZ=false, bool btag=false, bool Mjjside=false);
bool passJetSSstate(bool preselect=false, int nj=-1, int nb=-1, float Mjj=-1, float MjjL=-1, float Detajj=-1, bool is3lCR=false, int jec=0, bool btag=false, bool Mjjside=false);
int isSR3l(vector<int> tightlep,                                     bool preselect=false, int nj=-1, int nb=-1, LorentzVector MET = LorentzVector(0,0,0,0), int jec=0, bool btag=false);
int isAR3l(vector<int> tightlep, vector<int> looselep,               bool preselect=false, int nj=-1, int nb=-1, LorentzVector MET = LorentzVector(0,0,0,0), int jec=0, bool btag=false);
int isCR3l(vector<int> tightlep,                                     bool preselect=false, int nj=-1, int nb=-1, LorentzVector MET = LorentzVector(0,0,0,0), int jec=0, bool btag=false);
bool checkbothSRCR3l(int &isSR3l, int &isCR3l, vector<int> tightlep, bool preselect=false, int nj=-1, int nb=-1, LorentzVector MET = LorentzVector(0,0,0,0), int jec=0, bool btag=false);
vector<float> allMSFOS(vector<int> tightlep, vector<int> looselep={});
map<string, TH1D*> bookhistograms(string samplename, vector<string> histonames, vector<int> hbins, vector<float> hlow, vector<float> hup, TDirectory *rootdir, int splitWW=0);
bool fillSRhisto(map<string, TH1D*> histos, string histoname, string sn, string sn2, int SRSS, int SR3l, float weight, float weight3l=-2e6);
bool SaveHistosToFile(string filename, map<string, TH1D*> histos, bool addunderflow=true, bool addoverflow=true);

bool fileexists(string filename);
int loadlepSFfile(TFile *&f, TH2F *&hMuID, TH2F *&hMutrigger, TH2F *&hElID, TH2F *&hEltrigger, string filename="rootfiles/SF_TnP.root", string muIDname="muSF", string mutrigname="", string elIDname="elSFreco", string eltrigname="elSF_ID");
bool loadFakeRates(TFile *&f, TH2D *&hMuFR, TH2D *&hElFR, string filename="rootfiles/fakerate_pt_v_eta.root", string muname="muon_fakerate_conecorrpt_v_eta", string elname="elec_fakerate_conecorrpt_v_eta");
bool deleteFiles(bool SF, TFile *&fSF);

float getlepSFandError(float &error, int index, TH2F *hMuID, TH2F *hMutrigger, TH2F *hElID, TH2F *hEltrigger, float muIDptmin=20.1,float muIDptmax=199.9, float muIDetamin=0.01, float muIDetamax=2.49, float muTrptmin=20.1,float muTrptmax=199.9, float muTretamin=0.01, float muTretamax=2.49, float elIDptmin=10.1,float elIDptmax=199.9, float elIDetamin=0.01, float elIDetamax=2.49, float elTrptmin=25.1,float elTrptmax=199.9, float elTretamin=0.01, float elTretamax=2.49);
float getlepSFWeightandError(float &error, vector<int> tightlep, vector<int> looselep, TH2F *hMuID, TH2F *hMutrigger, TH2F *hElID, TH2F *hEltrigger, float muIDptmin=20.1,float muIDptmax=199.9, float muIDetamin=0.01, float muIDetamax=2.49, float muTrptmin=20.1,float muTrptmax=199.9, float muTretamin=0.01, float muTretamax=2.49, float elIDptmin=10.1,float elIDptmax=199.9, float elIDetamin=0.01, float elIDetamax=2.49, float elTrptmin=25.1,float elTrptmax=199.9, float elTretamin=0.01, float elTretamax=2.49);//need to feed looselep here :/
float getlepSFWeightandError(float &error, vector<float> efftight, vector<float> errtight, vector<float> effloose={}, vector<float> errloose={});
bool addeventtocheck(vector<myevt> &eventvector, unsigned int runnumber, unsigned int lumisection, long long eventnumber);
bool checkthisevent(vector<myevt> eventvector, unsigned int runnumber, unsigned int lumisection, long long eventnumber);
void addeventtolist(std::ostringstream *&stream);
void addeventtolist(int SS, std::ostringstream *&streamEE, std::ostringstream *&streamEM, std::ostringstream *&streamMM);
void addeventtolist(int SRSS, int SR3l, std::ostringstream *&streamEE, std::ostringstream *&streamEM, std::ostringstream *&streamMM, std::ostringstream *&stream0SFOS, std::ostringstream *&stream1SFOS, std::ostringstream *&stream2SFOS);
void storeeventlist(string output, string sample, std::ostringstream *&stream);
void storeeventlist(string output, string sample, bool isSS, std::ostringstream *&streamEE, std::ostringstream *&streamEM, std::ostringstream *&streamMM);
void storeeventlist(string output, string sample, std::ostringstream *&streamEE, std::ostringstream *&streamEM, std::ostringstream *&streamMM, std::ostringstream *&stream0SFOS, std::ostringstream *&stream1SFOS, std::ostringstream *&stream2SFOS);
