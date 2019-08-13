#include "Functions_test.h"
#include "LeptonAndTriggerSF.h"

using namespace std;
using namespace tas;

vector<string> split(string mystring, string c) {
  string s = mystring;
  vector<string> v;
  string::size_type i = 0;
  string::size_type j = s.find(c);

  while (j != string::npos) {
    v.push_back(s.substr(i, j-i));
    i = ++j;
    j = s.find(c, j);
    
    if (j == string::npos)
      v.push_back(s.substr(i, s.length()));
   }
  return v;
}

string leptype(unsigned lep_index){
  if(      lep_motherIdSS().at(lep_index) == 1)  return "real";
  else if( lep_motherIdSS().at(lep_index) == 2)  return "chargeflip";
  else if( lep_motherIdSS().at(lep_index) ==-1)  return "bfake";
  else if( lep_motherIdSS().at(lep_index) ==-2)  return "cfake";
  else if( lep_motherIdSS().at(lep_index) ==-3)  return "gamma";
  else if( lep_motherIdSS().at(lep_index) == 0)  return "lightfake";
  else                                           return "unknown";
}

int gentype_v2(unsigned lep1_index,unsigned lep2_index, int lep3_index){
  bool gammafake = false;
  bool jetfake   = false;
  unsigned int ngenlep = ngenLepFromTau()+ ngenLep();
  unsigned int nW(0), nZ(0);
  bool lep1_real = lep_motherIdSS().at(lep1_index) > 0;
  bool lep2_real = lep_motherIdSS().at(lep2_index) > 0;
  bool lep3_real = false;
  if(lep3_index>=0) lep3_real = lep_motherIdSS().at(lep3_index) > 0;
  vector<int> reallepindex;

  for (unsigned int lepindex = 0;lepindex<lep_p4().size();++lepindex){
      if(lep_motherIdSS().at(lepindex) > 0) reallepindex.push_back(lepindex);
      else if(lep_motherIdSS().at(lepindex) == -3) gammafake = true;
      else                                           jetfake = true;
      if(lep_isFromW().at(lepindex)) nW++;
      if(lep_isFromZ().at(lepindex)) nZ++;
  }
  //found two real leptons
  if(lep3_index<0){
    bool ischargeflip = false;
    bool isSS = false;
    if(lep1_real&&lep2_real) {
      int ilep1 =   lep_genPart_index().at(lep1_index);
      int ilep2 =   lep_genPart_index().at(lep2_index);
      bool lep1_chargeflip  =genPart_charge().at(ilep1)!= lep_charge().at(lep1_index);
      bool lep2_chargeflip  =genPart_charge().at(ilep2)!= lep_charge().at(lep2_index);
      if (     !lep1_chargeflip&&!lep2_chargeflip&&nW==2) return 0; // true SS
      else if (!lep1_chargeflip&&!lep2_chargeflip)        isSS = true; // true SS - but could be still lost lepton WZ
      if (      lep1_chargeflip|| lep2_chargeflip)        ischargeflip = true; // charge flip
    }
    
    if(ngenlep>2 || reallepindex.size()>2 || (nW>0 && nZ>0)) return 3; // lostlep
    if((ngenlep<2 ||!lep1_real||!lep2_real)&&       jetfake) return 4; // jetfake - if double fake with one jet fake and one gamma fake call it jet fake
    if((ngenlep<2 ||!lep1_real||!lep2_real)&&     gammafake) return 5; // gammafake
    if((ngenlep<2 ||!lep1_real||!lep2_real)&&    !gammafake) return 4; // call all without gamma fake jetfake - safety cut
    if(isSS)         return 0;
    if(ischargeflip) return 2;
    
    cout << "This event was not classified - 2 lepton event - v2" << endl;
    return 1;
  } else {
    //found three real leptons
    bool ischargeflip = false;
    bool isthreelep = false;
    if(lep1_real&&lep2_real&&lep3_real) {
      int ilep1 =   lep_genPart_index().at(lep1_index);
      int ilep2 =   lep_genPart_index().at(lep2_index);
      int ilep3 =   lep_genPart_index().at(lep3_index);
      bool lep1_chargeflip  =genPart_charge().at(ilep1)!= lep_charge().at(lep1_index);
      bool lep2_chargeflip  =genPart_charge().at(ilep2)!= lep_charge().at(lep2_index);
      bool lep3_chargeflip  =genPart_charge().at(ilep3)!= lep_charge().at(lep3_index);
      if (     !lep1_chargeflip&&!lep2_chargeflip&&!lep3_chargeflip&&nW==3) return 0; // true WWW
      else if (!lep1_chargeflip&&!lep2_chargeflip&&!lep3_chargeflip)        isthreelep = true; // true 3l, but could be lost lepton ZZ
      if (      lep1_chargeflip|| lep2_chargeflip|| lep3_chargeflip)        ischargeflip = true; // charge flip
    }
    if(ngenlep>3 || reallepindex.size()>3 || (nW>=2 && nZ>=1) || (nZ>=3)) return 3; // lostlep (2 lep from W and 2 from Z, or 4 from Z)
    //there is the case of having WZZ with two lost leptons --> ngenlep>3 - correctly put has lostlep
    if((ngenlep<3 ||!lep1_real||!lep2_real||!lep3_real)&&        jetfake) return 4; // jetfake
    if((ngenlep<3 ||!lep1_real||!lep2_real||!lep3_real)&&      gammafake) return 5; // gammafake
    if((ngenlep<3 ||!lep1_real||!lep2_real||!lep3_real)&&     !gammafake) return 4; // jetfake
    if(isthreelep)   return 1;
    if(ischargeflip) return 2;
    
    cout << "This event was not classified - 3 lepton event - v2" << endl;
    return 0;
  }
}

float dR(LorentzVector vec1,LorentzVector vec2 ){
  float dphi = std::min(float(fabs(vec1.Phi() - vec2.Phi())), float(2 * M_PI - fabs(vec1.Phi() - vec2.Phi())) );
  float deta = vec1.Eta() - vec2.Eta();
  return sqrt(dphi*dphi + deta*deta);
}

float dEta(LorentzVector vec1,LorentzVector vec2 ){
  return fabs(vec1.Eta() - vec2.Eta());
}

float dPhi(LorentzVector vec1,LorentzVector vec2 ){
  return fabs(std::min(float(fabs(vec1.Phi() - vec2.Phi())), float(2 * M_PI - fabs(vec1.Phi() - vec2.Phi())) ));
}

float deltaPhi(float phi1,float phi2 ){
  return fabs(std::min(float(fabs(phi1-phi2)),float(2*M_PI-fabs(phi1-phi2))));
}

float mT(LorentzVector p4, LorentzVector met){
  float phi1 = p4.Phi();
  float phi2 = met.Phi();
  float Et1  = p4.Et();
  float Et2  = met.Et();
  return sqrt(2*Et1*Et2*(1.0 - cos(phi1-phi2)));
}
bool sortMSFOS(float M1, float M2){
  return fabs(M1-MZ)<fabs(M2-MZ);
}
bool sortPt(int i1, int i2){
  return lep_p4()[i1].Pt()>=lep_p4()[i2].Pt();
}
bool sortDecreasing(double i1, double i2){
  return i1>=i2;
}

int numJ(float jetpt, float jeteta, float csv, int jec){
  int num = 0;
  if(jec==1){
    for(unsigned int n = 0; n<jets_up_csv().size();++n){
      if(fabs(jets_up_p4()[n].Eta())>jeteta) continue;
      if(     jets_up_p4()[n].Pt()  <jetpt)  continue;
      if(csv<0)                     ++num;
      else if(jets_up_csv()[n]>csv) ++num;
    }
  }
  else if(jec==(-1)){
    for(unsigned int n = 0; n<jets_dn_csv().size();++n){
      if(fabs(jets_dn_p4()[n].Eta())>jeteta) continue;
      if(     jets_dn_p4()[n].Pt()  <jetpt)  continue;
      if(csv<0)                     ++num;
      else if(jets_dn_csv()[n]>csv) ++num;
    }
  }
  else{ 
    for(unsigned int n = 0; n<jets_csv().size();++n){
      if(fabs(   jets_p4()[n].Eta())>jeteta) continue;
      if(        jets_p4()[n].Pt()  <jetpt)  continue;
      if(csv<0)                     ++num;
      else if(   jets_csv()[n]>csv) ++num;
    }
  }
  return num;
}

bool getalljetnumbers(int &nj, int &nj30, int &nb, int jec){
  int num  = 0;
  int num30 = 0;
  int numb  = 0;
  if(jec==1){
    for(unsigned int n = 0; n<jets_up_csv().size();++n){
      if(jets_up_p4()[n].Pt()>=30&&fabs(jets_up_p4()[n].Eta())<=5)                             ++num;
      if(jets_up_p4()[n].Pt()>=30&&fabs(jets_up_p4()[n].Eta())<=2.5)                           ++num30;
      if(jets_up_p4()[n].Pt()>=20&&fabs(jets_up_p4()[n].Eta())<=2.4&&jets_up_csv()[n]>=0.5426) ++numb;
    }
  }
  else if(jec==(-1)){
    for(unsigned int n = 0; n<jets_dn_csv().size();++n){
      if(jets_dn_p4()[n].Pt()>=30&&fabs(jets_dn_p4()[n].Eta())<=5)                             ++num;
      if(jets_dn_p4()[n].Pt()>=30&&fabs(jets_dn_p4()[n].Eta())<=2.5)                           ++num30;
      if(jets_dn_p4()[n].Pt()>=20&&fabs(jets_dn_p4()[n].Eta())<=2.4&&jets_dn_csv()[n]>=0.5426) ++numb;
    }
  }
  else{ 
    for(unsigned int n = 0; n<jets_csv().size();++n){
      if(   jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=5)                             ++num;
      if(   jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=2.5)                           ++num30;
      if(   jets_p4()[n].Pt()>=20&&fabs(   jets_p4()[n].Eta())<=2.4&&   jets_csv()[n]>=0.5426) ++numb;
    }
  }
  nj   = num;
  nj30 = num30;
  nb   = numb;
  return true;
}

bool istightlepton(int i, int version){
  bool istight = false;
  if(version==0) {
    if(lep_pass_VVV_cutbased_tight_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015&&lep_isTriggerSafe_v1()[i]) {
      bool isgoodelectron = abs(lep_pdgId()[i])==11&&lep_lostHits()[i]==0&&lep_tightCharge()[i]==2 &&( (fabs(lep_etaSC()[i])<=1.479&&lep_relIso03EAv2()[i]<0.0588) || (fabs(lep_etaSC()[i]) >1.479&&lep_relIso03EAv2()[i]<0.0571) );
      istight = abs(lep_pdgId()[i])==13 &&lep_relIso03EAv2()[i]<0.06 ? true : isgoodelectron;
    }
  }
  if(version==1) {
    bool isgoodelectron = abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&lep_isTriggerSafe_v1()[i]&&fabs(lep_ip3d()[i])<0.01&&lep_tightCharge()[i]==2&&( (fabs(lep_etaSC()[i])<=1.479&&lep_ptRatio()[i]>0.9&&lep_MVA()[i]>0.941) || (fabs(lep_etaSC()[i]) >1.479&&lep_ptRatio()[i]>0.9&&lep_MVA()[i]>0.925) );
    bool isgoodmuon     = abs(lep_pdgId()[i])==13&&lep_pass_VVV_cutbased_tight_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015&&lep_tightCharge()[i]==2&&lep_ptRatio()[i]>0.9;
    istight = abs(lep_pdgId()[i])==13 ? isgoodmuon : isgoodelectron;
  }
  if(version==-1) {
    bool isgoodelectron = abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&lep_isTriggerSafe_v1()[i]&&fabs(lep_ip3d()[i])<0.015&&( (fabs(lep_etaSC()[i])<=1.479&&lep_ptRatio()[i]>0.84&&lep_MVA()[i]>0.920) || (fabs(lep_etaSC()[i]) >1.479&&lep_ptRatio()[i]>0.84&&lep_MVA()[i]>0.880) );
    bool isgoodmuon     = abs(lep_pdgId()[i])==13&&lep_pass_VVV_cutbased_tight_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015&&lep_ptRatio()[i]>0.84;
    istight = abs(lep_pdgId()[i])==13 ? isgoodmuon : isgoodelectron;
  }
  return istight;
}

bool islooselepton(int i, int version){
   bool isloose = false;
   if(version ==0){
     isloose = abs(lep_pdgId()[i])==13 ? lep_pass_VVV_cutbased_fo_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015&&lep_isTriggerSafe_v1()[i] && lep_relIso03EAv2()[i]<0.4 :  lep_pass_VVV_cutbased_fo_noiso()[i]   &&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015&&lep_isTriggerSafe_v1()[i] && lep_lostHits()[i]==0&&lep_tightCharge()[i]==2 && lep_relIso03EAv2()[i]<0.2;
   }   
   if(version ==1){
     bool islooseelectron = abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&lep_isTriggerSafe_v1()[i]&&fabs(lep_ip3d()[i])<0.01&&lep_tightCharge()[i]==2&&( (fabs(lep_etaSC()[i])<=1.479&&lep_ptRatio()[i]>0.65&&lep_MVA()[i]>0.941) || (fabs(lep_etaSC()[i]) >1.479&&lep_ptRatio()[i]>0.65&&lep_MVA()[i]>0.925) );
     bool isloosemuon     = abs(lep_pdgId()[i])==13&&lep_pass_VVV_cutbased_tight_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015&&lep_tightCharge()[i]==2&&lep_ptRatio()[i]>0.65;
     isloose = abs(lep_pdgId()[i])==13 ? isloosemuon : islooseelectron;
   }   
   if(version ==-1){
     bool islooseelectron = abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&lep_isTriggerSafe_v1()[i]&&fabs(lep_ip3d()[i])<0.015&&( (fabs(lep_etaSC()[i])<=1.479&&lep_relIso04EAv2()[i]<0.4&&lep_MVA()[i]>0.920) || (fabs(lep_etaSC()[i]) >1.479&&lep_ptRatio()[i]>0.65&&lep_MVA()[i]>0.880) );
     bool isloosemuon     = abs(lep_pdgId()[i])==13&&lep_pass_VVV_cutbased_tight_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015&&lep_tightCharge()[i]==2&&lep_ptRatio()[i]>0.65;
     isloose = abs(lep_pdgId()[i])==13 ? isloosemuon : islooseelectron;
   }   
   return isloose;
}

bool getleptonindices_v0(vector<int> &iSS, vector<int> &i3l, vector<int> &iaSS, vector<int> &ia3l, vector<int> &vSS, vector<int> &v3l, vector<int> &vaSS, vector<int> &va3l, int version, float pTSS, float pT3l){
  iSS.clear();
  i3l.clear();
  iaSS.clear();
  ia3l.clear();
  vSS.clear();
  v3l.clear();
  vaSS.clear();
  va3l.clear();
  for(unsigned int i = 0; i<lep_pdgId().size();++i){
    bool isSS    = false; bool is3l    = false;
    bool isaSS   = false; bool isa3l   = false;
    //tight ID
    if(istightlepton(i,-1*version)&&lep_p4()[i].Pt()>pT3l) { i3l.push_back(i); is3l = true; }
    if(istightlepton(i,   version)&&lep_p4()[i].Pt()>pTSS) { iSS.push_back(i); isSS = true; }
   //loose ID
    if(islooselepton(i,-1*version)&&!is3l&&lep_p4()[i].Pt()>pT3l) { ia3l.push_back(i); isa3l = true; }
    if(islooselepton(i,   version)&&!isSS&&lep_p4()[i].Pt()>pTSS) { iaSS.push_back(i); isaSS = true; }
    //vetoID --simple, leave it here.
    if (version == 0)
    {
        if(lep_pass_VVV_cutbased_veto_noiso()[i] &&fabs(lep_p4()[i].Eta())<2.4&&lep_relIso03EAv2()[i]<=0.4) {
            if(!isSS &&          lep_p4()[i].Pt()>10) vSS  .push_back(i);
            if(!is3l &&          lep_p4()[i].Pt()>10) v3l  .push_back(i);
            if(!isSS && !isaSS &&lep_p4()[i].Pt()>10) vaSS .push_back(i);
            if(!is3l && !isa3l &&lep_p4()[i].Pt()>10) va3l .push_back(i);
        }
    }
    else if (version == 1)
    {
        if (abs(lep_pdgId()[i]) == 11)
        {
            bool isvetoel = true;
            if (fabs(lep_etaSC().at(i)) > 2.5) isvetoel = false;
            if (fabs(lep_dxy().at(i)) > 0.05)  isvetoel = false;
            if (fabs(lep_dz().at(i)) > 0.1)    isvetoel = false;
            if (lep_ptRatio()[i]<=0.58)        isvetoel = false;
	    if(fabs(lep_etaSC()[i])<=1.479){
	      if (lep_MVA()[i]<=0.6)           isvetoel = false;
	    } else {
	      if (lep_MVA()[i]<=0.0)           isvetoel = false;
	    }
            if(isvetoel) {
                if(!isSS &&          lep_p4()[i].Pt()>10) vSS  .push_back(i);
                if(!is3l &&          lep_p4()[i].Pt()>10) v3l  .push_back(i);
                if(!isSS && !isaSS &&lep_p4()[i].Pt()>10) vaSS .push_back(i);
                if(!is3l && !isa3l &&lep_p4()[i].Pt()>10) va3l .push_back(i);
            }
        }
        else if (abs(lep_pdgId()[i]) == 13)
        {
            if(lep_pass_VVV_cutbased_veto_noiso()[i] &&fabs(lep_p4()[i].Eta())<2.4&&lep_ptRatio()[i]>0.58) {
                if(!isSS &&          lep_p4()[i].Pt()>10) vSS  .push_back(i);
                if(!is3l &&          lep_p4()[i].Pt()>10) v3l  .push_back(i);
                if(!isSS && !isaSS &&lep_p4()[i].Pt()>10) vaSS .push_back(i);
                if(!is3l && !isa3l &&lep_p4()[i].Pt()>10) va3l .push_back(i);
            }
        }
    }
  }
  return true;
}

#ifdef USE_CMS3_WWW100
bool getleptonindices_v2(vector<int> &iSS, vector<int> &i3l, vector<int> &iaSS, vector<int> &ia3l, vector<int> &vSS, vector<int> &v3l, vector<int> &vaSS, vector<int> &va3l, float pTSS, float pT3l){
  iSS.clear();
  i3l.clear();
  iaSS.clear();
  ia3l.clear();
  vSS.clear();
  v3l.clear();
  vaSS.clear();
  va3l.clear();
  for(unsigned int i = 0; i<lep_pdgId().size();++i){
    bool isSS    = false; bool is3l    = false;
    bool isaSS   = false; bool isa3l   = false;
    //tight ID
    if(lep_pass_VVV_cutbased_3l_tight()[i]&&lep_p4()[i].Pt()>pT3l) { i3l.push_back(i); is3l = true; }
    if(lep_pass_VVV_cutbased_tight   ()[i]&&lep_p4()[i].Pt()>pTSS) { iSS.push_back(i); isSS = true; }
    //loose ID
    if(lep_pass_VVV_cutbased_3l_fo()[i]&&!is3l&&lep_p4()[i].Pt()>pT3l) { ia3l.push_back(i); isa3l = true; }
    if(lep_pass_VVV_cutbased_fo   ()[i]&&!isSS&&lep_p4()[i].Pt()>pTSS) { iaSS.push_back(i); isaSS = true; }
    //veto ID
    if(lep_pass_VVV_cutbased_veto()[i] && !isSS &&          lep_p4()[i].Pt()>10) vSS  .push_back(i);
    if(lep_pass_VVV_cutbased_veto()[i] && !is3l &&          lep_p4()[i].Pt()>10) v3l  .push_back(i);
    if(lep_pass_VVV_cutbased_veto()[i] && !isSS && !isaSS &&lep_p4()[i].Pt()>10) vaSS .push_back(i);
    if(lep_pass_VVV_cutbased_veto()[i] && !is3l && !isa3l &&lep_p4()[i].Pt()>10) va3l .push_back(i);
  }
  return true;
}

int NtightSS(float pTSS){
  int n = 0;
  for(unsigned int i = 0; i<lep_pdgId().size();++i){
    if(lep_pass_VVV_cutbased_tight   ()[i]&&lep_p4()[i].Pt()>pTSS) ++n;
  }
  return n;
}
int NlooseSS(float pTSS){
  int n = 0;
  for(unsigned int i = 0; i<lep_pdgId().size();++i){
    if(lep_pass_VVV_cutbased_fo      ()[i]&&lep_p4()[i].Pt()>pTSS) ++n;
  }
  return n;
}
int Ntight3l(float pT3l){
  int n = 0;
  for(unsigned int i = 0; i<lep_pdgId().size();++i){
    if(lep_pass_VVV_cutbased_3l_tight()[i]&&lep_p4()[i].Pt()>pT3l) ++n;
  }
  return n;
}
int Nloose3l(float pT3l){
  int n = 0;
  for(unsigned int i = 0; i<lep_pdgId().size();++i){
    if(lep_pass_VVV_cutbased_3l_fo   ()[i]&&lep_p4()[i].Pt()>pT3l) ++n;
  }
  return n;
}

#endif

float coneCorrPt(int lepi, int version){
  float coneptcorr = 0;
  float relIso = lep_relIso03EAv2().at(lepi);
  float relIso04 = lep_relIso04EAv2().at(lepi);
  float ptratio = lep_ptRatio().at(lepi);
  if (version == 0) {
      if(   abs(lep_pdgId().at(lepi)) == 11) {
          if(fabs(lep_p4().at(lepi).Eta()) <= 1.479) coneptcorr = std::max( 0., relIso - 0.0588 );
          else                                       coneptcorr = std::max( 0., relIso - 0.0571 );
      }
      if(   abs(lep_pdgId().at(lepi)) == 13)         coneptcorr = std::max( 0., relIso - 0.06   );
  }
  else if (version == 1) {
      if(   abs(lep_pdgId().at(lepi)) == 11) {
          if(fabs(lep_p4().at(lepi).Eta()) <= 1.479) coneptcorr = std::max( 0., relIso04 - 0.05 );
          else                                       coneptcorr = std::max( 0., relIso04 - 0.07 );
      }
      if(   abs(lep_pdgId().at(lepi)) == 13)         coneptcorr = std::max( 0., 0.9/ptratio - 1.);
  }
  else if (version == -1){
      if(   abs(lep_pdgId().at(lepi)) == 11)         coneptcorr = std::max( 0., relIso04 - 0.1   );
      if(   abs(lep_pdgId().at(lepi)) == 13)         coneptcorr = std::max( 0., 0.84/ptratio - 1.);
  }

  return coneptcorr;
}

float calcMjj(bool closestDR, int jec){
  float minDR=999.;
  int jDR1(-1), jDR2(-1);
  if(jec==1){
    for(unsigned int j1 = 0; j1<jets_up_p4().size();++j1){
      if(fabs(jets_up_p4()[j1].Eta())>2.5) continue;
      if(     jets_up_p4()[j1].Pt()  <30.) continue;
      for(unsigned int j2 = j1+1; j2<jets_up_p4().size();++j2){
	if(fabs(jets_up_p4()[j2].Eta())>2.5) continue;
	if(     jets_up_p4()[j2].Pt()  <30.) continue;
	if(!closestDR) return  (jets_up_p4()[j1]+jets_up_p4()[j2]).M();
	if(       dR(jets_up_p4()[j1], jets_up_p4()[j2])<minDR){
	  minDR = dR(jets_up_p4()[j1], jets_up_p4()[j2]);
	  jDR1 = j1; jDR2 = j2;
	}
      }
    }
    if(jDR1>=0&&jDR2>=2) return (jets_up_p4()[jDR1]+jets_up_p4()[jDR2]).M();
  }
  else if(jec==(-1)){
    for(unsigned int j1 = 0; j1<jets_dn_p4().size();++j1){
      if(fabs(jets_dn_p4()[j1].Eta())>2.5) continue;
      if(     jets_dn_p4()[j1].Pt()  <30.) continue;
      for(unsigned int j2 = j1+1; j2<jets_dn_p4().size();++j2){
	if(fabs(jets_dn_p4()[j2].Eta())>2.5) continue;
	if(     jets_dn_p4()[j2].Pt()  <30.) continue;
	if(!closestDR) return  (jets_dn_p4()[j1]+jets_dn_p4()[j2]).M();
	if(dR(       jets_dn_p4()[j1], jets_dn_p4()[j2])<minDR){
	  minDR = dR(jets_dn_p4()[j1], jets_dn_p4()[j2]);
	  jDR1 = j1; jDR2 = j2;
	}
      }
    }
    if(jDR1>=0&&jDR2>=2) return (jets_dn_p4()[jDR1]+jets_dn_p4()[jDR2]).M();
  }
  else {
    for(unsigned int j1 = 0; j1<jets_p4().size();++j1){
      if(fabs(jets_p4()[j1].Eta())>2.5) continue;
      if(     jets_p4()[j1].Pt()  <30.) continue;
      for(unsigned int j2 = j1+1; j2<jets_p4().size();++j2){
	if(fabs(jets_p4()[j2].Eta())>2.5) continue;
	if(     jets_p4()[j2].Pt()  <30.) continue;
	if(!closestDR) return  (jets_p4()[j1]+jets_p4()[j2]).M();
	if(dR(       jets_p4()[j1], jets_p4()[j2])<minDR){
	  minDR = dR(jets_p4()[j1], jets_p4()[j2]);
	  jDR1 = j1; jDR2 = j2;
	}
      }
    }
    if(jDR1>=0&&jDR2>=2) return (jets_p4()[jDR1]+jets_p4()[jDR2]).M();
  }
  return -1;
}

float Detajj(int jec){
  if(jec==1){
    for(unsigned int j1 = 0; j1<jets_up_p4().size();++j1){
      if(fabs(jets_up_p4()[j1].Eta())>2.5) continue;
      if(     jets_up_p4()[j1].Pt()  <30.) continue;
      for(unsigned int j2 = j1+1; j2<jets_up_p4().size();++j2){
	if(fabs(jets_up_p4()[j2].Eta())>2.5) continue;
	if(     jets_up_p4()[j2].Pt()  <30.) continue;
	return dEta(jets_up_p4()[j1],jets_up_p4()[j2]);
      }
    }
  }
  else if(jec==(-1)){
    for(unsigned int j1 = 0; j1<jets_dn_p4().size();++j1){
      if(fabs(jets_dn_p4()[j1].Eta())>2.5) continue;
      if(     jets_dn_p4()[j1].Pt()  <30.) continue;
      for(unsigned int j2 = j1+1; j2<jets_dn_p4().size();++j2){
	if(fabs(jets_dn_p4()[j2].Eta())>2.5) continue;
	if(     jets_dn_p4()[j2].Pt()  <30.) continue;
	return dEta(jets_dn_p4()[j1],jets_dn_p4()[j2]);
      }
    }
  }
  else {
    for(unsigned int j1 = 0; j1<jets_p4().size();++j1){
      if(fabs(jets_p4()[j1].Eta())>2.5) continue;
      if(     jets_p4()[j1].Pt()  <30.) continue;
      for(unsigned int j2 = j1+1; j2<jets_p4().size();++j2){
	if(fabs(jets_p4()[j2].Eta())>2.5) continue;
	if(     jets_p4()[j2].Pt()  <30.) continue;
	return dEta(jets_p4()[j1],jets_p4()[j2]);
      }
    }
  }
  return -1;
}

bool getMjjAndDeta(float &Mjj, float &MjjL, float &Detajj, int jec){//xxx
  float Massjj     = -1;
  float MassjjL    = -1;
  float Deltaetajj = -1;
  float minDR=999.;
  int jDR1(-1), jDR2(-1);
  if(jec==1){
    for(unsigned int j1 = 0; j1<jets_up_p4().size();++j1){
      if(fabs(jets_up_p4()[j1].Eta())>2.5) continue;
      if(     jets_up_p4()[j1].Pt()  <30.) continue;
      for(unsigned int j2 = j1+1; j2<jets_up_p4().size();++j2){
	if(fabs(jets_up_p4()[j2].Eta())>2.5) continue;
	if(     jets_up_p4()[j2].Pt()  <30.) continue;
	if(MassjjL<0) {
	  MassjjL    =     (jets_up_p4()[j1]+jets_up_p4()[j2]).M();
	  Deltaetajj = dEta(jets_up_p4()[j1],jets_up_p4()[j2]);
	}
	if(       dR(jets_up_p4()[j1], jets_up_p4()[j2])<minDR){
	  minDR = dR(jets_up_p4()[j1], jets_up_p4()[j2]);
	  jDR1 = j1; jDR2 = j2;
	}
      }
    }
    if(jDR1>=0&&jDR2>=0) Massjj = (jets_up_p4()[jDR1]+jets_up_p4()[jDR2]).M();
  }
  else if(jec==(-1)){
    for(unsigned int j1 = 0; j1<jets_dn_p4().size();++j1){
      if(fabs(jets_dn_p4()[j1].Eta())>2.5) continue;
      if(     jets_dn_p4()[j1].Pt()  <30.) continue;
      for(unsigned int j2 = j1+1; j2<jets_dn_p4().size();++j2){
	if(fabs(jets_dn_p4()[j2].Eta())>2.5) continue;
	if(     jets_dn_p4()[j2].Pt()  <30.) continue;
	if(MassjjL<0) {
	  MassjjL    =     (jets_dn_p4()[j1]+jets_dn_p4()[j2]).M();
	  Deltaetajj = dEta(jets_dn_p4()[j1],jets_dn_p4()[j2]);
	}
	if(       dR(jets_dn_p4()[j1], jets_dn_p4()[j2])<minDR){
	  minDR = dR(jets_dn_p4()[j1], jets_dn_p4()[j2]);
	  jDR1 = j1; jDR2 = j2;
	}
      }
    }
    if(jDR1>=0&&jDR2>=0) Massjj = (jets_dn_p4()[jDR1]+jets_dn_p4()[jDR2]).M();
  }
  else {
    for(unsigned int j1 = 0; j1<jets_p4().size();++j1){
      if(fabs(jets_p4()[j1].Eta())>2.5) continue;
      if(     jets_p4()[j1].Pt()  <30.) continue;
      for(unsigned int j2 = j1+1; j2<jets_p4().size();++j2){
	if(fabs(jets_p4()[j2].Eta())>2.5) continue;
	if(     jets_p4()[j2].Pt()  <30.) continue;
	if(MassjjL<0) {
	  MassjjL    =     (jets_p4()[j1]+jets_p4()[j2]).M();
	  Deltaetajj = dEta(jets_p4()[j1],jets_p4()[j2]);
	}
	if(       dR(jets_p4()[j1], jets_p4()[j2])<minDR){
	  minDR = dR(jets_p4()[j1], jets_p4()[j2]);
	  jDR1 = j1; jDR2 = j2;
	}
      }
    }
    if(jDR1>=0&&jDR2>=0) Massjj = (jets_p4()[jDR1]+jets_p4()[jDR2]).M();
  }
  Mjj    = Massjj;
  MjjL   = MassjjL;
  Detajj = Deltaetajj;
  return true;
}

float calcMTmax(vector<int> lepindex, LorentzVector MET, bool compareSSpairs){
  float MT = -1;
  if(!compareSSpairs){
    for(unsigned int n = 0; n<lepindex.size(); ++n){
      float myMT = mT(lep_p4()[lepindex[n] ],MET);
      if(MT<myMT) MT = myMT;
    }
    return MT;
  }
  float MT1(-1),MT2(-1);
  int ns1(0), ns2(0);
  for(unsigned int n = 0; n<lepindex.size(); ++n){
    float myMT = mT(lep_p4()[lepindex[n] ],MET);
    if((lep_pdgId()[lepindex[n] ])>0){
      if(MT1<myMT) MT1 = myMT;
      ++ns1;
    }
    else {
      if(MT2<myMT) MT2 = myMT;
      ++ns2;
    }
  }
  if(ns1>ns2) return MT1;
  else        return MT2;

}
float calcMTmax(int index1, int index2, LorentzVector MET){
  if(index1==index2)               return -1.;
  if(index1<0)                     return -1.;
  if(index1>=(int)lep_p4().size()) return -1.;
  if(index2>=(int)lep_p4().size()) return -1.;
  float MT1         = mT(lep_p4()[index1],MET);
  float MT2 = -1.;
  if(index2>0)  MT2 = mT(lep_p4()[index2],MET);
  if(MT1>MT2)  return MT1;
  return MT2;
}

bool passTriggers(bool MCoffline, bool MCtrigger, float leadleppt, float trailleppt){
  if(!MCoffline && !MCtrigger && !isData()) return true;//no offline requirements, and no trigger emulation
  int nel25 = 0;
  int nel = 0;
  int nmu25 = 0;
  int nmu = 0;
  for(unsigned int i = 0; i<lep_pdgId().size();++i){
    if(lep_pass_VVV_cutbased_3l_fo()[i]&&fabs(lep_p4()[i].Eta())<2.5&&lep_p4()[i].Pt()>leadleppt) {
      if(     abs(lep_pdgId()[i])==11) ++nel25;
      else if(abs(lep_pdgId()[i])==13) ++nmu25;
    }
    if(lep_pass_VVV_cutbased_3l_fo()[i]&&fabs(lep_p4()[i].Eta())<2.5&&lep_p4()[i].Pt()>trailleppt) {
      if(     abs(lep_pdgId()[i])==11) ++nel;
      else if(abs(lep_pdgId()[i])==13) ++nmu;
    }
  }
  if(MCoffline && !MCtrigger && !isData()){
    if(nmu>=2)           return true;
    if(nmu25>=1&&nel>=1) return true;
    if(nel25>=1&&nmu>=1) return true;
    if(nel25>=1&&nel>=2) return true;
    return false;
  }
  if(nmu>=2&&(HLT_DoubleMu()) )                             return true;
  if(nmu25>=1&&nel>=1&&HLT_MuEG())                          return true;
  if(nel25>=1&&nmu>=1&&HLT_MuEG())                          return true;
  if(nel25>=1&&nel>=2&&(HLT_DoubleEl()||HLT_DoubleEl_DZ())) return true;
  return false;
}

bool passofflineTriggers(vector<int> tightlep, vector<int> looselep){
  bool passofflineforTrigger = false;
  int nel25 = 0;
  int nel = 0;
  int nmu25 = 0;
  int nmu = 0;
  for(unsigned int i = 0; i<tightlep.size(); ++i){
    if(abs(lep_pdgId()[tightlep[i] ])==11){
      ++nel;
      if(lep_p4()[tightlep[i] ].Pt()>25) ++nel25;
    } else if(abs(lep_pdgId()[tightlep[i] ])==13){
      ++nmu;
      if(lep_p4()[tightlep[i] ].Pt()>25) ++nmu25;
    }
  }
  for(unsigned int i = 0; i<looselep.size(); ++i){
    if(abs(lep_pdgId()[looselep[i] ])==11){
      ++nel;
      if(lep_p4()[looselep[i] ].Pt()>25) ++nel25;
    } else if(abs(lep_pdgId()[looselep[i] ])==13){
      ++nmu;
      if(lep_p4()[looselep[i] ].Pt()>25) ++nmu25;
    }
  }
  if(nmu>=2)           return true;
  if(nmu25>=1&&nel>=1) return true;
  if(nel25>=1&&nmu>=1) return true;
  if(nel25>=1&&nel>=2) return true;
  return false;
}

bool passonlineTriggers(vector<int> tightlep, vector<int> looselep){
  if(!isData()) return true;//no trigger emulation
  int nel25 = 0;
  int nel = 0;
  int nmu25 = 0;
  int nmu = 0;
  for(unsigned int i = 0; i<tightlep.size(); ++i){
    if(abs(lep_pdgId()[tightlep[i] ])==11){
      ++nel;
      if(lep_p4()[tightlep[i] ].Pt()>25) ++nel25;
    } else if(abs(lep_pdgId()[tightlep[i] ])==13){
      ++nmu;
      if(lep_p4()[tightlep[i] ].Pt()>25) ++nmu25;
    }
  }
  for(unsigned int i = 0; i<looselep.size(); ++i){
    if(abs(lep_pdgId()[looselep[i] ])==11){
      ++nel;
      if(lep_p4()[looselep[i] ].Pt()>25) ++nel25;
    } else if(abs(lep_pdgId()[looselep[i] ])==13){
      ++nmu;
      if(lep_p4()[looselep[i] ].Pt()>25) ++nmu25;
    }
  }
  if(nmu>=2&&(HLT_DoubleMu()) )                             return true;
  if(nmu25>=1&&nel>=1&&HLT_MuEG())                          return true;
  if(nel25>=1&&nmu>=1&&HLT_MuEG())                          return true;
  if(nel25>=1&&nel>=2&&(HLT_DoubleEl()||HLT_DoubleEl_DZ())) return true;
  return false;
}

bool passFilters_v0(){
  if(!isData()) return true;//don't apply filters on simulation
  if (!Flag_EcalDeadCellTriggerPrimitiveFilter() ) return false;
  if (!Flag_badChargedCandidateFilter         () ) return false;
  if (!Flag_HBHENoiseFilter                   () ) return false;
  if (!Flag_HBHEIsoNoiseFilter                () ) return false;
  if (!Flag_goodVertices                      () ) return false;
  if (!Flag_eeBadScFilter                     () ) return false;
  if (!Flag_globalTightHalo2016               () ) return false;
  if (!Flag_badMuonFilter                     () ) return false;
  return true;
}

bool passFilters(){
  if(!isData()) return true;//don't apply filters on simulation
  if (!Flag_AllEventFilters() ) return false;
  return true;
}

bool splitVH(string filename){
  if(filename.find("vh_nonbb_amcnlo")==string::npos) return false;//file is certainly no WHtoWWW
  bool isHtoWW = false;
  bool isWnotFromH = false;
  for(unsigned int i = 0; i<genPart_pdgId().size();++i){
    if(abs(genPart_pdgId()[i])!=24) continue;
    if(    genPart_status()[i]!=22) continue;
    if(abs(genPart_motherId()[i])==25) { isHtoWW     = true; }
    if(abs(genPart_motherId()[i])!=2)  { isWnotFromH = true; }
    if(isHtoWW&&isWnotFromH) break;
  }
  return isHtoWW&&isWnotFromH;
}

string process(string filename, bool SS, vector<int> tightlep, vector<int> looselep){
  if(isData())                                        return "Data";
  if(splitVH(filename))                               return "WHtoWWW";
  if(filename.find("www_2l_")         !=string::npos) return "WWW";
  if(filename.find("www_incl_amcnlo_")!=string::npos) return "WWWv2";
  if(filename.find("data_")           !=string::npos) return "Data"; 
  if(SS){
    if((tightlep.size()+looselep.size())<2)           return "not2l";
    int l1, l2;
    if(tightlep.size()>=2)      { l1 = tightlep[0]; l2 = tightlep[1]; }
    else if(tightlep.size()==1) { l1 = tightlep[0]; l2 = looselep[0]; }
    else                        { l1 = looselep[0]; l2 = looselep[1]; }
    int gentype = gentype_v2(l1,l2,-1);
    if(     gentype==0) return "trueSS";
    else if(gentype==2) return "chargeflips";
    else if(gentype==3) return "SSLL";
    else if(gentype==4) return "fakes";
    else if(gentype==5) return "photonfakes";
    else                return "others";
  }
  //this is 3l
  if((tightlep.size()+looselep.size())<3) return "not3l";
  int l1, l2, l3;
  if(tightlep.size()>=3)      { l1 = tightlep[0]; l2 = tightlep[1]; l3 = tightlep[2]; }
  else if(tightlep.size()==2) { l1 = tightlep[0]; l2 = tightlep[1]; l3 = looselep[0]; }
  else if(tightlep.size()==1) { l1 = tightlep[0]; l2 = looselep[0]; l3 = looselep[1]; }
  else                        { l1 = looselep[0]; l2 = looselep[1]; l3 = looselep[2]; }
  int gentype = gentype_v2(l1,l2,l3);
  if(     gentype==0) return "trueWWW";
  else if(gentype==1) return "true3L";
  else if(gentype==2) return "chargeflips";
  else if(gentype==3) return "3lLL";
  else if(gentype==4) return "fakes";
  else if(gentype==5) return "photonfakes";
  else                return "others";
}

bool vetophotonprocess(string filename, bool process){
  if((filename.find("wjets_" )!=string::npos||filename.find("dy_"    )!=string::npos||filename.find("ttbar_" )!=string::npos||filename.find("ww_" )!=string::npos||filename.find("wz_" )!=string::npos)&&(process)) return true;
  if((filename.find("wgjets_")!=string::npos||filename.find("wgstar_")!=string::npos||filename.find("zgamma_")!=string::npos||filename.find("ttg_")!=string::npos||filename.find("wwg_")!=string::npos ||filename.find("wzg_")!=string::npos)&&(!process)) return true;
  return false;
}

bool passJetSSstate_v0(bool preselect, int nj, int nb, float Mjj, float MjjL, float Detajj, bool is3lCR, int jec, bool btag, bool Mjjside, int version){
  int numj = nj;
  int numb = nb;
  int x;
  if(nb<0||nj<0) getalljetnumbers(x,numj,numb, jec);
  if(numj<2)         return false;
  if( btag&&numb==0) return false;
  if(!btag&&numb>0)  return false;
  if(preselect)      return true;
  float Massjj  = Mjj;
  float MassjjL = MjjL;
  float Deta    = Detajj;
  if(Mjj<0||MjjL<0||Detajj<0) getMjjAndDeta(Massjj,MassjjL,Deta,jec);
  if(fabs(Deta)>1.5)                         return false;
  if(fabs(MassjjL)>400.)                     return false;
  if(is3lCR)                                 return true;
  float masswindow = 20.;
  if(abs(version)==1) masswindow = 15.;
  if( Mjjside&&fabs(Massjj-80.) <masswindow) return false;
  if(!Mjjside&&fabs(Massjj-80.)>=masswindow) return false;
  return true;
}

bool passJetSSstate(bool preselect, bool is3lCR, int jec, bool btag, bool Mjjside, int version){
    float masswindow = 20.;
    if(abs(version)==1) masswindow = 15.;
    if(jec==1){
      if(nj30_up()<2)                              return false;
      if( btag&&nb_up()==0)                        return false;
      if(!btag&&nb_up()>0)                         return false;
      if(preselect)                              return true;
      if(fabs(DetajjL_up())>1.5)                   return false;
      if(fabs(MjjL_up())>400.)                     return false;
      if(is3lCR)                                 return true;
      if( Mjjside&&fabs(Mjj_up()-80.) <masswindow) return false;
      if(!Mjjside&&fabs(Mjj_up()-80.)>=masswindow) return false;
      return true;
    } else if(jec==(-1)){
      if(nj30_dn()<2)                              return false;
      if( btag&&nb_dn()==0)                        return false;
      if(!btag&&nb_dn()>0)                         return false;
      if(preselect)                              return true;
      if(fabs(DetajjL_dn())>1.5)                   return false;
      if(fabs(MjjL_dn())>400.)                     return false;
      if(is3lCR)                                 return true;
      if( Mjjside&&fabs(Mjj_dn()-80.) <masswindow) return false;
      if(!Mjjside&&fabs(Mjj_dn()-80.)>=masswindow) return false;
      return true;
    } else {
      if(nj30()<2)                              return false;
      if( btag&&nb()==0)                        return false;
      if(!btag&&nb()>0)                         return false;
      if(preselect)                           return true;
      if(fabs(DetajjL())>1.5)                   return false;
      if(fabs(MjjL())>400.)                     return false;
      if(is3lCR)                              return true;
      if( Mjjside&&fabs(Mjj()-80.) <masswindow) return false;
      if(!Mjjside&&fabs(Mjj()-80.)>=masswindow) return false;
      return true;
    }
  return false;//should not happen
}

int isSRSS(bool preselect, int jec, bool btag, bool Mjjside, int version){
  if(!passJetSSstate(preselect, false, jec, btag, Mjjside,version)) return -1;
  float met = met_pt();
  if(     jec==  1 ) met = met_up_pt();
  else if(jec==(-1)) met = met_dn_pt();
  if(NtightSS()!=2)                                return -1;
  if(nVlep()!=2)                                   return -1;
  if(nisoTrack_mt2_cleaned_VVV_cutbased_veto()!=0) return -1;
  if(preselect){
    if(passSSee() && fabs(MeeSS()-MZ)>10.)         return  0;
    if(passSSem())                                 return  1;
    if(passSSmm())                                 return  2;
    return -1;
  }
  float metcut = 40.;
  if(abs(version)==1) metcut = 60.;
  float maxMT = MTmax();
  if(jec==  1)  maxMT = MTmax_up();
  if(jec==(-1)) maxMT = MTmax_dn();
  if(passSSee() && fabs(MeeSS()-MZ)>10.&&met>metcut&&MllSS()>40.)           return 0;
  if(passSSem() &&                       met>metcut&&MllSS()>30.&&maxMT>90) return 1;
  if(passSSmm() &&                                   MllSS()>40.)           return 2;
  return -1;
}
int isARSS(bool preselect, int jec, bool btag, bool Mjjside, int version){
  if(!passJetSSstate(preselect, false, jec, btag, Mjjside,version)) return -1;
  float met = met_pt();
  if(     jec==  1 ) met = met_up_pt();
  else if(jec==(-1)) met = met_dn_pt();
  if(NtightSS()!=1)                                return -1;
  if(NlooseSS()!=2)                                return -1;
  if(nVlep()!=2)                                   return -1;
  if(nisoTrack_mt2_cleaned_VVV_cutbased_veto()!=0) return -1;
  if(preselect){
    if(passSSee() && fabs(MeeSS()-MZ)>10.)         return  0;
    if(passSSem())                                 return  1;
    if(passSSmm())                                 return  2;
    return -1;
  }
  float metcut = 40.;
  if(abs(version)==1) metcut = 60.;
  float maxMT = MTmax();
  if(jec==  1)  maxMT = MTmax_up();
  if(jec==(-1)) maxMT = MTmax_dn();
  if(passSSee() && fabs(MeeSS()-MZ)>10.&&met>metcut&&MllSS()>40.)           return 0;
  if(passSSem() &&                       met>metcut&&MllSS()>30.&&maxMT>90) return 1;
  if(passSSmm() &&                                   MllSS()>40.)           return 2;
  return -1;
}
int isCRSS(bool preselect, int jec, bool noZ, bool btag, bool Mjjside, int version){
  if(!passJetSSstate(preselect, true, jec, btag, Mjjside,version)) return -1;
  float met = met_pt();
  if(     jec==  1 ) met = met_up_pt();
  else if(jec==(-1)) met = met_dn_pt();
  if(NtightSS()<2)                                 return -1;
  if(Ntight3l()!=3)                                return -1;
  if(nVlep()!=3)                                   return -1;
  if(nisoTrack_mt2_cleaned_VVV_cutbased_veto()!=0) return -1;
  int lep3 = -1; int SS1 = lep_idx0_SS(); int SS2 = lep_idx1_SS();
  if(SS1<0||SS2<0) return -1;
  if(SS1==0&&SS2==1) lep3 = 2;
  else if(SS1==0&&SS2==2) lep3 = 1;
  else if(SS1==1&&SS2==2) lep3 = 0;
  if(!noZ){
    bool hasSFOSZ = false;
    if((lep_pdgId()[SS1]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[SS1]+lep_p4()[lep3]).M()-MZ)<=10)) hasSFOSZ = true;
    if((lep_pdgId()[SS2]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[SS2]+lep_p4()[lep3]).M()-MZ)<=10)) hasSFOSZ = true;
    if(!hasSFOSZ) return -1;
  } else {
    bool hasSFOS = false;
    if((lep_pdgId()[SS1]==(-lep_pdgId()[lep3]))) hasSFOS = true;
    if((lep_pdgId()[SS2]==(-lep_pdgId()[lep3]))) hasSFOS = true;
   if(!hasSFOS)  return -1;
  }
  if(preselect){
    if(passSSee() && fabs(MeeSS()-MZ)>10.)         return  0;
    if(passSSem())                                 return  1;
    if(passSSmm())                                 return  2;
    return -1;
  }
  float metcut = 40.;
  if(abs(version)==1) metcut = 60.;
  float maxMT = MTmax();
  if(jec==  1)  maxMT = MTmax_up();
  if(jec==(-1)) maxMT = MTmax_dn();
  if(passSSee() && fabs(MeeSS()-MZ)>10.&&met>metcut&&MllSS()>40.)           return 0;
  if(passSSem() &&                       met>metcut&&MllSS()>30.&&maxMT>90) return 1;
  if(passSSmm() &&                                   MllSS()>40.)           return 2;
  return -1;
}
bool passAnySS(int &SR, int &AR, int &CR, bool preselect, int jec, bool noZ, bool btag, bool Mjjside, int version){
  SR = -1; AR = -1; CR = -1;
  bool changeSR = (NtightSS()==2                  && nVlep()==2);  
  bool changeAR = (NtightSS()==1 && NlooseSS()==2 && nVlep()==2);  
  bool changeCR = (NtightSS()>=2 && Ntight3l()==3 && nVlep()==3);
  bool mychange = changeSR;
  if(nisoTrack_mt2_cleaned_VVV_cutbased_veto()!=0) return true;
  if(!passJetSSstate(preselect, true, jec, btag, Mjjside,version)) return true;//don't pass any selection
  if(!passJetSSstate(preselect, false, jec, btag, Mjjside,version)) { changeSR = false; changeAR = false; }
  float met = met_pt();
  if(     jec==  1 ) met = met_up_pt();
  else if(jec==(-1)) met = met_dn_pt();
  if(changeCR){
    int lep3 = -1; int SS1 = lep_idx0_SS(); int SS2 = lep_idx1_SS();
    if(SS1<0||SS2<0) return -1;
    if(SS1==0&&SS2==1) lep3 = 2;
    else if(SS1==0&&SS2==2) lep3 = 1;
    else if(SS1==1&&SS2==2) lep3 = 0;
    else { cout << "This should not happen SS1 " << SS1 << " SS2 " << SS2 << endl; return -1; }
    if(!noZ){
      bool hasSFOSZ = false;
      if((lep_pdgId()[SS1]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[SS1]+lep_p4()[lep3]).M()-MZ)<=10)) hasSFOSZ = true;
      if((lep_pdgId()[SS2]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[SS2]+lep_p4()[lep3]).M()-MZ)<=10)) hasSFOSZ = true;
      if(!hasSFOSZ) return true;
    } else {
      bool hasSFOS = false;
      if((lep_pdgId()[SS1]==(-lep_pdgId()[lep3]))) hasSFOS = true;
      if((lep_pdgId()[SS2]==(-lep_pdgId()[lep3]))) hasSFOS = true;
      if(!hasSFOS)  return true;
    }
  }
  if(preselect){
    if(passSSee() && fabs(MeeSS()-MZ)>10.) { SR = changeSR ? 0 : -1; AR = changeAR ? 0 : -1; CR = changeCR ? 0 : -1; }
    if(passSSem())                         { SR = changeSR ? 1 : -1; AR = changeAR ? 1 : -1; CR = changeCR ? 1 : -1; }
    if(passSSmm())                         { SR = changeSR ? 2 : -1; AR = changeAR ? 2 : -1; CR = changeCR ? 2 : -1; }
 
    return true;
  }
  float metcut = 40.;
  if(abs(version)==1) metcut = 60.;
  float maxMT = MTmax();
  if(jec==  1)  maxMT = MTmax_up();
  if(jec==(-1)) maxMT = MTmax_dn();
  if(passSSee() && fabs(MeeSS()-MZ)>10.&&met>metcut&&MllSS()>40.)           { SR = changeSR ? 0 : -1; AR = changeAR ? 0 : -1; CR = changeCR ? 0 : -1; }
  if(passSSem() &&                       met>metcut&&MllSS()>30.&&maxMT>90) { SR = changeSR ? 1 : -1; AR = changeAR ? 1 : -1; CR = changeCR ? 1 : -1; }
  if(passSSmm() &&                                   MllSS()>40.)           { SR = changeSR ? 2 : -1; AR = changeAR ? 2 : -1; CR = changeCR ? 2 : -1; }
  return true;
}

int isSRSS_v0(vector<int> tightlep, vector<int> vetolep,  bool preselect, float maxMT, int nj, int nb, float Mjj, float MjjL, float Detajj, LorentzVector MET, int jec, bool btag, bool Mjjside, int version){
  if(!passJetSSstate_v0(preselect, nj, nb, Mjj, MjjL, Detajj, false, jec, btag, Mjjside,version)) return -1;
  if(tightlep.size()!=2)                                                                          return -1;
  if(vetolep.size()!=0)                                                                           return -1;
  if(nisoTrack_mt2_cleaned_VVV_cutbased_veto()!=0)                                                return -1;
  if((lep_pdgId()[tightlep[0] ]*lep_pdgId()[tightlep[1] ])<0)                                     return -1;
  bool ee = false; bool em = false; bool mm = false;
LorentzVector met = LorentzVector(0,0,0,0);
  if(MET.Pt()>0)     met = MET;
  else if(jec==1)    met.SetPxPyPzE(met_up_pt()*std::cos(met_up_phi()),met_up_pt()*std::sin(met_up_phi()),0,met_up_pt());
  else if(jec==(-1)) met.SetPxPyPzE(met_dn_pt()*std::cos(met_dn_phi()),met_dn_pt()*std::sin(met_dn_phi()),0,met_dn_pt());
  else               met.SetPxPyPzE(met_pt()*std::cos(met_phi()),met_pt()*std::sin(met_phi()),0,met_pt());
 
  if(((lep_pdgId()[tightlep[0] ])*(lep_pdgId()[tightlep[1] ]))==121&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()-MZ)>10.) ee = true;
  if(((lep_pdgId()[tightlep[0] ])*(lep_pdgId()[tightlep[1] ]))==143)  em = true;
  if(((lep_pdgId()[tightlep[0] ])*(lep_pdgId()[tightlep[1] ]))==169)  mm = true;
  if(preselect){
    if(ee) return 0;
    if(em) return 1;
    if(mm) return 2;
    return -1;
  }
  float metcut = 40.;
  if(abs(version)==1) metcut = 60.;
  if(ee&&met.Pt()>metcut&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()>40.) return 0;
  if(mm&&                 (lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()>40.) return 2;
  if(em&&met.Pt()>metcut&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()>30.) {
    if(      maxMT>90.) return 1;
    else if(maxMT<0){
      float tempmaxMT = calcMTmax(tightlep[0],tightlep[1],met);
      if(tempmaxMT>90.) return 1;
    }
  }
  return -1;
}

int isARSS_v0(vector<int> tightlep, vector<int> looselep, vector<int> vetolep,  bool preselect, float maxMT, int nj, int nb, float Mjj, float MjjL, float Detajj, LorentzVector MET, int jec, bool btag, bool Mjjside, int version){
  if(!passJetSSstate_v0(preselect, nj, nb, Mjj, MjjL, Detajj, false, jec, btag, Mjjside,version)) return -1;
  if(tightlep.size()!=1)                                                                          return -1;
  if(looselep.size()!=1)                                                                          return -1;
  if(vetolep.size()!=0)                                                                           return -1;
  if(nisoTrack_mt2_cleaned_VVV_cutbased_veto()!=0)                                                return -1;
  if((lep_pdgId()[tightlep[0] ]*lep_pdgId()[looselep[0] ])<0)                                     return -1;
   LorentzVector met = LorentzVector(0,0,0,0);
  if(MET.Pt()>0)     met = MET;
  else if(jec==1)    met.SetPxPyPzE(met_up_pt()*std::cos(met_up_phi()),met_up_pt()*std::sin(met_up_phi()),0,met_up_pt());
  else if(jec==(-1)) met.SetPxPyPzE(met_dn_pt()*std::cos(met_dn_phi()),met_dn_pt()*std::sin(met_dn_phi()),0,met_dn_pt());
  else               met.SetPxPyPzE(met_pt()*std::cos(met_phi()),met_pt()*std::sin(met_phi()),0,met_pt());
  bool ee = false; bool em = false; bool mm = false;
  if(((lep_pdgId()[tightlep[0] ])*(lep_pdgId()[looselep[0] ]))==121&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[looselep[0] ]).M()-MZ)>10.) ee = true;
  if(((lep_pdgId()[tightlep[0] ])*(lep_pdgId()[looselep[0] ]))==143)  em = true;
  if(((lep_pdgId()[tightlep[0] ])*(lep_pdgId()[looselep[0] ]))==169)  mm = true;
  if(preselect){
    if(ee) return 0;
    if(em) return 1;
    if(mm) return 2;
    return -1;
  }
  float metcut = 40.;
  if(abs(version)==1) metcut = 60.;
  if(ee&&met.Pt()>metcut&&(lep_p4()[tightlep[0] ]+lep_p4()[looselep[0] ]).M()>40.) return 0;
  if(mm&&                 (lep_p4()[tightlep[0] ]+lep_p4()[looselep[0] ]).M()>40.) return 2;
  if(em&&met.Pt()>metcut&&(lep_p4()[tightlep[0] ]+lep_p4()[looselep[0] ]).M()>30.) {
    if(      maxMT>90.) return 1;
    else if(maxMT<0){
      float tempmaxMT = calcMTmax(tightlep[0],looselep[0],met);
      if(tempmaxMT>90.) return 1;
    }
  }
  return -1;
}

int isCRSS_v0(vector<int> tightlep, vector<int> selectlep, vector<int> vetolep, bool preselect, float maxMT, int nj, int nb, float Mjj, float MjjL, float Detajj, LorentzVector MET, int jec, bool noZ, bool btag, bool Mjjside, int version){
  if(!passJetSSstate_v0(preselect, nj, nb, Mjj, MjjL, Detajj, true, jec, btag, Mjjside,version)) return -1;
  if(tightlep.size()  <2)                                                                        return -1;
  if(selectlep.size()!=3)                                                                        return -1;
  if(vetolep.size()!=0)                                                                          return -1;
  if(nisoTrack_mt2_cleaned_VVV_cutbased_veto()!=0)                                               return -1;
  //SS is just tighter pt, so can loop over i3l
  bool isSS = false;
  int lep3 = -1; int SS1 = -1; int SS2 = -1;
  if((lep_pdgId()[selectlep[0] ]*lep_pdgId()[selectlep[1] ])>0) { lep3 = selectlep[2]; SS1 = selectlep[0]; SS2 = selectlep[1]; isSS = true; }
  if((lep_pdgId()[selectlep[0] ]*lep_pdgId()[selectlep[2] ])>0) { lep3 = selectlep[1]; SS1 = selectlep[0]; SS2 = selectlep[2]; isSS = true; }
  if((lep_pdgId()[selectlep[1] ]*lep_pdgId()[selectlep[2] ])>0) { lep3 = selectlep[0]; SS1 = selectlep[1]; SS2 = selectlep[2]; isSS = true; }
  if(!lep_pass_VVV_cutbased_tight()[SS1]) isSS = false;//changed from original function to agree with ntuples
  if(!lep_pass_VVV_cutbased_tight()[SS2]) isSS = false;//changed from original function to agree with ntuples
  if(!isSS)       return -1;
  if(!noZ){
    bool hasSFOSZ = false;
    if((lep_pdgId()[SS1]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[SS1]+lep_p4()[lep3]).M()-MZ)<=10)) hasSFOSZ = true;
    if((lep_pdgId()[SS2]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[SS2]+lep_p4()[lep3]).M()-MZ)<=10)) hasSFOSZ = true;
    if(!hasSFOSZ) return -1;
  } else {
    bool hasSFOS = false;
    if((lep_pdgId()[SS1]==(-lep_pdgId()[lep3]))) hasSFOS = true;
    if((lep_pdgId()[SS2]==(-lep_pdgId()[lep3]))) hasSFOS = true;
    if(!hasSFOS)  return -1;
  }
  bool ee = false; bool em = false; bool mm = false;
  if(((lep_pdgId()[SS1])*(lep_pdgId()[SS2]))==121&&fabs((lep_p4()[SS1]+lep_p4()[SS2]).M()-MZ)>10.) ee = true;
  if(((lep_pdgId()[SS1])*(lep_pdgId()[SS2]))==143)  em = true;
  if(((lep_pdgId()[SS1])*(lep_pdgId()[SS2]))==169)  mm = true;
  if(preselect){
    if(ee) return 0;
    if(em) return 1;
    if(mm) return 2;
    return -1;
  }
  LorentzVector met = LorentzVector(0,0,0,0);
  if(MET.Pt()>0)     met = MET;
  else if(jec==1)    met.SetPxPyPzE(met_up_pt()*std::cos(met_up_phi()),met_up_pt()*std::sin(met_up_phi()),0,met_up_pt());
  else if(jec==(-1)) met.SetPxPyPzE(met_dn_pt()*std::cos(met_dn_phi()),met_dn_pt()*std::sin(met_dn_phi()),0,met_dn_pt());
  else               met.SetPxPyPzE(met_pt()*std::cos(met_phi()),met_pt()*std::sin(met_phi()),0,met_pt());
  float metcut = 40.;
  if(abs(version)==1) metcut = 60.;
  //changed Mll definition from original function to agree with ntuples
  if(ee&&met.Pt()>metcut&&(lep_p4()[SS1]+lep_p4()[SS2]).M()>40.) return 0;
  if(mm&&                 (lep_p4()[SS1]+lep_p4()[SS2]).M()>40.) return 2;
  if(em&&met.Pt()>metcut&&(lep_p4()[SS1]+lep_p4()[SS2]).M()>30.) {
    if(      maxMT>90.) return 1;
    else if(maxMT<0){
      float tempmaxMT = calcMTmax(SS1,SS2,met);
      if(tempmaxMT>90.) return 1;
    }
  }
  return -1;
}

bool passJet3lstate(int jec, bool btag){
  if(jec==1){
    if(nj_up() >1)        return false;
    if( btag&&nb_up()==0) return false;
    if(!btag&&nb_up()!=0) return false; 
  } else if(jec==(-1)){
    if(nj_dn() >1)        return false;
    if( btag&&nb_dn()==0) return false;
    if(!btag&&nb_dn()!=0) return false; 
  } else {
    if(nj() >1)           return false;
    if( btag&&nb()==0)    return false;
    if(!btag&&nb()!=0)    return false; 
  }
  return true;
}

int isSR3l(bool preselect, int jec, bool btag, int version){
  if(!passJet3lstate(jec,btag))               return -1;
  if(Ntight3l()!=3)                           return -1;
  if(nVlep()!=3)                              return -1;
  float met = met_pt();
  if(jec==1)    met = met_up_pt();
  if(jec==(-1)) met = met_dn_pt();
  float dphi = DPhi3lMET();
  if(jec==1)    dphi = DPhi3lMET_up();
  if(jec==(-1)) dphi = DPhi3lMET_dn();
  if(abs(lep_charge()[0]+lep_charge()[1]+lep_charge()[2])==3) return -1;
  if(            Mll3L ()>=0.&&Mll3L ()<20.)  return -1;//is MSF for 0SFOS, else MSFOS
  if(nSFOS()>=2&&Mll3L1()>=0.&&Mll3L1()<20.)  return -1;//2nd MSFOS for 2 SFOS
  if(nSFOS()==0){
    if(fabs(Mee3L()-MZ)<15.)                  return -1;//Mee - is -999 if no ee --> this is fine
    if(preselect)                           return  0;
    if(version==0&&Pt3l()<60.)                return -1;
    if(abs(version)==1&&met<30)               return -1;
  } else if(nSFOS()==1){
    if(version==0&&Mll3L()>55.&&Mll3L()<110.) return -1;
    if(version!=0&&fabs(Mll3L()-MZ)<20.)      return -1;
    if(preselect)                           return  1;
    float metcut = 45.;
    if(abs(version)==1) metcut = 40.;
    if(met<metcut)                            return -1;
    float MTthird = MT3rd();
    if(jec==  1)  MTthird = MT3rd_up();
    if(jec==(-1)) MTthird = MT3rd_dn();
    if(MTthird>=0&&MTthird<90.)               return -1;
  } else if(nSFOS()==2){
    if(fabs(Mll3L ()-MZ)<20.)                 return -1;
    if(fabs(Mll3L1()-MZ)<20.)                 return -1;
    if(preselect)                           return  2;
    if(met<55.)                               return -1;
  }
  if(dphi<2.5)                                return -1;
  if(fabs(M3l()-MZ)<10.)                      return -1;
  if(nSFOS()==0)                            return 0;
  if(Pt3l()<60.)                              return -1;
  if(nSFOS()==1)                            return 1;
  if(nSFOS()==2)                            return 2;
  return -1;
}

int isAR3l(bool preselect, int jec, bool btag, int version){
  if(!passJet3lstate(jec,btag)) return -1;
  if(Ntight3l()!=2)             return -1;
  if(Nloose3l()!=3)             return -1;
  if(nVlep()!=3)                return -1;
  float met = met_pt();
  if(jec==1)    met = met_up_pt();
  if(jec==(-1)) met = met_dn_pt();
  float dphi = DPhi3lMET();
  if(jec==1)    dphi = DPhi3lMET_up();
  if(jec==(-1)) dphi = DPhi3lMET_dn();
  if(abs(lep_charge()[0]+lep_charge()[1]+lep_charge()[2])==3) return -1;
  if(            Mll3L ()>=0.&&Mll3L ()<20.)  return -1;
  if(nSFOS()>=2&&Mll3L1()>=0.&&Mll3L1()<20.)  return -1;
  if(nSFOS()==0){
    if(fabs(Mee3L()-MZ)<15.)                  return -1;
    if(preselect)                           return  0;
    if(version==0&&Pt3l()<60.)                return -1;
    if(abs(version)==1&&met<30)               return -1;
  } else if(nSFOS()==1){
    if(version==0&&Mll3L()>55.&&Mll3L()<110.) return -1;
    if(version!=0&&fabs(Mll3L()-MZ)<20.)      return -1;
    if(preselect)                           return  1;
    float metcut = 45.;
    if(abs(version)==1) metcut = 40.;
    if(met<metcut)                            return -1;
    float MTthird = MT3rd();
    if(jec==  1)  MTthird = MT3rd_up();
    if(jec==(-1)) MTthird = MT3rd_dn();
    if(MTthird>=0&&MTthird<90.)               return -1;
  } else if(nSFOS()==2){
    if(fabs(Mll3L ()-MZ)<20.)                 return -1;
    if(fabs(Mll3L1()-MZ)<20.)                 return -1;
    if(preselect)                           return  2;
    if(met<55.)                               return -1;
  }
  if(dphi<2.5)                                return -1;
  if(fabs(M3l()-MZ)<10.)                      return -1;
  if(nSFOS()==0)                            return 0;
  if(Pt3l()<60.)                              return -1;
  if(nSFOS()==1)                            return 1;
  if(nSFOS()==2)                            return 2;
  return -1;
}

int isCR3l(bool preselect, int jec, bool btag, int version){
  if(!passJet3lstate(jec,btag)) return -1;
  if(Ntight3l()!=3)             return -1;
  if(nVlep()!=3)                return -1;
  float met = met_pt();
  if(jec==1)    met = met_up_pt();
  if(jec==(-1)) met = met_dn_pt();
  float dphi = DPhi3lMET();
  if(jec==1)    dphi = DPhi3lMET_up();
  if(jec==(-1)) dphi = DPhi3lMET_dn();
  if(abs(lep_charge()[0]+lep_charge()[1]+lep_charge()[2])==3) return -1;
  if(            Mll3L ()>=0.&&Mll3L ()<20.)           return -1;
  if(nSFOS()>=2&&Mll3L1()>=0.&&Mll3L1()<20.)           return -1;
  if(nSFOS()==0){
                                                       return -1;
  } else if(nSFOS()==1){
    if(version==0&&!(Mll3L()>55.&&Mll3L()<110.))       return -1;
    if(version!=0&&fabs(Mll3L()-MZ)>=20.)              return -1;
    if(preselect)                                    return  1;
    float metcut = 45.;
    if(abs(version)==1) metcut = 40.;
    if(met<metcut)                                     return -1;
    float MTthird = MT3rd();
    if(jec==  1)  MTthird = MT3rd_up();
    if(jec==(-1)) MTthird = MT3rd_dn();
    if(MTthird>=0&&MTthird<90.)                        return -1;
  } else if(nSFOS()==2){
    if(fabs(Mll3L ()-MZ)>=20.&&fabs(Mll3L1()-MZ)>=20.) return -1;
    if(preselect)                                    return  2;
    if(met<55.)                                        return -1;
  }
  if(dphi<2.5)                                         return -1;
  if(fabs(M3l()-MZ)<10.)                               return -1;
  if(Pt3l()<60.)                                       return -1;
  if(nSFOS()==1)                                      return 1;
  if(nSFOS()==2)                                      return 2;
  return -1;
}

bool passAny3l(int &SR, int &AR, int &CR, bool preselect, int jec, bool btag, int version){
  SR = -1; AR = -1; CR = -1;
  bool isSR = false; bool isAR = false; bool isCR = false; 
  if(!passJet3lstate(jec,btag)) return true;
  if(nVlep()!=3)                return true;
  if(Ntight3l()==3)                { isSR = true; isCR = true; }
  if(Ntight3l()==2&&Nloose3l()==3)   isAR = true;
  float met = met_pt();
  if(jec==1)    met = met_up_pt();
  if(jec==(-1)) met = met_dn_pt();
  float dphi = DPhi3lMET();
  if(jec==1)    dphi = DPhi3lMET_up();
  if(jec==(-1)) dphi = DPhi3lMET_dn();
  if(abs(lep_charge()[0]+lep_charge()[1]+lep_charge()[2])==3) return true;
  if(            Mll3L ()>=0.&&Mll3L ()<20.) return true;
  if(nSFOS()>=2&&Mll3L1()>=0.&&Mll3L1()<20.) return true;
  if(nSFOS()==0){
    if(fabs(Mee3L()-MZ)<15.)    return true;
    if(preselect){
      if(isSR) SR = 0;
      if(isAR) AR = 0;
      return true;
    }
    if(version==0&&Pt3l()<60.)  return true;
    if(abs(version)==1&&met<30) return true;
  } else if(nSFOS()==1){
    if(version==0&&Mll3L()>55.&&Mll3L()<110.) { isSR = false; isAR = false; }
    else if(version==0)                         isCR = false;
    if(version!=0&&fabs(Mll3L()-MZ)<=20.)     { isSR = false; isAR = false; }
    else if(version!=0)                         isCR = false;
    if(preselect) {
      if(isSR) SR = 1;
      if(isAR) AR = 1;
      if(isCR) CR = 1;
      return true;
    }
    float metcut = 45.;
    if(abs(version)==1) metcut = 40.;
    if(met<metcut)                            return true;
    float MTthird = MT3rd();
    if(jec==  1)  MTthird = MT3rd_up();
    if(jec==(-1)) MTthird = MT3rd_dn();
    if(MTthird>=0&&MTthird<90.)               return true;
  } else if(nSFOS()==2){
    if(fabs(Mll3L ()-MZ)>20.&&fabs(Mll3L1()-MZ)>20.)   isCR = false;
    else                                             { isSR = false; isAR = false; }
    if(preselect) {
      if(isSR) SR = 2;
      if(isAR) AR = 2;
      if(isCR) CR = 2;
      return true;
    }
    if(met<55.)                                        return true;
  }
  if(dphi<2.5)           return true;
  if(fabs(M3l()-MZ)<10.) return true;
  if(nSFOS()==0){
      if(isSR) SR = 0;
      if(isAR) AR = 0;
      return true;
    }
  if(Pt3l()<60.)         return true;
  if(nSFOS()==1){
      if(isSR) SR = 1;
      if(isAR) AR = 1;
      if(isCR) CR = 1;
      return true;
    }
  if(nSFOS()==2){
      if(isSR) SR = 2;
      if(isAR) AR = 2;
      if(isCR) CR = 2;
      return true;
    }
  return false;
}
		      
int isSR3l_v0(vector<int> tightlep, bool preselect, int nj, int nb, LorentzVector MET, int jec, bool btag, int version){
  int numj = nj;
  int numb = nb;
  int x;
  if(numj<0||numb<0) getalljetnumbers(numj,x,numb, jec);
  if(numj >1)            return -1;
  if( btag&&numb==0)     return -1;
  if(!btag&&numb!=0)     return -1;
  if(tightlep.size()!=3) return -1;
  float DPhilllMET = -1;
  float pTlll = -1;
  LorentzVector met = LorentzVector(0,0,0,0);
  if(!preselect){
    if(MET.Pt()>0)     met = MET;
    else if(jec==1)    met.SetPxPyPzE(met_up_pt()*std::cos(met_up_phi()),met_up_pt()*std::sin(met_up_phi()),0,met_up_pt());
    else if(jec==(-1)) met.SetPxPyPzE(met_dn_pt()*std::cos(met_dn_phi()),met_dn_pt()*std::sin(met_dn_phi()),0,met_dn_pt());
    else               met.SetPxPyPzE(met_pt()*std::cos(met_phi()),met_pt()*std::sin(met_phi()),0,met_pt());
    if(fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()-MZ)<10.) return -1;
    pTlll = (lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).Pt();
    DPhilllMET = deltaPhi((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).Phi(),met.Phi());
  }
  else DPhilllMET = 3.2;//pass DPhilllMET if preselect
  bool OS01 = (lep_pdgId()[tightlep[0] ]*lep_pdgId()[tightlep[1] ]<0); bool SF01 = (abs(lep_pdgId()[tightlep[0] ])==abs(lep_pdgId()[tightlep[1] ]));
  bool OS02 = (lep_pdgId()[tightlep[0] ]*lep_pdgId()[tightlep[2] ]<0); bool SF02 = (abs(lep_pdgId()[tightlep[0] ])==abs(lep_pdgId()[tightlep[2] ]));
  bool OS12 = (lep_pdgId()[tightlep[1] ]*lep_pdgId()[tightlep[2] ]<0); bool SF12 = (abs(lep_pdgId()[tightlep[1] ])==abs(lep_pdgId()[tightlep[2] ]));
  int SFOScounter = 0;
  if(OS01&&SF01) ++SFOScounter;
  if(OS02&&SF02) ++SFOScounter;
  if(OS12&&SF12) ++SFOScounter;
  bool pass0(false), pass1(false), pass2(false);   //SR
  if(OS01&&SF01&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()<20.) SFOScounter = -1;
  if(OS02&&SF02&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()<20.) SFOScounter = -1;
  if(OS12&&SF12&&(lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()<20.) SFOScounter = -1;
  if(abs(lep_charge()[tightlep[0] ]+lep_charge()[tightlep[1] ]+lep_charge()[tightlep[2] ])==3) SFOScounter = -1;
  if(SFOScounter==0){
    pass0 = true;
    if(SF01&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()<20.) { pass0 = false; }
    if(SF02&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()<20.) { pass0 = false; }
    if(SF12&&(lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()<20.) { pass0 = false; }
    if(SF01&&abs(lep_pdgId()[tightlep[0] ])==11&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()-MZ)<15.) { pass0 = false; }
    if(SF02&&abs(lep_pdgId()[tightlep[0] ])==11&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()-MZ)<15.) { pass0 = false; }
    if(SF12&&abs(lep_pdgId()[tightlep[1] ])==11&&fabs((lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()-MZ)<15.) { pass0 = false; }
    if(version==0){
      if(!preselect&&pTlll<60.) return -1;
    }
    if(abs(version)==1){
      if(!preselect&&met.Pt()<30) return -1;
    }
  }
  if(SFOScounter==1){
    if(!preselect&&pTlll<60.) return -1;
    pass1 = true;
    float metcut = 45.;
    if(abs(version)==1) metcut = 40.;
    if(!preselect&&met.Pt()<metcut) return -1;
    if(version==0){
      if(OS01&&SF01&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()>55.&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()<110.) { pass1 = false; }
      if(OS02&&SF02&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()>55.&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()<110.) { pass1 = false; }
      if(OS12&&SF12&&(lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()>55.&&(lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()<110.) { pass1 = false; }
    }
    if(abs(version)==1){
      int thirdlep = -1;
      if(OS01&&SF01) thirdlep = tightlep[2];
      if(OS02&&SF02) thirdlep = tightlep[1];
      if(OS12&&SF12) thirdlep = tightlep[0];
      if(OS01&&SF01&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()-MZ)<20.) { pass1 = false; }
      if(OS02&&SF02&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()-MZ)<20.) { pass1 = false; }
      if(OS12&&SF12&&fabs((lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()-MZ)<20.) { pass1 = false; }
      if(!preselect&&thirdlep>=0&&mT(lep_p4()[thirdlep],met)<90.) return -1;
    }
  }
  if(SFOScounter==2){
    if(!preselect&&pTlll<60.) return -1;
    pass2 = true;
    if(!preselect&&met.Pt()<55) return -1;
    if(OS01&&SF01&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()-MZ)<20.) { pass2 = false; }
    if(OS02&&SF02&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()-MZ)<20.) { pass2 = false; }
    if(OS12&&SF12&&fabs((lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()-MZ)<20.) { pass2 = false; }
  }
  if(DPhilllMET>=2.5&&pass0 ) return 0;//2.7 --> 2.5 changed from original function
  if(DPhilllMET>=2.5&&pass1 ) return 1;
  if(DPhilllMET>=2.5&&pass2 ) return 2;
  return -1;
}

int isAR3l_v0(vector<int> tightlep, vector<int> looselep, bool preselect, int nj, int nb, LorentzVector MET, int jec, bool btag, int version){
  int numj = nj;
  int numb = nb;
  int x;
  if(numj<0||numb<0) getalljetnumbers(numj,x,numb, jec);
  if(numj >1)            return -1;
  if( btag&&numb==0)     return -1;
  if(!btag&&numb!=0)     return -1;
  if(tightlep.size()!=2) return -1;
  if(looselep.size()!=1) return -1;
  float DPhilllMET = -1;
  float pTlll = -1;
  LorentzVector met = LorentzVector(0,0,0,0);
  if(!preselect){
    if(MET.Pt()>0)     met = MET;
    else if(jec==1)    met.SetPxPyPzE(met_up_pt()*std::cos(met_up_phi()),met_up_pt()*std::sin(met_up_phi()),0,met_up_pt());
    else if(jec==(-1)) met.SetPxPyPzE(met_dn_pt()*std::cos(met_dn_phi()),met_dn_pt()*std::sin(met_dn_phi()),0,met_dn_pt());
    else               met.SetPxPyPzE(met_pt()*std::cos(met_phi()),met_pt()*std::sin(met_phi()),0,met_pt());
    if(fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]+lep_p4()[looselep[0] ]).M()-MZ)<10.) return -1;
    pTlll = (lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]+lep_p4()[looselep[0] ]).Pt();
    DPhilllMET = deltaPhi((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]+lep_p4()[looselep[0] ]).Phi(),met.Phi());
  }
  else DPhilllMET = 3.2;//pass DPhilllMET if preselect
  bool OS01 = (lep_pdgId()[tightlep[0] ]*lep_pdgId()[tightlep[1] ]<0); bool SF01 = (abs(lep_pdgId()[tightlep[0] ])==abs(lep_pdgId()[tightlep[1] ]));
  bool OS02 = (lep_pdgId()[tightlep[0] ]*lep_pdgId()[looselep[0] ]<0); bool SF02 = (abs(lep_pdgId()[tightlep[0] ])==abs(lep_pdgId()[looselep[0] ]));
  bool OS12 = (lep_pdgId()[tightlep[1] ]*lep_pdgId()[looselep[0] ]<0); bool SF12 = (abs(lep_pdgId()[tightlep[1] ])==abs(lep_pdgId()[looselep[0] ]));
  int SFOScounter = 0;
  if(OS01&&SF01) ++SFOScounter;
  if(OS02&&SF02) ++SFOScounter;
  if(OS12&&SF12) ++SFOScounter;
  bool pass0(false), pass1(false), pass2(false);   //SR
  if(OS01&&SF01&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()<20.) SFOScounter = -1;
  if(OS02&&SF02&&(lep_p4()[tightlep[0] ]+lep_p4()[looselep[0] ]).M()<20.) SFOScounter = -1;
  if(OS12&&SF12&&(lep_p4()[tightlep[1] ]+lep_p4()[looselep[0] ]).M()<20.) SFOScounter = -1;
  if(abs(lep_charge()[tightlep[0] ]+lep_charge()[tightlep[1] ]+lep_charge()[looselep[0] ])==3) SFOScounter = -1;
  if(SFOScounter==0){
    pass0 = true;
    if(SF01&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()<20.) { pass0 = false; }
    if(SF02&&(lep_p4()[tightlep[0] ]+lep_p4()[looselep[0] ]).M()<20.) { pass0 = false; }
    if(SF12&&(lep_p4()[tightlep[1] ]+lep_p4()[looselep[0] ]).M()<20.) { pass0 = false; }
    if(SF01&&abs(lep_pdgId()[tightlep[0] ])==11&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()-MZ)<15.) { pass0 = false; }
    if(SF02&&abs(lep_pdgId()[tightlep[0] ])==11&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[looselep[0] ]).M()-MZ)<15.) { pass0 = false; }
    if(SF12&&abs(lep_pdgId()[tightlep[1] ])==11&&fabs((lep_p4()[tightlep[1] ]+lep_p4()[looselep[0] ]).M()-MZ)<15.) { pass0 = false; }
    if(version==0){
      if(!preselect&&pTlll<60.) return -1;
    }
    if(abs(version)==1){
      if(!preselect&&met.Pt()<30) return -1;
    }
  }
  if(SFOScounter==1){
    if(!preselect&&pTlll<60.) return -1;
    pass1 = true;
    float metcut = 45.;
    if(abs(version)==1) metcut = 40.;
    if(!preselect&&met.Pt()<metcut) return -1;
    if(version==0){
      if(OS01&&SF01&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()>55.&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()<110.) { pass1 = false; }
      if(OS02&&SF02&&(lep_p4()[tightlep[0] ]+lep_p4()[looselep[0] ]).M()>55.&&(lep_p4()[tightlep[0] ]+lep_p4()[looselep[0] ]).M()<110.) { pass1 = false; }
      if(OS12&&SF12&&(lep_p4()[tightlep[1] ]+lep_p4()[looselep[0] ]).M()>55.&&(lep_p4()[tightlep[1] ]+lep_p4()[looselep[0] ]).M()<110.) { pass1 = false; }
    }
    if(abs(version)==1){
      int thirdlep = -1;
      if(OS01&&SF01) thirdlep = looselep[0];
      if(OS02&&SF02) thirdlep = tightlep[1];
      if(OS12&&SF12) thirdlep = tightlep[0];
      if(OS01&&SF01&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()-MZ)<20.) { pass1 = false; }
      if(OS02&&SF02&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[looselep[0] ]).M()-MZ)<20.) { pass1 = false; }
      if(OS12&&SF12&&fabs((lep_p4()[tightlep[1] ]+lep_p4()[looselep[0] ]).M()-MZ)<20.) { pass1 = false; }
      if(!preselect&&thirdlep>=0&&mT(lep_p4()[thirdlep],met)<90.) return -1;
   }
  }
  if(SFOScounter==2){
    if(!preselect&&pTlll<60.) return -1;
    pass2 = true;
    if(!preselect&&met.Pt()<55) return -1;
    if(OS01&&SF01&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()-MZ)<20.) { pass2 = false; }
    if(OS02&&SF02&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[looselep[0] ]).M()-MZ)<20.) { pass2 = false; }
    if(OS12&&SF12&&fabs((lep_p4()[tightlep[1] ]+lep_p4()[looselep[0] ]).M()-MZ)<20.) { pass2 = false; }
  }
  if(DPhilllMET>=2.5&&pass0 ) return 0;//2.7 --> 2.5 changed from original function
  if(DPhilllMET>=2.5&&pass1 ) return 1;
  if(DPhilllMET>=2.5&&pass2 ) return 2;
  return -1;
}

int isCR3l_v0(vector<int> tightlep, bool preselect, int nj, int nb, LorentzVector MET, int jec, bool btag, int version){
  int numj = nj;
  int numb = nb;
  int x;
  if(numj<0||numb<0) getalljetnumbers(numj,x,numb, jec);
  if(numj >1)            return -1;
  if( btag&&numb==0)     return -1;
  if(!btag&&numb!=0)     return -1;
  if(tightlep.size()!=3) return -1;
  float DPhilllMET = -1;
  float pTlll = -1;
  LorentzVector met = LorentzVector(0,0,0,0);
  if(!preselect){
    if(MET.Pt()>0)     met = MET;
    else if(jec==1)    met.SetPxPyPzE(met_up_pt()*std::cos(met_up_phi()),met_up_pt()*std::sin(met_up_phi()),0,met_up_pt());
    else if(jec==(-1)) met.SetPxPyPzE(met_dn_pt()*std::cos(met_dn_phi()),met_dn_pt()*std::sin(met_dn_phi()),0,met_dn_pt());
    else               met.SetPxPyPzE(met_pt()*std::cos(met_phi()),met_pt()*std::sin(met_phi()),0,met_pt());
    if(fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()-MZ)<10.) return -1;
    pTlll = (lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).Pt();
    DPhilllMET = deltaPhi((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).Phi(),met.Phi());
  }
  else DPhilllMET = 3.2;//pass DPhilllMET if preselect
  bool OS01 = (lep_pdgId()[tightlep[0] ]*lep_pdgId()[tightlep[1] ]<0); bool SF01 = (abs(lep_pdgId()[tightlep[0] ])==abs(lep_pdgId()[tightlep[1] ]));
  bool OS02 = (lep_pdgId()[tightlep[0] ]*lep_pdgId()[tightlep[2] ]<0); bool SF02 = (abs(lep_pdgId()[tightlep[0] ])==abs(lep_pdgId()[tightlep[2] ]));
  bool OS12 = (lep_pdgId()[tightlep[1] ]*lep_pdgId()[tightlep[2] ]<0); bool SF12 = (abs(lep_pdgId()[tightlep[1] ])==abs(lep_pdgId()[tightlep[2] ]));
  int SFOScounter = 0;
  if(OS01&&SF01) ++SFOScounter;
  if(OS02&&SF02) ++SFOScounter;
  if(OS12&&SF12) ++SFOScounter;
  bool pass0X(false), pass1X(false), pass2X(false);//CR-WZ
  if(OS01&&SF01&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()<20.) SFOScounter = -1;
  if(OS02&&SF02&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()<20.) SFOScounter = -1;
  if(OS12&&SF12&&(lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()<20.) SFOScounter = -1;
  if(abs(lep_charge()[tightlep[0] ]+lep_charge()[tightlep[1] ]+lep_charge()[tightlep[2] ])==3) SFOScounter = -1;
  //SFOScounter==0 has no CR
  if(SFOScounter==1){
    if(!preselect&&pTlll<60.) return -1;
    pass1X = true;
    float metcut = 45.;
    if(abs(version)==1) metcut = 40.;
    if(!preselect&&met.Pt()<metcut) return -1;
    bool atleastoneSFOSZ = false;
    if(version==0){
      if(OS01&&SF01&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()>55.&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()<110.) { atleastoneSFOSZ = true; }
      if(OS02&&SF02&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()>55.&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()<110.) { atleastoneSFOSZ = true; }
      if(OS12&&SF12&&(lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()>55.&&(lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()<110.) { atleastoneSFOSZ = true; }
    }
    if(abs(version)==1){
      int thirdlep = -1;
      if(OS01&&SF01) thirdlep = tightlep[2];
      if(OS02&&SF02) thirdlep = tightlep[1];
      if(OS12&&SF12) thirdlep = tightlep[0];
      if(OS01&&SF01&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()-MZ)<20.) { atleastoneSFOSZ = true; }
      if(OS02&&SF02&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()-MZ)<20.) { atleastoneSFOSZ = true; }
      if(OS12&&SF12&&fabs((lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()-MZ)<20.) { atleastoneSFOSZ = true; }
      if(!preselect&&thirdlep>=0&&mT(lep_p4()[thirdlep],met)<90.) return -1;
    }
    if(!atleastoneSFOSZ) pass1X = false;
  }
  if(SFOScounter==2){
    if(!preselect&&pTlll<60.) return -1;
    pass2X = true;
    if(!preselect&&met.Pt()<55) { pass2X=false; }
    bool atleastoneSFOSZ = false;
    if(OS01&&SF01&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()-MZ)<20.) { atleastoneSFOSZ = true; }
    if(OS02&&SF02&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()-MZ)<20.) { atleastoneSFOSZ = true; }
    if(OS12&&SF12&&fabs((lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()-MZ)<20.) { atleastoneSFOSZ = true; }
    if(!atleastoneSFOSZ) pass2X = false;
  }
  if(DPhilllMET>=2.5&&pass0X ) return 0;//2.7 --> 2.5 changed from original function
  if(DPhilllMET>=2.5&&pass1X ) return 1;
  if(DPhilllMET>=2.5&&pass2X ) return 2;
  return -1;
}

bool checkbothSRCR3l_v0(int &isSR3l, int &isCR3l, vector<int> tightlep, bool preselect, int nj, int nb, LorentzVector MET, int jec, bool btag, int version){
  isSR3l = -1; isCR3l = -1;//assume no match
  int numj = nj;
  int numb = nb;
  int x;
  if(numj<0||numb<0) getalljetnumbers(numj,x,numb, jec);
  if(numj >1)            return true;
  if( btag&&numb==0)     return true;
  if(!btag&&numb!=0)     return true;
  if(tightlep.size()!=3) return true;
  float DPhilllMET = -1;
  float pTlll = -1;
  LorentzVector met = LorentzVector(0,0,0,0);
  if(!preselect){
    if(MET.Pt()>0)     met = MET;
    else if(jec==1)    met.SetPxPyPzE(met_up_pt()*std::cos(met_up_phi()),met_up_pt()*std::sin(met_up_phi()),0,met_up_pt());
    else if(jec==(-1)) met.SetPxPyPzE(met_dn_pt()*std::cos(met_dn_phi()),met_dn_pt()*std::sin(met_dn_phi()),0,met_dn_pt());
    else               met.SetPxPyPzE(met_pt()*std::cos(met_phi()),met_pt()*std::sin(met_phi()),0,met_pt());
    if(fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()-MZ)<10.) return true;
    pTlll = (lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).Pt();
    DPhilllMET = deltaPhi((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).Phi(),met.Phi());
  }
  else DPhilllMET = 3.2;//pass DPhilllMET if preselect
  bool OS01 = (lep_pdgId()[tightlep[0] ]*lep_pdgId()[tightlep[1] ]<0); bool SF01 = (abs(lep_pdgId()[tightlep[0] ])==abs(lep_pdgId()[tightlep[1] ]));
  bool OS02 = (lep_pdgId()[tightlep[0] ]*lep_pdgId()[tightlep[2] ]<0); bool SF02 = (abs(lep_pdgId()[tightlep[0] ])==abs(lep_pdgId()[tightlep[2] ]));
  bool OS12 = (lep_pdgId()[tightlep[1] ]*lep_pdgId()[tightlep[2] ]<0); bool SF12 = (abs(lep_pdgId()[tightlep[1] ])==abs(lep_pdgId()[tightlep[2] ]));
  int SFOScounter = 0;
  if(OS01&&SF01) ++SFOScounter;
  if(OS02&&SF02) ++SFOScounter;
  if(OS12&&SF12) ++SFOScounter;
  bool pass0(false), pass1(false), pass2(false);   //SR
  bool pass0X(false), pass1X(false), pass2X(false);//CR-WZ
  if(OS01&&SF01&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()<20.) SFOScounter = -1;
  if(OS02&&SF02&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()<20.) SFOScounter = -1;
  if(OS12&&SF12&&(lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()<20.) SFOScounter = -1;
  if(abs(lep_charge()[tightlep[0] ]+lep_charge()[tightlep[1] ]+lep_charge()[tightlep[2] ])==3) SFOScounter = -1;
  if(SFOScounter==0){
    pass0 = true; pass0X = true;
    if(SF01&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()<20.) { pass0 = false; pass0X = false; }
    if(SF02&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()<20.) { pass0 = false; pass0X = false; }
    if(SF12&&(lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()<20.) { pass0 = false; pass0X = false; }
    bool atleastoneSFOSZ = false;
    if(SF01&&abs(lep_pdgId()[tightlep[0] ])==11&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()-MZ)<15.) { pass0 = false; if(OS01) atleastoneSFOSZ = true; }
    if(SF02&&abs(lep_pdgId()[tightlep[0] ])==11&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()-MZ)<15.) { pass0 = false; if(OS02) atleastoneSFOSZ = true; }
    if(SF12&&abs(lep_pdgId()[tightlep[1] ])==11&&fabs((lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()-MZ)<15.) { pass0 = false; if(OS12) atleastoneSFOSZ = true; }
    if(!atleastoneSFOSZ) pass0X = false;//this is always false :)
    if(version==0){
      if(!preselect&&pTlll<60.) return false;
    }
    if(abs(version)==1){
      if(!preselect&&met.Pt()<30) return false;
    }
  }
  if(SFOScounter==1){
    if(!preselect&&pTlll<60.) return false;
    float metcut = 45.;
    if(abs(version)==1) metcut = 40.;
    if(!preselect&&met.Pt()<metcut) return false;
    pass1 = true; pass1X = true;
    bool atleastoneSFOSZ = false;
    if(version==0){
      if(OS01&&SF01&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()>55.&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
      if(OS02&&SF02&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()>55.&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
      if(OS12&&SF12&&(lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()>55.&&(lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
    }
    if(abs(version)==1){
      int thirdlep = -1;
      if(OS01&&SF01) thirdlep = tightlep[2];
      if(OS02&&SF02) thirdlep = tightlep[1];
      if(OS12&&SF12) thirdlep = tightlep[0];
      if(OS01&&SF01&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()-MZ)<20.) { pass1 = false; atleastoneSFOSZ = true; }
      if(OS02&&SF02&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()-MZ)<20.) { pass1 = false; atleastoneSFOSZ = true; }
      if(OS12&&SF12&&fabs((lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()-MZ)<20.) { pass1 = false; atleastoneSFOSZ = true; }
      if(!preselect&&thirdlep>=0&&mT(lep_p4()[thirdlep],met)<90.) return -1;
    }
    if(!atleastoneSFOSZ) pass1X = false;
  }
  if(SFOScounter==2){
    if(!preselect&&pTlll<60.) return false;
    pass2 = true; pass2X = true;
    if(!preselect&&met.Pt()<55) { pass2=false; pass2X=false; }
    bool atleastoneSFOSZ = false;
    if(OS01&&SF01&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()-MZ)<20.) { pass2 = false; atleastoneSFOSZ = true; }
    if(OS02&&SF02&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()-MZ)<20.) { pass2 = false; atleastoneSFOSZ = true; }
    if(OS12&&SF12&&fabs((lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()-MZ)<20.) { pass2 = false; atleastoneSFOSZ = true; }
    if(!atleastoneSFOSZ) pass2X = false;
  }
  //SR
  if(DPhilllMET>=2.5&&pass0 ) isSR3l = 0;//2.7 --> 2.5 changed from original function
  if(DPhilllMET>=2.5&&pass1 ) isSR3l = 1;
  if(DPhilllMET>=2.5&&pass2 ) isSR3l = 2;
  //CR
  if(DPhilllMET>=2.5&&pass0X) isCR3l = 0;//2.7 --> 2.5 changed from original function
  if(DPhilllMET>=2.5&&pass1X) isCR3l = 1;
  if(DPhilllMET>=2.5&&pass2X) isCR3l = 2;
  return true;
}

vector<float> allMSFOS(vector<int> tightlep, vector<int> looselep){
  vector<int> l = tightlep;
  l.insert( l.end(), looselep.begin(), looselep.end() );
  vector<float> MSFOS;
  for(unsigned int i = 0; i<l.size();++i){
    for(unsigned int j = i+1; j<l.size(); ++j){
      if((lep_pdgId()[l[i] ])!=(-lep_pdgId()[l[j] ])) continue;
      MSFOS.push_back((lep_p4()[l[i] ]+lep_p4()[l[j] ]).M());
    }
  }
  sort(MSFOS.begin(),MSFOS.end(),sortMSFOS);
  return MSFOS;
}
vector<float> allMSFOS(){
  vector<float> MSFOS;
  if(isSFOS01()==1) MSFOS.push_back(M01());
  if(isSFOS02()==1) MSFOS.push_back(M02());
  if(isSFOS12()==1) MSFOS.push_back(M12());
  sort(MSFOS.begin(),MSFOS.end(),sortMSFOS);
  return MSFOS;
}

map<string, TH1D*> bookhistclosure(string samplename, vector<string> histonames, vector<int> hbins, vector<float> hlow, vector<float> hup, TDirectory *rootdir){
  map<string, TH1D*> histos;
  if(histonames.size()!=hbins.size()) return histos;
  if(histonames.size()!=hlow.size())  return histos;
  if(histonames.size()!=hup.size())   return histos;
  for(unsigned int i = 0; i<histonames.size(); ++i){
    string mapname = histonames[i];
    mapname = histonames[i] + "_"+samplename;
    if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    mapname = histonames[i] + "_bfake";
    if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    mapname = histonames[i] + "_cfake";
    if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    mapname = histonames[i] + "_lightfake";
    if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    mapname = histonames[i] + "_gamma";
    if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    mapname = histonames[i] + "_real";
    if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    mapname = histonames[i] + "_chargeflip";
    if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    mapname = histonames[i] + "_unknown";
    if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
  }
  for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
    h->second->Sumw2(); h->second->SetDirectory(rootdir);
  }
  return histos;
}

map<string, TH1D*> bookhistograms(string samplename, vector<string> histonames, vector<int> hbins, vector<float> hlow, vector<float> hup, TDirectory *rootdir, int splitWW){
  map<string, TH1D*> histos;
  if(histonames.size()!=hbins.size()) return histos;
  if(histonames.size()!=hlow.size())  return histos;
  if(histonames.size()!=hup.size())   return histos;
  for(unsigned int i = 0; i<histonames.size(); ++i){
    string mapname = histonames[i];
    mapname = histonames[i] + "_"+samplename;
    if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    mapname = histonames[i] + "_fakesPred";
    if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    mapname = histonames[i] + "_Other";
    if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    mapname = histonames[i] + "_WHtoWWW";
    if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    mapname = histonames[i] + "_trueSS";
    if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    mapname = histonames[i] + "_chargeflips";
    if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    mapname = histonames[i] + "_SSLL";
    if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    mapname = histonames[i] + "_fakes";
    if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    mapname = histonames[i] + "_photonfakes";
    if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    mapname = histonames[i] + "_others";
    if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    mapname = histonames[i] + "_trueWWW";
    if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    mapname = histonames[i] + "_3lLL";
    if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    mapname = histonames[i] + "_true3L";
    if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    if(splitWW==1){
      mapname = histonames[i] + "_WWRest";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_otherWWSS";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_WWVBS";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_WWDPS";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_ttVSS";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_OtherttV";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_ttV3l";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    }
    else if(splitWW==2){
      mapname = histonames[i] + "_WWVBS";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_ttW";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_ttZ";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    }
  }
  for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
    h->second->Sumw2(); h->second->SetDirectory(rootdir);
  }
  return histos;
}

bool deleteHistograms(map<string, TH1D*> histos){
  for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h) delete h->second;
  return true;
}

bool fillSRhisto(map<string, TH1D*> histos, string histoname, string sample, string sn, string sn2, int SRSS, int SR3l, float weight, float weight3l, bool fillsample){
  float weight2 = weight3l;
  if(weight3l<-1e6) weight2 = weight;
  if(fillsample){
    if(SRSS==0)   histos[histoname+"_"+sample]->Fill(0.,weight);
    if(SRSS==1)   histos[histoname+"_"+sample]->Fill(1.,weight);
    if(SRSS==2)   histos[histoname+"_"+sample]->Fill(2.,weight);
    if(SR3l==0)   histos[histoname+"_"+sample]->Fill(3.,weight2);
    if(SR3l==1)   histos[histoname+"_"+sample]->Fill(4.,weight2);
    if(SR3l==2)   histos[histoname+"_"+sample]->Fill(5.,weight2);
  }
  if(sample!=sn){
    if(SRSS==0) histos[histoname+"_"+sn    ]->Fill(0.,weight);
    if(SRSS==1) histos[histoname+"_"+sn    ]->Fill(1.,weight);
    if(SRSS==2) histos[histoname+"_"+sn    ]->Fill(2.,weight);
  }
  if(sample!=sn2){
    if(SR3l==0) histos[histoname+"_"+sn2   ]->Fill(3.,weight2);
    if(SR3l==1) histos[histoname+"_"+sn2   ]->Fill(4.,weight2);
    if(SR3l==2) histos[histoname+"_"+sn2   ]->Fill(5.,weight2);
  }
  if(SRSS>=0||SR3l>=0) return true;
  else return false;
}

bool fillhisto(map<string, TH1D*> histos, string histoname, string sample, string sn, float value, float weight, bool fillsample){
  if(fillsample) histos[  histoname+"_"+sample]->Fill(value,weight);
  if(sample!=sn)
    histos[histoname+"_"+sn    ]->Fill(value,weight);
  return true;
}

bool SaveHistosToFile(string filename, map<string, TH1D*> histos, bool addunderflow, bool addoverflow, bool deletefile){
  for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
    //add overflow
    if(addoverflow){
      h->second->SetBinContent(h->second->GetNbinsX(),          h->second->GetBinContent(h->second->GetNbinsX() )   +    h->second->GetBinContent(h->second->GetNbinsX()+1) );
      h->second->SetBinError(  h->second->GetNbinsX(), sqrt(pow(h->second->GetBinError(  h->second->GetNbinsX() ),2)+pow(h->second->GetBinError(  h->second->GetNbinsX()+1),2) ) );
    }
    //add underflow
    if(addunderflow){
      h->second->SetBinContent(1,          h->second->GetBinContent(1)   +    h->second->GetBinContent(0) );
      h->second->SetBinError(  1, sqrt(pow(h->second->GetBinError(  1),2)+pow(h->second->GetBinError(  0),2) ) );
    }
  }
  
  TFile *f;
  if(deletefile) f = new TFile(filename.c_str(),"RECREATE");
  else           f = new TFile(filename.c_str(),"UPDATE");
  f->cd();
  for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
    if(f->GetListOfKeys()->Contains(h->first.c_str())){
      h->second->SetName((h->first+"new").c_str());
      TH1D *hold = (TH1D*)f->Get(h->first.c_str());
      hold->SetName((h->first+"old").c_str());
      TH1D *hnew = (TH1D*)hold->Clone(h->first.c_str());
      hnew->Add(h->second);
      hnew->Write(h->first.c_str(),TObject::kOverwrite);
    }
    else
      h->second->Write(h->first.c_str(),TObject::kOverwrite);
  }
  f->Close();
  cout << "Saved histos in " << f->GetName() << endl;
  delete f;
  return true;
}

bool fileexists(string filename){
  ifstream f(filename.c_str());
  return f.good();
}

//only used for tight right now - loose region is completely data driven
float getlepSFandError(float &error, int index){
  error = 0.;
  if(index<0)                        return 0;
  if(index>=(int)lep_pdgId().size()) return 0;
  if(abs(lep_pdgId()[index])==11){
    float effreco  = lepsf_EGammaReco(                  lep_p4()[index].Pt(),lep_etaSC()[index]);
    float efftight = lepsf_EGammaTightID(               lep_p4()[index].Pt(),lep_etaSC()[index]);
    float effWWW   = lepsf_EGammaTightPOG_EGammaVVV(    lep_p4()[index].Pt(),lep_etaSC()[index]);
    float errreco  = lepsf_EGammaReco_unc(              lep_p4()[index].Pt(),lep_etaSC()[index]);
    float errtight = lepsf_EGammaTightID_unc(           lep_p4()[index].Pt(),lep_etaSC()[index]);
    float errWWW   = 0.01*fabs(1.-effWWW);//add 1% of difference to 1.
    error = sqrt(pow(errreco*efftight*effWWW,2)+pow(errtight*effreco*effWWW,2)+pow(errWWW*effreco*efftight,2));
    return effreco*efftight*effWWW;
  }
  else {//muon
    float effreco   = lepsf_MuReco(                lep_p4()[index].Pt(),lep_p4()[index].Eta());
    float efftight1 = lepsf_MuMediumID_BtoF(       lep_p4()[index].Pt(),lep_p4()[index].Eta());
    float efftight2 = lepsf_MuMediumID_GH(         lep_p4()[index].Pt(),lep_p4()[index].Eta());
    float efftight  = efftight1*0.549833+efftight2*0.450167;//luminosity weights B-F vs. G+H
    float effWWW    = lepsf_MuMediumPOG_MuTightVVV(lep_p4()[index].Pt(),lep_p4()[index].Eta());
    float errreco   = lepsf_MuReco_unc(            lep_p4()[index].Pt(),lep_p4()[index].Eta());
    //recommendation from https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults and linked twiki
    float errtight1 = sqrt(pow(lepsf_MuMediumID_BtoF_unc(   lep_p4()[index].Pt(),lep_p4()[index].Eta()),2)+pow(0.01*fabs(1.-efftight1),2));
    float errtight2 = sqrt(pow(lepsf_MuMediumID_GH_unc(     lep_p4()[index].Pt(),lep_p4()[index].Eta()),2)+pow(0.01*fabs(1.-efftight2),2));
    float errtight  = sqrt(pow(errtight1*0.549833,2)+pow(errtight2*0.450167,2));
    float errWWW    = 0.01*fabs(1.-effWWW);//add 1% of difference to 1.
    error = sqrt(pow(errreco*efftight*effWWW,2)+pow(errtight*effreco*effWWW,2)+pow(errWWW*effreco*efftight,2));
    return effreco*efftight*effWWW;
  }
}

float getlepSFWeightandError(float &error, vector<int> tightlep, vector<int> looselep){//loose not implemented yet
  //currently I use only tightlep SF, I don't use looselep at all!
  error = 0;
  if(tightlep.size()==0) return 1;
  vector<float> myeff, myerr;
  float tempeff(0),temperr(0);
  float SF = 1.;
  for(unsigned int i = 0; i<tightlep.size();++i){
    tempeff = getlepSFandError(temperr,tightlep[i]);
    myeff.push_back(tempeff);
    myerr.push_back(temperr);
    SF *= tempeff;
  }
  if(SF==0) return 0;//don't compute error
  for(unsigned int i = 0; i<tightlep.size();++i) error += pow(SF/myeff[i]*myerr[i],2);
  error = sqrt(error);
  return SF;
}

float getTriggerWeightandError(float &error, vector<int> tightlep, vector<int> looselep){
  error = 0;
  if(tightlep.size()+looselep.size()<2) return 0.;
  vector<int> lep = tightlep;
  lep.insert(    lep.end(),    looselep.begin(), looselep.end() );
  sort(lep.begin(),lep.end(),sortPt);
  float eff1(0.), eff2(0.);
  float err1(0.), err2(0.);
  if(abs(lep_pdgId()[lep[0] ])==11) eff1 = lepsf_EGammaVVV_EGammaVVVEleLead(lep_p4()[lep[0] ].Pt(), lep_etaSC()[lep[0] ]);
  else                              eff1 = lepsf_MuTightVVV_MuTightVVVMu17( lep_p4()[lep[0] ].Pt(), lep_p4()[   lep[0] ].Eta());
  if(abs(lep_pdgId()[lep[1] ])==11) eff2 = lepsf_EGammaVVV_EGammaVVVEle12(  lep_p4()[lep[1] ].Pt(), lep_etaSC()[lep[1] ]);
  else                              eff2 = lepsf_MuTightVVV_MuTightVVVMu8(  lep_p4()[lep[1] ].Pt(), lep_p4()[   lep[1] ].Eta());
  if(abs(lep_pdgId()[lep[0] ])==11) err1 = lepsf_EGammaVVV_EGammaVVVEleLead_unc(lep_p4()[lep[0] ].Pt(), lep_etaSC()[lep[0] ]);
  else                              err1 = lepsf_MuTightVVV_MuTightVVVMu17_unc( lep_p4()[lep[0] ].Pt(), lep_p4()[   lep[0] ].Eta());
  if(abs(lep_pdgId()[lep[1] ])==11) err2 = lepsf_EGammaVVV_EGammaVVVEle12_unc(  lep_p4()[lep[1] ].Pt(), lep_etaSC()[lep[1] ]);
  else                              err2 = lepsf_MuTightVVV_MuTightVVVMu8_unc(  lep_p4()[lep[1] ].Pt(), lep_p4()[   lep[1] ].Eta());
  error = sqrt(pow(err1*eff2,2)+pow(eff1*err2,2));
  return eff1*eff2;
}

float getlepFakeRateandError(float &error, int index, bool data, bool conecorr, int version){
  error = 0;
  if(index<0) return 0;
  if(index>=(int)lep_pdgId().size()) return 0;
  float correction = 1.;
  if(conecorr) correction = 1.+coneCorrPt(index, version);
  float conept = correction*lep_p4()[index].Pt();
  if(data){
    if(abs(lep_pdgId()[index])==11) {
      error = fakerate_el_data_unc(lep_p4()[index].Eta(),conept,version);
      return  fakerate_el_data(    lep_p4()[index].Eta(),conept,version);
    }
    else {
      error = fakerate_mu_data_unc(lep_p4()[index].Eta(),conept,version);
      return  fakerate_mu_data(    lep_p4()[index].Eta(),conept,version);
    }
  }
  else {
    if(abs(lep_pdgId()[index])==11) {
      error = fakerate_el_qcd_unc(lep_p4()[index].Eta(),conept,version);
      return  fakerate_el_qcd(    lep_p4()[index].Eta(),conept,version);
    }
    else {
      error = fakerate_mu_qcd_unc(lep_p4()[index].Eta(),conept,version);
      return  fakerate_mu_qcd(    lep_p4()[index].Eta(),conept,version);
    }
  }
}

float getlepFRWeightandError(float &error, int index, bool data, bool conecorr, bool addclosureunc, int version){
  error = 0;
  if(index<0) return 0;
  if(index>=(int)lep_pdgId().size()) return 0;
  float err = 0;
  float FR = getlepFakeRateandError(err, index, data, conecorr, version);
  if(FR>=1.){
    error = 0;
    return 0;
  }
  error = err/pow(1.-FR,2);
  if(addclosureunc){
    if(abs(lep_pdgId()[index])==11) error = sqrt(pow(error,2)+pow(0.17*FR/(1.-FR),2));
    else                            error = sqrt(pow(error,2)+pow(0.31*FR/(1.-FR),2));
  }

  return FR/(1.-FR);
}

float getlepFRClosureError(int index, bool data, bool conecorr, int version){
  if(index<0) return 0;
  if(index>=(int)lep_pdgId().size()) return 0;
  float err = 0;
  float FR = getlepFakeRateandError(err, index, data, conecorr, version);
  if(abs(lep_pdgId()[index])==11) return 0.17*FR/(1.-FR);
  else                            return 0.31*FR/(1.-FR);
}

float getPUWeight(int variation){
  return puWeight(nTrueInt(),variation);
}

float getPUWeightAndError(float &errorup, float &errordown){
  errorup   = getPUWeight( 1);
  errordown = getPUWeight(-1);
  return      getPUWeight();
}

bool addeventtocheck(vector<myevt> &eventvector, unsigned int runnumber, unsigned int lumisection, long long eventnumber){
  myevt event;
  event.run = runnumber;
  event.ls  = lumisection;
  event.evt = eventnumber;
  eventvector.push_back(event);
  return true;
}

bool checkthisevent(vector<myevt> eventvector, unsigned int runnumber, unsigned int lumisection, long long eventnumber){
  for(unsigned int i = 0; i<eventvector.size();++i){
    if(eventvector[i].run!=runnumber  ) continue;
    if(eventvector[i].ls !=lumisection) continue;
    if(eventvector[i].evt!=eventnumber) continue;
    return true;
  }
  return false;
}
void addeventtolist(std::ostringstream *&stream){
  *stream << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
}
void addeventtolist(int SS, std::ostringstream *&streamEE, std::ostringstream *&streamEM, std::ostringstream *&streamMM){
  if(SS==0) *streamEE << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
  if(SS==1) *streamEM << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
  if(SS==2) *streamMM << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
}
void addeventtolist(int SRSS, int SR3l, std::ostringstream *&streamEE, std::ostringstream *&streamEM, std::ostringstream *&streamMM, std::ostringstream *&stream0SFOS, std::ostringstream *&stream1SFOS, std::ostringstream *&stream2SFOS){
  if(SRSS==0) *streamEE    << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
  if(SRSS==1) *streamEM    << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
  if(SRSS==2) *streamMM    << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
  if(SR3l==0) *stream0SFOS << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
  if(SR3l==1) *stream1SFOS << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
  if(SR3l==2) *stream2SFOS << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
}
void storeeventlist(string output, string sample, std::ostringstream *&stream){
  ofstream log    ((string(output+"_"+sample+".log")).c_str(), ios::trunc);
  log    << stream   ->str();
  delete stream;
}
void storeeventlist(string output, string sample, bool isSS, std::ostringstream *&streamEE, std::ostringstream *&streamEM, std::ostringstream *streamMM){
  if(isSS){
    ofstream logEE   ((string(output+   "EE_"+sample+".log")).c_str(), ios::trunc);
    logEE   << streamEE->str();
    ofstream logEM   ((string(output+   "EM_"+sample+".log")).c_str(), ios::trunc);
    logEM   << streamEM->str();
    ofstream logMM   ((string(output+   "MM_"+sample+".log")).c_str(), ios::trunc);
    logMM   << streamMM->str();
  }
  else{
    ofstream log0SFOS((string(output+"0SFOS_"+sample+".log")).c_str(), ios::trunc);
    log0SFOS<< streamEE->str();
    ofstream log1SFOS((string(output+"1SFOS_"+sample+".log")).c_str(), ios::trunc);
    log1SFOS<< streamEM->str();
    ofstream log2SFOS((string(output+"2SFOS_"+sample+".log")).c_str(), ios::trunc);
    log2SFOS<< streamMM->str();
  }
  delete streamEE;
  delete streamEM;
  delete streamMM;

}
void storeeventlist(string output, string sample, std::ostringstream *&streamEE, std::ostringstream *&streamEM, std::ostringstream *&streamMM, std::ostringstream *&stream0SFOS, std::ostringstream *&stream1SFOS, std::ostringstream *&stream2SFOS){
  ofstream logEE   ((string(output+   "EE_"+sample+".log")).c_str(), ios::trunc);
  logEE   << streamEE   ->str();
  ofstream logEM   ((string(output+   "EM_"+sample+".log")).c_str(), ios::trunc);
  logEM   << streamEM   ->str();
  ofstream logMM   ((string(output+   "MM_"+sample+".log")).c_str(), ios::trunc);
  logMM   << streamMM   ->str();
  ofstream log0SFOS((string(output+"0SFOS_"+sample+".log")).c_str(), ios::trunc);
  log0SFOS<< stream0SFOS->str();
  ofstream log1SFOS((string(output+"1SFOS_"+sample+".log")).c_str(), ios::trunc);
  log1SFOS<< stream1SFOS->str();
  ofstream log2SFOS((string(output+"2SFOS_"+sample+".log")).c_str(), ios::trunc);
  log2SFOS<< stream2SFOS->str();
  delete streamEE;
  delete streamEM;
  delete streamMM;
  delete stream0SFOS;
  delete stream1SFOS;
  delete stream2SFOS;
}

int calcNSFOS(std::vector<int> tightlep)
{
  bool OS01 = (lep_pdgId()[tightlep[0] ]*lep_pdgId()[tightlep[1] ]<0); bool SF01 = (abs(lep_pdgId()[tightlep[0] ])==abs(lep_pdgId()[tightlep[1] ]));
  bool OS02 = (lep_pdgId()[tightlep[0] ]*lep_pdgId()[tightlep[2] ]<0); bool SF02 = (abs(lep_pdgId()[tightlep[0] ])==abs(lep_pdgId()[tightlep[2] ]));
  bool OS12 = (lep_pdgId()[tightlep[1] ]*lep_pdgId()[tightlep[2] ]<0); bool SF12 = (abs(lep_pdgId()[tightlep[1] ])==abs(lep_pdgId()[tightlep[2] ]));
  int SFOScounter = 0;
  if(OS01&&SF01) ++SFOScounter;
  if(OS02&&SF02) ++SFOScounter;
  if(OS12&&SF12) ++SFOScounter;
  if(OS01&&SF01&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()<20.) SFOScounter = -1;
  if(OS02&&SF02&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()<20.) SFOScounter = -1;
  if(OS12&&SF12&&(lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()<20.) SFOScounter = -1;
  if(abs(lep_charge()[tightlep[0] ]+lep_charge()[tightlep[1] ]+lep_charge()[tightlep[2] ])==3) SFOScounter = -1;
  return SFOScounter;
}

float get0SFOSMll(std::vector<int> lepidx)
{
    if (lepidx.size() != 3)
    { return -999; }
    int pdgid0 = tas::lep_pdgId()[lepidx[0]];
    int pdgid1 = tas::lep_pdgId()[lepidx[1]];
    int pdgid2 = tas::lep_pdgId()[lepidx[2]];
    if (pdgid0 == pdgid1)
        return (tas::lep_p4()[lepidx[0]]
                + tas::lep_p4()[lepidx[1]]).mass();
    else if (pdgid0 == pdgid2)
        return (tas::lep_p4()[lepidx[0]]
                + tas::lep_p4()[lepidx[2]]).mass();
    else if (pdgid1 == pdgid2)
        return (tas::lep_p4()[lepidx[1]]
                + tas::lep_p4()[lepidx[2]]).mass();
    std::cout <<
              "Warning: Shouldn't be here if function call are at the right places."
              << std::endl;
    return -999;
}

float get0SFOSMee(std::vector<int> lepidx)
{
    if (lepidx.size() != 3)
    { return -999; }
    int pdgid0 = tas::lep_pdgId()[lepidx[0]];
    int pdgid1 = tas::lep_pdgId()[lepidx[1]];
    int pdgid2 = tas::lep_pdgId()[lepidx[2]];
    if (pdgid0 == pdgid1 && abs(pdgid0) == 11)
        return (tas::lep_p4()[lepidx[0]]
                + tas::lep_p4()[lepidx[1]]).mass();
    else if (pdgid0 == pdgid2 && abs(pdgid0) == 11)
        return (tas::lep_p4()[lepidx[0]]
                + tas::lep_p4()[lepidx[2]]).mass();
    else if (pdgid1 == pdgid2 && abs(pdgid1) == 11)
        return (tas::lep_p4()[lepidx[1]]
                + tas::lep_p4()[lepidx[2]]).mass();
    return -999;
}

float get1SFOSMll(std::vector<int> lepidx)
{
    if (lepidx.size() != 3)
    { return -999; }
    int pdgid0 = tas::lep_pdgId()[lepidx[0]];
    int pdgid1 = tas::lep_pdgId()[lepidx[1]];
    int pdgid2 = tas::lep_pdgId()[lepidx[2]];
    if (pdgid0 == -pdgid1)
        return (tas::lep_p4()[lepidx[0]]
                + tas::lep_p4()[lepidx[1]]).mass();
    else if (pdgid0 == -pdgid2)
        return (tas::lep_p4()[lepidx[0]]
                + tas::lep_p4()[lepidx[2]]).mass();
    else if (pdgid1 == -pdgid2)
        return (tas::lep_p4()[lepidx[1]]
                + tas::lep_p4()[lepidx[2]]).mass();
    std::cout <<
              "Warning: Shouldn't be here if function call are at the right places."
              << std::endl;
    return -999;
}

float get2SFOSMll0(std::vector<int> lepidx)
{
    if (lepidx.size() != 3)
    { return -999; }
    int pdgid0 = tas::lep_pdgId()[lepidx[0]];
    int pdgid1 = tas::lep_pdgId()[lepidx[1]];
    int pdgid2 = tas::lep_pdgId()[lepidx[2]];
    if (pdgid0 == -pdgid1)
        return (tas::lep_p4()[lepidx[0]]
                + tas::lep_p4()[lepidx[1]]).mass();
    else if (pdgid0 == -pdgid2)
        return (tas::lep_p4()[lepidx[0]]
                + tas::lep_p4()[lepidx[2]]).mass();
    else if (pdgid1 == -pdgid2)
        return (tas::lep_p4()[lepidx[1]]
                + tas::lep_p4()[lepidx[2]]).mass();
    std::cout <<
              "Warning: Shouldn't be here if function call are at the right places."
              << std::endl;
    return -999;
}

float get2SFOSMll1(std::vector<int> lepidx)
{
    if (lepidx.size() != 3)
    { return -999; }
    int pdgid0 = tas::lep_pdgId()[lepidx[0]];
    int pdgid1 = tas::lep_pdgId()[lepidx[1]];
    int pdgid2 = tas::lep_pdgId()[lepidx[2]];
    if (pdgid2 == -pdgid1)
        return (tas::lep_p4()[lepidx[1]]
                + tas::lep_p4()[lepidx[2]]).mass();
    else if (pdgid0 == -pdgid2)
        return (tas::lep_p4()[lepidx[0]]
                + tas::lep_p4()[lepidx[2]]).mass();
    else if (pdgid1 == -pdgid0)
        return (tas::lep_p4()[lepidx[0]]
                + tas::lep_p4()[lepidx[1]]).mass();
    std::cout <<
              "Warning: Shouldn't be here if function call are at the right places."
              << std::endl;
    return -999;
}

// returns if the given igen'th gen particle is a W to be considered.
bool isW(int igen)
{
    if (!(abs(tas::genPart_pdgId()[igen]) == 24)) return false;
    if (!(abs(tas::genPart_status()[igen]) == 62)) return false;
    return true;
}

// returns if the given igen'th gen particle is a desecendent of W. (The fermions from W boson)
bool isParentW(int igen)
{
    if (!(abs(tas::genPart_motherId()[igen]) == 24)) return false;
    if (!(abs(tas::genPart_pdgId()[igen]) != 24)) return false;
    return true;
}

// returns if the given igen'th gen particle is a desecendent of W. (The fermions from W boson)
bool isParentTau(int igen)
{
    if (!(abs(tas::genPart_motherId()[igen]) == 15)) return false;
    if (!(abs(tas::genPart_pdgId()[igen]) != 15)) return false;
    return true;
}

// returns if the given igen'th gen particle is a lepton.
bool isLepton(int igen)
{
    if (!(abs(tas::genPart_pdgId()[igen]) == 11 ||
          abs(tas::genPart_pdgId()[igen]) == 13 ||
          abs(tas::genPart_pdgId()[igen]) == 15 )) return false;
    return true;
}

// returns if the given igen'th gen particle is a lepton.
bool isQuark(int igen)
{
    if (!(abs(tas::genPart_pdgId()[igen]) <= 5 )) return false;
    return true;
}

// returns if the given igen'th gen particle is a lepton.
bool isLightLepton(int igen)
{
    if (!(abs(tas::genPart_pdgId()[igen]) == 11 ||
          abs(tas::genPart_pdgId()[igen]) == 13 )) return false;
    return true;
}

// returns if the given igen'th gen particle is a lepton.
bool isTau(int igen)
{
    if (!(abs(tas::genPart_pdgId()[igen]) == 15)) return false;
    return true;
}

// Loops over gen particles and returns W and W's desecendents
std::tuple<vIdx, vIdx> getGenIndices()
{
//    vIdx iws;
//    vIdx ifs;
//    vIdx ils;
    vIdx iqs;
//    vIdx ills;
//    vIdx ilts;
//    vIdx its;
    vIdx ifls;
    for (unsigned int igen = 0; igen < tas::genPart_pdgId().size(); ++igen)
    {
//        /* W bosons */ if (isW(igen)) iws.push_back(igen);
//        /* Fermions */ if (isParentW(igen)) ifs.push_back(igen);
//        /* Leptons  */ if (isParentW(igen) && isLepton(igen)) ils.push_back(igen);
        /* Quarks   */ if (isParentW(igen) && isQuark(igen)) iqs.push_back(igen);
//        /* e, muons */ if (isParentW(igen) && isLightLepton(igen)) ills.push_back(igen);
//        /* taus     */ if (isParentW(igen) && isTau(igen)) its.push_back(igen);
//        /* e,m <- t */ if (isParentTau(igen) && isLepton(igen)) ilts.push_back(igen);
        /* e, muons */ if ((isParentTau(igen) && isLepton(igen)) || (isParentW(igen) && isLightLepton(igen))) ifls.push_back(igen);
    }
//    return std::make_tuple(iws, ifs, ils, iqs, ills, ilts, its, ifls);
    return std::make_tuple(iqs, ifls);
}

