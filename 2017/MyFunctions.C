#include "MyFunctions.h"

using namespace std;
using namespace tas;

int gentype_v2(unsigned lep1_index,unsigned lep2_index, int lep3_index){
  bool gammafake = false;
  bool jetfake   = false;
  unsigned int ngenlep = ngenLepFromTau()+ ngenLep();
  unsigned int nW(0), nZ(0);
  bool lep1_real = lep_motherIdSS().at(lep1_index) > 0;
  bool lep2_real = lep_motherIdSS().at(lep2_index) > 0;
  bool lep3_real = false;
  if(lep3_index>0) lep3_real = lep_motherIdSS().at(lep3_index) > 0;
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

bool getleptonindices(vector<int> &iSS, vector<int> &i3l, vector<int> &iaSS, vector<int> &ia3l, vector<int> &vSS, vector<int> &v3l, vector<int> &vaSS, vector<int> &va3l){
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
    if(lep_pass_VVV_cutbased_tight_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015&&lep_isTriggerSafe_v1()[i]) {
      if(abs(lep_pdgId()[i])==11&&lep_lostHits()[i]==0&&lep_tightCharge()[i]==2){
	if(fabs(lep_etaSC()[i])<=1.479&&lep_relIso03EAv2()[i]<0.0588){
	  if(lep_p4()[i].Pt()>20) { i3l.push_back(i); is3l = true; }
	  if(lep_p4()[i].Pt()>30) { iSS.push_back(i); isSS = true; }
	} else if(fabs(lep_etaSC()[i]) >1.479&&lep_relIso03EAv2()[i]<0.0571){
	  if(lep_p4()[i].Pt()>20) { i3l.push_back(i); is3l = true; }
	  if(lep_p4()[i].Pt()>30) { iSS.push_back(i); isSS = true; }
	}
      } else if(abs(lep_pdgId()[i])==13&&lep_relIso03EAv2()[i]<0.06){
	if(  lep_p4()[i].Pt()>20) { i3l.push_back(i); is3l = true; }
	if(  lep_p4()[i].Pt()>30) { iSS.push_back(i); isSS = true; }
      }
    }
    //loose ID
    if(lep_pass_VVV_cutbased_fo_noiso()[i]   &&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015&&lep_isTriggerSafe_v1()[i]) {
      if(abs(lep_pdgId()[i])==11&&lep_lostHits()[i]==0&&lep_tightCharge()[i]==2){
	if(fabs(lep_etaSC()[i])<=1.479&&lep_relIso03EAv2()[i]<0.2){
	  if(!is3l&&lep_p4()[i].Pt()>20) { ia3l.push_back(i); isa3l = true; }
	  if(!isSS&&lep_p4()[i].Pt()>30) { iaSS.push_back(i); isaSS = true; }
	} else if(fabs(lep_etaSC()[i]) >1.479&&lep_relIso03EAv2()[i]<0.2){
	  if(!is3l&&lep_p4()[i].Pt()>20) { ia3l.push_back(i); isa3l = true; }
	  if(!isSS&&lep_p4()[i].Pt()>30) { iaSS.push_back(i); isaSS = true; }
	}
      } else if(abs(lep_pdgId()[i])==13&&lep_relIso03EAv2()[i]<0.4){
	if(  !is3l&&lep_p4()[i].Pt()>20) { ia3l.push_back(i); isa3l = true; }
	if(  !isSS&&lep_p4()[i].Pt()>30) { iaSS.push_back(i); isaSS = true; }
      }
    }
    //vetoID
    if(lep_pass_VVV_cutbased_veto_noiso()[i] &&fabs(lep_p4()[i].Eta())<2.4&&lep_relIso03EAv2()[i]<=0.4) {
      if(!isSS &&          lep_p4()[i].Pt()>10) vSS  .push_back(i);
      if(!is3l &&          lep_p4()[i].Pt()>10) v3l  .push_back(i);
      if(!isSS && !isaSS &&lep_p4()[i].Pt()>10) vaSS .push_back(i);
      if(!is3l && !isa3l &&lep_p4()[i].Pt()>10) va3l .push_back(i);
    }
  }
  return true;
}

float coneCorrPt(int lepi){
  float coneptcorr = 0;
  float relIso = lep_relIso03EAv2().at(lepi);
  if(   abs(lep_pdgId().at(lepi)) == 11) {
    if(fabs(lep_eta().at(lepi)) <= 1.479) coneptcorr = std::max( 0., relIso - 0.0588 );
    else                                  coneptcorr = std::max( 0., relIso - 0.0571 );
  }
  if(   abs(lep_pdgId().at(lepi)) == 13)  coneptcorr = std::max( 0., relIso - 0.06   );

  return coneptcorr;
}
float loadFR(float &FRSSerr,int index, TH2D *hMuFR, TH2D *hElFR, float muFRptmin, float muFRptmax, float muFRetamin, float muFRetamax, float elFRptmin, float elFRptmax, float elFRetamin, float elFRetamax, bool conecorrected){
  float FR = 0.;
  float FRerr = 0;
  if(index<0) return -1.;
  if((unsigned int)index>=lep_pdgId().size()) return -1;
  float correction = 1.;
  if(conecorrected) correction = (1.+coneCorrPt(index));
  int bin;
  float binx,biny;
  if(abs(lep_pdgId()[index])==11){
    if(!(hElFR==0)){
      if(elFRetamin<0) biny = std::max(elFRetamin,std::min(elFRetamax,     lep_etaSC()[index]) );
      else             biny = std::max(elFRetamin,std::min(elFRetamax,fabs(lep_etaSC()[index])));
      binx                  = std::max(elFRptmin, std::min(elFRptmax, correction*lep_p4()[index].Pt()));
      bin   = hElFR->FindBin(binx, biny);
      FR    = hElFR->GetBinContent(bin);
      FRerr = hElFR->GetBinError(  bin);
    }
  }
  if(abs(lep_pdgId()[index])==13){
    if(!(hMuFR==0)){
      if(muFRetamin<0) biny = std::max(muFRetamin,std::min(muFRetamax,     lep_p4()[index].Eta()) );
      else             biny = std::max(muFRetamin,std::min(muFRetamax,fabs(lep_p4()[index].Eta())));
      binx                  = std::max(muFRptmin, std::min(muFRptmax, correction*lep_p4()[index].Pt()));
      bin   = hMuFR->FindBin(binx, biny);
      FR    = hMuFR->GetBinContent(bin);
      FRerr = hMuFR->GetBinError(  bin);
    }
  }
  FRSSerr = FRerr;
 
  return FR;
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
  if(index2>0)  MT2 = mT(lep_p4()[index1],MET);
  if(MT1>MT2)  return MT1;
  return MT2;
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

bool passFilters(){
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
  if(splitVH(filename))                               return "WHtoWWW";//be careful not to double-count (for Hannsjoerg only)
  if(filename.find("www_2l_")         !=string::npos) return "WWW";    //be careful not to double-count (for Hannsjoerg only)
  if(filename.find("www_incl_amcnlo_")!=string::npos) return "WWWv2";  //be careful not to double-count (for Hannsjoerg only)
  if(filename.find("data_")           !=string::npos) return "Data";  //be careful not to double-count (for Hannsjoerg only)
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
  //if((filename.find("wjets_")!=string::npos||filename.find("dy_")!=string::npos||filename.find("ttbar_")!=string::npos||filename.find("ww_")!=string::npos||filename.find("wz_")!=string::npos)&&(process==5)) return true;
  //if((filename.find("wgjets_")!=string::npos||filename.find("wgstar_")!=string::npos||filename.find("zgamma_")!=string::npos||filename.find("ttg_")!=string::npos||filename.find("wwg_")!=string::npos ||filename.find("wzg_")!=string::npos)&&(process!=5)) return true;
  if((filename.find("wjets_" )!=string::npos||filename.find("dy_"    )!=string::npos||filename.find("ttbar_" )!=string::npos||filename.find("ww_" )!=string::npos||filename.find("wz_" )!=string::npos)&&(process)) return true;
  if((filename.find("wgjets_")!=string::npos||filename.find("wgstar_")!=string::npos||filename.find("zgamma_")!=string::npos||filename.find("ttg_")!=string::npos||filename.find("wwg_")!=string::npos ||filename.find("wzg_")!=string::npos)&&(!process)) return true;
  return false;
}

bool passJetSSstate(bool preselect, int nj, int nb, float Mjj, float MjjL, float Detajj, bool is3lCR, int jec){
  int numj = nj;
  int numb = nb;
  int x;
  if(nb<0||nj<0) getalljetnumbers(x,numj,numb, jec);
  if(numj<2)    return false;
  if(numb>0)    return false;
  if(preselect) return true;
  float Massjj  = Mjj;
  float MassjjL = MjjL;
  float Deta    = Detajj;
  if(Mjj<0||MjjL<0||Detajj<0) getMjjAndDeta(Massjj,MassjjL,Deta,jec);
  if(fabs(Deta)>1.5)       return false;
  if(fabs(MassjjL)>400.)   return false;
  if(is3lCR)               return true;
  if(fabs(Massjj-80.)>20.) return false;
  return true;
}
      
int isSRSS(vector<int> tightlep, vector<int> vetolep,  bool preselect, float maxMT, int nj, int nb, float Mjj, float MjjL, float Detajj, LorentzVector MET, int jec){
  if(!passJetSSstate(preselect, nj, nb, Mjj, MjjL, Detajj, false, jec)) return -1;
  if(tightlep.size()!=2)                                                return -1;
  if(vetolep.size()!=0)                                                 return -1;
  if(nisoTrack_mt2_cleaned_VVV_cutbased_veto()!=0)                      return -1;
  if((lep_pdgId()[tightlep[0] ]*lep_pdgId()[tightlep[1] ])<0)           return -1;
  bool ee = false; bool em = false; bool mm = false;
  if(((lep_pdgId()[tightlep[0] ])*(lep_pdgId()[tightlep[1] ]))==121&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()-MZ)>10.) ee = true;
  if(((lep_pdgId()[tightlep[0] ])*(lep_pdgId()[tightlep[1] ]))==143)  em = true;
  if(((lep_pdgId()[tightlep[0] ])*(lep_pdgId()[tightlep[1] ]))==169)  mm = true;
  if(preselect){
    if(ee) return 0;
    if(em) return 1;
    if(mm) return 2;
    return -1;
  }
  LorentzVector met = LorentzVector(0,0,0,0);
  if(MET.Pt()>0)     met = MET;
  else if(jec==1)    met.SetPxPyPzE(met_T1CHS_miniAOD_CORE_up_pt()*std::cos(met_T1CHS_miniAOD_CORE_up_phi()),met_T1CHS_miniAOD_CORE_up_pt()*std::sin(met_T1CHS_miniAOD_CORE_up_phi()),0,met_T1CHS_miniAOD_CORE_up_pt());
  else if(jec==(-1)) met.SetPxPyPzE(met_T1CHS_miniAOD_CORE_dn_pt()*std::cos(met_T1CHS_miniAOD_CORE_dn_phi()),met_T1CHS_miniAOD_CORE_dn_pt()*std::sin(met_T1CHS_miniAOD_CORE_dn_phi()),0,met_T1CHS_miniAOD_CORE_dn_pt());
  else               met.SetPxPyPzE(met_pt()*std::cos(met_phi()),met_pt()*std::sin(met_phi()),0,met_pt());
  if(ee&&met.Pt()>40.&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()>40.) return 0;
  if(mm&&              (lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()>40.) return 2;
  if(em&&met.Pt()>40.&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()>30.) {
    if(      maxMT>90.) return 1;
    else if(maxMT<0){
      float tempmaxMT = calcMTmax(tightlep[0],tightlep[1],met);
      if(tempmaxMT>90.) return 1;
    }
  }
  return -1;
}

int isARSS(vector<int> tightlep, vector<int> looselep, vector<int> vetolep,  bool preselect, float maxMT, int nj, int nb, float Mjj, float MjjL, float Detajj, LorentzVector MET, int jec){
  if(!passJetSSstate(preselect, nj, nb, Mjj, MjjL, Detajj, false, jec)) return -1;
  if(tightlep.size()!=1)                                                return -1;
  if(looselep.size()!=1)                                                return -1;
  if(vetolep.size()!=0)                                                 return -1;
  if(nisoTrack_mt2_cleaned_VVV_cutbased_veto()!=0)                      return -1;
  if((lep_pdgId()[tightlep[0] ]*lep_pdgId()[looselep[0] ])<0)           return -1;
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
  LorentzVector met = LorentzVector(0,0,0,0);
  if(MET.Pt()>0)     met = MET;
  else if(jec==1)    met.SetPxPyPzE(met_T1CHS_miniAOD_CORE_up_pt()*std::cos(met_T1CHS_miniAOD_CORE_up_phi()),met_T1CHS_miniAOD_CORE_up_pt()*std::sin(met_T1CHS_miniAOD_CORE_up_phi()),0,met_T1CHS_miniAOD_CORE_up_pt());
  else if(jec==(-1)) met.SetPxPyPzE(met_T1CHS_miniAOD_CORE_dn_pt()*std::cos(met_T1CHS_miniAOD_CORE_dn_phi()),met_T1CHS_miniAOD_CORE_dn_pt()*std::sin(met_T1CHS_miniAOD_CORE_dn_phi()),0,met_T1CHS_miniAOD_CORE_dn_pt());
  else               met.SetPxPyPzE(met_pt()*std::cos(met_phi()),met_pt()*std::sin(met_phi()),0,met_pt());
  if(ee&&met.Pt()>40.&&(lep_p4()[tightlep[0] ]+lep_p4()[looselep[0] ]).M()>40.) return 0;
  if(mm&&              (lep_p4()[tightlep[0] ]+lep_p4()[looselep[0] ]).M()>40.) return 2;
  if(em&&met.Pt()>40.&&(lep_p4()[tightlep[0] ]+lep_p4()[looselep[0] ]).M()>30.) {
    if(      maxMT>90.) return 1;
    else if(maxMT<0){
      float tempmaxMT = calcMTmax(tightlep[0],looselep[0],met);
      if(tempmaxMT>90.) return 1;
    }
  }
  return -1;
}

int isCRSS(vector<int> tightlep, vector<int> selectlep, vector<int> vetolep, bool preselect, float maxMT, int nj, int nb, float Mjj, float MjjL, float Detajj, LorentzVector MET, int jec, bool noZ){
  if(!passJetSSstate(preselect, nj, nb, Mjj, MjjL, Detajj, true, jec)) return -1;
  if(tightlep.size()  <2)                                              return -1;
  if(selectlep.size()!=3)                                              return -1;
  if(vetolep.size()!=0)                                                return -1;
  if(nisoTrack_mt2_cleaned_VVV_cutbased_veto()!=0)                     return -1;
  //SS is just tighter pt, so can loop over i3l
  bool isSS = false;
  int lep3 = -1; int SS1 = -1; int SS2 = -1;
  if((lep_pdgId()[selectlep[0] ]*lep_pdgId()[selectlep[1] ])>0) { lep3 = selectlep[2]; SS1 = selectlep[0]; SS2 = selectlep[1]; isSS = true; }
  if((lep_pdgId()[selectlep[0] ]*lep_pdgId()[selectlep[2] ])>0) { lep3 = selectlep[1]; SS1 = selectlep[0]; SS2 = selectlep[2]; isSS = true; }
  if((lep_pdgId()[selectlep[1] ]*lep_pdgId()[selectlep[2] ])>0) { lep3 = selectlep[0]; SS1 = selectlep[1]; SS2 = selectlep[2]; isSS = true; }
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
  else if(jec==1)    met.SetPxPyPzE(met_T1CHS_miniAOD_CORE_up_pt()*std::cos(met_T1CHS_miniAOD_CORE_up_phi()),met_T1CHS_miniAOD_CORE_up_pt()*std::sin(met_T1CHS_miniAOD_CORE_up_phi()),0,met_T1CHS_miniAOD_CORE_up_pt());
  else if(jec==(-1)) met.SetPxPyPzE(met_T1CHS_miniAOD_CORE_dn_pt()*std::cos(met_T1CHS_miniAOD_CORE_dn_phi()),met_T1CHS_miniAOD_CORE_dn_pt()*std::sin(met_T1CHS_miniAOD_CORE_dn_phi()),0,met_T1CHS_miniAOD_CORE_dn_pt());
  else               met.SetPxPyPzE(met_pt()*std::cos(met_phi()),met_pt()*std::sin(met_phi()),0,met_pt());
  if(ee&&met.Pt()>40.&&(lep_p4()[selectlep[0] ]+lep_p4()[selectlep[1] ]).M()>40.) return 0;
  if(mm&&              (lep_p4()[selectlep[0] ]+lep_p4()[selectlep[1] ]).M()>40.) return 2;
  if(em&&met.Pt()>40.&&(lep_p4()[selectlep[0] ]+lep_p4()[selectlep[1] ]).M()>30.) {
    if(      maxMT>90.) return 1;
    else if(maxMT<0){
      float tempmaxMT = calcMTmax(SS1,SS2,met);
      if(tempmaxMT>90.) return 1;
    }
  }
  return -1;
}
		      
int isSR3l(vector<int> tightlep, bool preselect, int nj, int nb, LorentzVector MET, int jec){
  int numj = nj;
  int numb = nb;
  int x;
  if(numj<0||numb<0) getalljetnumbers(numj,x,numb, jec);
  if(numj >1) return -1;
  if(numb!=0) return -1;
  if(tightlep.size()!=3) return -1;
  float DPhilllMET = -1;
  LorentzVector met = LorentzVector(0,0,0,0);
  if(!preselect){
    if(MET.Pt()>0)     met = MET;
    else if(jec==1)    met.SetPxPyPzE(met_T1CHS_miniAOD_CORE_up_pt()*std::cos(met_T1CHS_miniAOD_CORE_up_phi()),met_T1CHS_miniAOD_CORE_up_pt()*std::sin(met_T1CHS_miniAOD_CORE_up_phi()),0,met_T1CHS_miniAOD_CORE_up_pt());
    else if(jec==(-1)) met.SetPxPyPzE(met_T1CHS_miniAOD_CORE_dn_pt()*std::cos(met_T1CHS_miniAOD_CORE_dn_phi()),met_T1CHS_miniAOD_CORE_dn_pt()*std::sin(met_T1CHS_miniAOD_CORE_dn_phi()),0,met_T1CHS_miniAOD_CORE_dn_pt());
    else               met.SetPxPyPzE(met_pt()*std::cos(met_phi()),met_pt()*std::sin(met_phi()),0,met_pt());
    if(fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()-MZ)<10.) return -1;
    if(     (lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).Pt()   <60.) return -1;
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
  }
  if(SFOScounter==1){
    pass1 = true;
    if(!preselect&&met.Pt()<45) { pass1=false; }
    if(OS01&&SF01&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()>55.&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()<110.) { pass1 = false; }
    if(OS02&&SF02&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()>55.&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()<110.) { pass1 = false; }
    if(OS12&&SF12&&(lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()>55.&&(lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()<110.) { pass1 = false; }
  }
  if(SFOScounter==2){
    pass2 = true;
    if(!preselect&&met.Pt()<55) { pass2=false; }
    if(OS01&&SF01&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()-MZ)<20.) { pass2 = false; }
    if(OS02&&SF02&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()-MZ)<20.) { pass2 = false; }
    if(OS12&&SF12&&fabs((lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()-MZ)<20.) { pass2 = false; }
  }
  if(DPhilllMET>=2.7&&pass0 ) return 0;
  if(DPhilllMET>=2.5&&pass1 ) return 1;
  if(DPhilllMET>=2.5&&pass2 ) return 2;
  return -1;
}

int isAR3l(vector<int> tightlep, vector<int> looselep, bool preselect, int nj, int nb, LorentzVector MET, int jec){
  int numj = nj;
  int numb = nb;
  int x;
  if(numj<0||numb<0) getalljetnumbers(numj,x,numb, jec);
  if(numj >1) return -1;
  if(numb!=0) return -1;
  if(tightlep.size()!=2) return -1;
  if(looselep.size()!=1) return -1;
  float DPhilllMET = -1;
  LorentzVector met = LorentzVector(0,0,0,0);
  if(!preselect){
    if(MET.Pt()>0)     met = MET;
    else if(jec==1)    met.SetPxPyPzE(met_T1CHS_miniAOD_CORE_up_pt()*std::cos(met_T1CHS_miniAOD_CORE_up_phi()),met_T1CHS_miniAOD_CORE_up_pt()*std::sin(met_T1CHS_miniAOD_CORE_up_phi()),0,met_T1CHS_miniAOD_CORE_up_pt());
    else if(jec==(-1)) met.SetPxPyPzE(met_T1CHS_miniAOD_CORE_dn_pt()*std::cos(met_T1CHS_miniAOD_CORE_dn_phi()),met_T1CHS_miniAOD_CORE_dn_pt()*std::sin(met_T1CHS_miniAOD_CORE_dn_phi()),0,met_T1CHS_miniAOD_CORE_dn_pt());
    else               met.SetPxPyPzE(met_pt()*std::cos(met_phi()),met_pt()*std::sin(met_phi()),0,met_pt());
    if(fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]+lep_p4()[looselep[0] ]).M()-MZ)<10.) return -1;
    if(     (lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]+lep_p4()[looselep[0] ]).Pt()   <60.) return -1;
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
  }
  if(SFOScounter==1){
    pass1 = true;
    if(!preselect&&met.Pt()<45) { pass1=false; }
    if(OS01&&SF01&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()>55.&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()<110.) { pass1 = false; }
    if(OS02&&SF02&&(lep_p4()[tightlep[0] ]+lep_p4()[looselep[0] ]).M()>55.&&(lep_p4()[tightlep[0] ]+lep_p4()[looselep[0] ]).M()<110.) { pass1 = false; }
    if(OS12&&SF12&&(lep_p4()[tightlep[1] ]+lep_p4()[looselep[0] ]).M()>55.&&(lep_p4()[tightlep[1] ]+lep_p4()[looselep[0] ]).M()<110.) { pass1 = false; }
  }
  if(SFOScounter==2){
    pass2 = true;
    if(!preselect&&met.Pt()<55) { pass2=false; }
    if(OS01&&SF01&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()-MZ)<20.) { pass2 = false; }
    if(OS02&&SF02&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[looselep[0] ]).M()-MZ)<20.) { pass2 = false; }
    if(OS12&&SF12&&fabs((lep_p4()[tightlep[1] ]+lep_p4()[looselep[0] ]).M()-MZ)<20.) { pass2 = false; }
  }
  if(DPhilllMET>=2.7&&pass0 ) return 0;
  if(DPhilllMET>=2.5&&pass1 ) return 1;
  if(DPhilllMET>=2.5&&pass2 ) return 2;
  return -1;
}

int isCR3l(vector<int> tightlep, bool preselect, int nj, int nb, LorentzVector MET, int jec){
  int numj = nj;
  int numb = nb;
  int x;
  if(numj<0||numb<0) getalljetnumbers(numj,x,numb, jec);
  if(numj >1) return -1;
  if(numb!=0) return -1;
  if(tightlep.size()!=3) return -1;
  float DPhilllMET = -1;
  LorentzVector met = LorentzVector(0,0,0,0);
  if(!preselect){
    if(MET.Pt()>0)     met = MET;
    else if(jec==1)    met.SetPxPyPzE(met_T1CHS_miniAOD_CORE_up_pt()*std::cos(met_T1CHS_miniAOD_CORE_up_phi()),met_T1CHS_miniAOD_CORE_up_pt()*std::sin(met_T1CHS_miniAOD_CORE_up_phi()),0,met_T1CHS_miniAOD_CORE_up_pt());
    else if(jec==(-1)) met.SetPxPyPzE(met_T1CHS_miniAOD_CORE_dn_pt()*std::cos(met_T1CHS_miniAOD_CORE_dn_phi()),met_T1CHS_miniAOD_CORE_dn_pt()*std::sin(met_T1CHS_miniAOD_CORE_dn_phi()),0,met_T1CHS_miniAOD_CORE_dn_pt());
    else               met.SetPxPyPzE(met_pt()*std::cos(met_phi()),met_pt()*std::sin(met_phi()),0,met_pt());
    if(fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()-MZ)<10.) return -1;
    if(     (lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).Pt()   <60.) return -1;
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
    pass1X = true;
    if(!preselect&&met.Pt()<45) { pass1X=false; }
    bool atleastoneSFOSZ = false;
    if(OS01&&SF01&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()>55.&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()<110.) { atleastoneSFOSZ = true; }
    if(OS02&&SF02&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()>55.&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()<110.) { atleastoneSFOSZ = true; }
    if(OS12&&SF12&&(lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()>55.&&(lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()<110.) { atleastoneSFOSZ = true; }
    if(!atleastoneSFOSZ) pass1X = false;
  }
  if(SFOScounter==2){
    pass2X = true;
    if(!preselect&&met.Pt()<55) { pass2X=false; }
    bool atleastoneSFOSZ = false;
    if(OS01&&SF01&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()-MZ)<20.) { atleastoneSFOSZ = true; }
    if(OS02&&SF02&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()-MZ)<20.) { atleastoneSFOSZ = true; }
    if(OS12&&SF12&&fabs((lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()-MZ)<20.) { atleastoneSFOSZ = true; }
    if(!atleastoneSFOSZ) pass2X = false;
  }
  if(DPhilllMET>=2.7&&pass0X ) return 0;
  if(DPhilllMET>=2.5&&pass1X ) return 1;
  if(DPhilllMET>=2.5&&pass2X ) return 2;
  return -1;
}

bool checkbothSRCR3l(int &isSR3l, int &isCR3l, vector<int> tightlep, bool preselect, int nj, int nb, LorentzVector MET, int jec){
  isSR3l = -1; isCR3l = -1;//assume no match
  int numj = nj;
  int numb = nb;
  int x;
  if(numj<0||numb<0) getalljetnumbers(numj,x,numb, jec);
  if(numj >1) return true;
  if(numb!=0) return true;
  if(tightlep.size()!=3) return true;
  float DPhilllMET = -1;
  LorentzVector met = LorentzVector(0,0,0,0);
  if(!preselect){
    if(MET.Pt()>0)     met = MET;
    else if(jec==1)    met.SetPxPyPzE(met_T1CHS_miniAOD_CORE_up_pt()*std::cos(met_T1CHS_miniAOD_CORE_up_phi()),met_T1CHS_miniAOD_CORE_up_pt()*std::sin(met_T1CHS_miniAOD_CORE_up_phi()),0,met_T1CHS_miniAOD_CORE_up_pt());
    else if(jec==(-1)) met.SetPxPyPzE(met_T1CHS_miniAOD_CORE_dn_pt()*std::cos(met_T1CHS_miniAOD_CORE_dn_phi()),met_T1CHS_miniAOD_CORE_dn_pt()*std::sin(met_T1CHS_miniAOD_CORE_dn_phi()),0,met_T1CHS_miniAOD_CORE_dn_pt());
    else               met.SetPxPyPzE(met_pt()*std::cos(met_phi()),met_pt()*std::sin(met_phi()),0,met_pt());
    if(fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()-MZ)<10.) return true;
    if(     (lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).Pt()   <60.) return true;
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
  }
  if(SFOScounter==1){
    pass1 = true; pass1X = true;
    if(!preselect&&met.Pt()<45) { pass1=false; pass1X=false; }
    bool atleastoneSFOSZ = false;
    if(OS01&&SF01&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()>55.&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
    if(OS02&&SF02&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()>55.&&(lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
    if(OS12&&SF12&&(lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()>55.&&(lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
    if(!atleastoneSFOSZ) pass1X = false;
  }
  if(SFOScounter==2){
    pass2 = true; pass2X = true;
    if(!preselect&&met.Pt()<55) { pass2=false; pass2X=false; }
    bool atleastoneSFOSZ = false;
    if(OS01&&SF01&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[1] ]).M()-MZ)<20.) { pass2 = false; atleastoneSFOSZ = true; }
    if(OS02&&SF02&&fabs((lep_p4()[tightlep[0] ]+lep_p4()[tightlep[2] ]).M()-MZ)<20.) { pass2 = false; atleastoneSFOSZ = true; }
    if(OS12&&SF12&&fabs((lep_p4()[tightlep[1] ]+lep_p4()[tightlep[2] ]).M()-MZ)<20.) { pass2 = false; atleastoneSFOSZ = true; }
    if(!atleastoneSFOSZ) pass2X = false;
  }
  //SR
  if(DPhilllMET>=2.7&&pass0 ) isSR3l = 0;
  if(DPhilllMET>=2.5&&pass1 ) isSR3l = 1;
  if(DPhilllMET>=2.5&&pass2 ) isSR3l = 2;
  //CR
  if(DPhilllMET>=2.7&&pass0X) isCR3l = 0;
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

map<string, TH1D*> bookhistograms(string samplename, vector<string> histonames, vector<int> hbins, vector<float> hlow, vector<float> hup, TDirectory *rootdir, bool splitWW){
  map<string, TH1D*> histos;
  if(histonames.size()!=hbins.size()) return histos;
  if(histonames.size()!=hlow.size()) return histos;
  if(histonames.size()!=hup.size()) return histos;
  for(unsigned int i = 0; i<histonames.size(); ++i){
    string mapname = histonames[i];
    if(splitWW){
      if(samplename=="WW"){
	mapname = histonames[i] + "_WWRest";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	mapname = histonames[i] + "_otherWWSS";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	mapname = histonames[i] + "_WWVBS";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	mapname = histonames[i] + "_WWDPS";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      }
      if(samplename.find("ttV")!=string::npos){
	mapname = histonames[i] + "_ttVSS";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	mapname = histonames[i] + "_OtherttV";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	mapname = histonames[i] + "_ttV3l";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      }
    }
    if(samplename.find("Other")!=string::npos){
      mapname = histonames[i] + "_Other";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_WHtoWWW";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    }
    else if(samplename.find("Background")!=string::npos){
      mapname = histonames[i] + "_trueSS";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_chargeflips";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_SSLL";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_fakes";//jetfakes
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
    } else {
      mapname = histonames[i] + "_"+samplename;
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    }
  }
  for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
    h->second->Sumw2(); h->second->SetDirectory(rootdir);
  }
  return histos;
}

bool fillSRhisto(map<string, TH1D*> histos, string histoname, string sn, string sn2, int SRSS, int SR3l, float weight, float weight3l){
  float weight2 = weight3l;
  if(weight3l<-1e6) weight2 = weight;
  if(SRSS==0) histos[histoname+"_"+sn ]->Fill(0.,weight);
  if(SRSS==1) histos[histoname+"_"+sn ]->Fill(1.,weight);
  if(SRSS==2) histos[histoname+"_"+sn ]->Fill(2.,weight);
  if(SR3l==0) histos[histoname+"_"+sn2]->Fill(3.,weight2);
  if(SR3l==1) histos[histoname+"_"+sn2]->Fill(4.,weight2);
  if(SR3l==2) histos[histoname+"_"+sn2]->Fill(5.,weight2);
  if(SRSS>=0||SR3l>=0) return true;
  else return false;
}

bool SaveHistosToFile(string filename, map<string, TH1D*> histos, bool addunderflow, bool addoverflow){
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
  
  TFile *f = new TFile(filename.c_str(),"UPDATE");
  f->cd();
  for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
    h->second->Write(h->first.c_str(),TObject::kOverwrite);
  }
  f->Close();
  cout << "Saved histos in " << f->GetName() << endl;
  return true;
}

bool fileexists(string filename){
  ifstream f(filename.c_str());
  return f.good();
}

int loadlepSFfile(TFile *&f, TH2F *&hMuID, TH2F *&hMutrigger, TH2F *&hElID, TH2F *&hEltrigger, string filename, string muIDname, string mutrigname, string elIDname, string eltrigname){

  if(filename=="") return -1;
  if(!fileexists(filename)) return -1;
  f = new TFile(filename.c_str(),"read");
  int id = 0;
  if(f->GetListOfKeys()->Contains(muIDname.c_str()  ))  hMuID      = (TH2F*)f->Get(muIDname.c_str());
  else id += 10;
  //if(f->GetListOfKeys()->Contains(mutrigname.c_str())) hMutrigger = (TH2F*)f->Get(mutrigname.c_str());
  //else id += 100;
  if(f->GetListOfKeys()->Contains(elIDname.c_str()  ))  hElID      = (TH2F*)f->Get(elIDname.c_str());
  else id +=1000;
  if(f->GetListOfKeys()->Contains(eltrigname.c_str()))  hEltrigger = (TH2F*)f->Get(eltrigname.c_str());
  else id +=10000;
  return id;

}
bool loadFakeRates(TFile *&f, TH2D *&hMuFR, TH2D *&hElFR, string filename, string muname, string elname){

  if(filename=="") {
    cout << "Filename for FR not specified - return false" << endl;
    return false;
  }
  if(!fileexists(filename)){
    cout << "FR File " << filename << " does not exist - couldn't load FR" << endl;
    return false;
  }
  f = new TFile(filename.c_str(),"read");
  if(!(f->GetListOfKeys()->Contains(muname.c_str()))){
    cout << "FRHistogram " << muname << " does not exist - couldn't load FR" << endl;
    delete f;
    return false;
  }
  if(!(f->GetListOfKeys()->Contains(elname.c_str()))){
    cout << "FRHistogram " << elname << " does not exist - couldn't load FR" << endl;
    delete f;
    return false;
  }
  hMuFR = (TH2D*)f->Get(muname.c_str());
  hElFR = (TH2D*)f->Get(elname.c_str());
  return true;

}
bool deleteFiles(bool SF, TFile *&fSF){
  if(SF) { fSF->Close(); delete fSF; }
  return true;
}

float getlepSFandError(float &error, int index, TH2F *hMuID, TH2F *hMutrigger, TH2F *hElID, TH2F *hEltrigger, float muIDptmin,float muIDptmax, float muIDetamin, float muIDetamax, float muTrptmin,float muTrptmax, float muTretamin, float muTretamax, float elIDptmin,float elIDptmax, float elIDetamin, float elIDetamax, float elTrptmin,float elTrptmax, float elTretamin, float elTretamax){
  float SF1 = 1.;
  float SF1err2 = 0;
  float SF2 = 1.;
  float SF2err2 = 0;
  if(index<0) return -1.;
  if((unsigned int)index>=lep_pdgId().size()) return -1;
  int bin;
  float binx,biny;
  if(abs(lep_pdgId()[index])==11){
    if(!(hElID==0)){
      if(elIDetamin<0)   binx = std::max(elIDetamin,std::min(elIDetamax,lep_etaSC()[index]));
      else               binx = std::max(elIDetamin,std::min(elIDetamax,fabs(lep_etaSC()[index])));
      biny = std::max(elIDptmin,std::min(elIDptmax, lep_p4()[index].Pt()));
      bin = hElID->FindBin(binx, biny);
      SF1 = hElID->GetBinContent(bin);
      SF1err2 = pow(hElID->GetBinError(bin),2);
    }
    if(!(hEltrigger==0)){
      if(elTretamin<0)   binx = std::max(elTretamin,std::min(elTretamax,lep_etaSC()[index]));
      else               binx = std::max(elTretamin,std::min(elTretamax,fabs(lep_etaSC()[index])));
      biny = std::max(elTrptmin,std::min(elTrptmax, lep_p4()[index].Pt()));
      bin = hEltrigger->FindBin(binx, biny);
      SF2 = hEltrigger->GetBinContent(bin);
      SF2err2 = pow(hEltrigger->GetBinError(bin),2);
    }
  }
  if(abs(lep_pdgId()[index])==13){
    if(!(hMuID==0)){
      if(muIDetamin<0)   binx = std::max(muIDetamin,std::min(muIDetamax,lep_p4()[index].Eta()));
      else               binx = std::max(muIDetamin,std::min(muIDetamax,fabs(lep_p4()[index].Eta())));
      biny = std::max(muIDptmin,std::min(muIDptmax, lep_p4()[index].Pt()));
      bin = hMuID->FindBin(binx, biny);
      SF1 = hMuID->GetBinContent(bin);
      SF1err2 = pow(hMuID->GetBinError(bin),2);
    }
    //if(!(hMutrigger==0)){
    //  if(muTretamin<0)   binx = std::max(muTretamin,std::min(muTretamax,lep_p4()[index].Eta()));
    //  else               binx = std::max(muTretamin,std::min(muTretamax,fabs(lep_p4()[index].Eta())));
    //  biny = std::max(muTrptmin,std::min(muTrptmax, lep_p4()[index].Pt()));
    //  bin = hMutrigger->FindBin(binx, biny);
    //  SF2 = hMutrigger->GetBinContent(bin);
    //  SF2err2 = pow(hMutrigger->GetBinError(bin),2);
    //}
  }
  error = sqrt(pow(SF1,2)*SF2err2+pow(SF2,2)*SF1err2);
 
  return SF1*SF2;
}

float getlepSFWeightandError(float &error, vector<int> tightlep, vector<int> looselep, TH2F *hMuID, TH2F *hMutrigger, TH2F *hElID, TH2F *hEltrigger, float muIDptmin,float muIDptmax, float muIDetamin, float muIDetamax, float muTrptmin,float muTrptmax, float muTretamin, float muTretamax, float elIDptmin,float elIDptmax, float elIDetamin, float elIDetamax, float elTrptmin,float elTrptmax, float elTretamin, float elTretamax){
  vector <float> SF;
  vector <float> SFerr;
  float tempSF, temperr;
  for(unsigned int i = 0; i<tightlep.size(); ++i){
    tempSF = getlepSFandError(temperr, tightlep[i], hMuID, hMutrigger, hElID, hEltrigger, muIDptmin, muIDptmax, muIDetamin, muIDetamax, muTrptmin, muTrptmax, muTretamin, muTretamax, elIDptmin, elIDptmax, elIDetamin, elIDetamax, elTrptmin, elTrptmax, elTretamin, elTretamax);
    SF.push_back(tempSF);
    SFerr.push_back(temperr);
  }
  for(unsigned int i = 0; i<looselep.size(); ++i){
    tempSF = getlepSFandError(temperr, looselep[i], hMuID, hMutrigger, hElID, hEltrigger, muIDptmin, muIDptmax, muIDetamin, muIDetamax, muTrptmin, muTrptmax, muTretamin, muTretamax, elIDptmin, elIDptmax, elIDetamin, elIDetamax, elTrptmin, elTrptmax, elTretamin, elTretamax);
    SF.push_back(tempSF);
    SFerr.push_back(temperr);
  }
  return getlepSFWeightandError(error,SF,SFerr);
}

float getlepSFWeightandError(float &error, vector<float> efftight, vector<float> errtight, vector<float> effloose, vector<float> errloose){

  vector<float> SF = efftight;
  SF.insert(    SF.end(),    effloose.begin(), effloose.end() );
  vector<float> SFerr = errtight;
  SFerr.insert( SFerr.end(), errloose.begin(), errloose.end() );
  float SFsum = 1;
  float errsum2 = 0;
  for(unsigned int i = 0; i<SF   .size();++i) { if(SF[i]>0) SFsum *=SF[i]; }
  for(unsigned int i = 0; i<SFerr.size();++i) { if(SF[i]>0) errsum2 +=pow(SFsum/SF[i]*SFerr[i],2); }
  error = sqrt(errsum2);
  return SFsum;

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
