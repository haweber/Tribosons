float Zbi_(float sig, float bg, float bgrelunc=0.3){
    
    double bgunc = bgrelunc*bg;
    double tau = bg/pow(bgunc,2);//bgunc is absolute
    double n_on = sig+bg;//total yield in SR = sig + bg
    double n_off = tau*bg;
    double P_Bi = TMath::BetaIncomplete(1./(1.+tau),n_on,n_off+1);
    double Z_Bi = sqrt(2.)*TMath::ErfInverse(1.-2.*P_Bi);
    return Z_Bi;
}


void GetSensitivityPerBin(){

    bool data = false;
    bool twosig = false;
    bool addsig = true;
    bool addlostlepto3l = true;
    map<string, TH1D*> h;
    map<string, float> tZbi;
    map<string, float> tSB;
    vector<string> samples;
    vector<string> histonames;

    bool SFOS = true; bool SS = true;
    TFile *f = TFile::Open("InvestigateSR.root");
    
    
    bool mainSR  = true;
    bool softjet = true;
    bool onejet  = true;
    bool twojet  = true;
    bool Mjjside = true;
    
    if(SS){
        if(mainSR){
            histonames.push_back("SR_addone_SSee");
            histonames.push_back("SR_addone_dropMET_SSee");
            histonames.push_back("SR_addone_dropMETMjj_SSee");
            histonames.push_back("SR_addone_dropMETMll_SSee");
            histonames.push_back("SR_addone_dropMETMjjMll_SSee");
            histonames.push_back("SR_addone_dropMjj_SSee");
            histonames.push_back("SR_addone_dropMll_SSee");
            histonames.push_back("SR_addone_dropDeta_SSee");
            histonames.push_back("SR_addone_dropDetaMjj_SSee");
            histonames.push_back("SR_addone_SSem");
            histonames.push_back("SR_addone_dropMET_SSem");
            histonames.push_back("SR_addone_dropMETMjj_SSem");
            histonames.push_back("SR_addone_dropMETMll_SSem");
            histonames.push_back("SR_addone_dropMjj_SSem");
            histonames.push_back("SR_addone_dropMll_SSem");
            histonames.push_back("SR_addone_dropDeta_SSem");
            histonames.push_back("SR_addone_dropDetaMjj_SSem");
            histonames.push_back("SR_addone_dropMTmax_SSem");
            histonames.push_back("SR_addone_dropMTmaxMll_SSem");
            histonames.push_back("SR_addone_dropMTmaxMjj_SSem");
            histonames.push_back("SR_addone_SSmm");
            histonames.push_back("SR_addone_dropMET_SSmm");
            histonames.push_back("SR_addone_dropMETMjj_SSmm");
            histonames.push_back("SR_addone_dropMETMll_SSmm");
            histonames.push_back("SR_addone_dropMjj_SSmm");
            histonames.push_back("SR_addone_dropMll_SSmm");
            histonames.push_back("SR_addone_dropDeta_SSmm");
            histonames.push_back("SR_addone_dropDetaMjj_SSmm");
        } if(onejet){
            histonames.push_back("SRonejet_addone_SSee");
            histonames.push_back("SRonejet_addone_dropMET_SSee");
            histonames.push_back("SRonejet_addone_dropMETMll_SSee");
            histonames.push_back("SRonejet_addone_dropMll_SSee");
            histonames.push_back("SRonejet_addone_SSem");
            histonames.push_back("SRonejet_addone_dropMET_SSem");
            histonames.push_back("SRonejet_addone_dropMETMll_SSem");
            histonames.push_back("SRonejet_addone_dropMll_SSem");
            histonames.push_back("SRonejet_addone_dropMTmax_SSem");
            histonames.push_back("SRonejet_addone_dropMTmaxMll_SSem");
            histonames.push_back("SRonejet_addone_SSmm");
            histonames.push_back("SRonejet_addone_dropMET_SSmm");
            histonames.push_back("SRonejet_addone_dropMETMll_SSmm");
            histonames.push_back("SRonejet_addone_dropMll_SSmm");
        } if(softjet){
            histonames.push_back("SRsoftjet_addone_SSee");
            histonames.push_back("SRsoftjet_addone_dropMET_SSee");
            histonames.push_back("SRsoftjet_addone_dropMETMjj_SSee");
            histonames.push_back("SRsoftjet_addone_dropMETMll_SSee");
            histonames.push_back("SRsoftjet_addone_dropMETMjjMll_SSee");
            histonames.push_back("SRsoftjet_addone_dropMjj_SSee");
            histonames.push_back("SRsoftjet_addone_dropMll_SSee");
            histonames.push_back("SRsoftjet_addone_dropDeta_SSee");
            histonames.push_back("SRsoftjet_addone_dropDetaMjj_SSee");
            histonames.push_back("SRsoftjet_addone_SSem");
            histonames.push_back("SRsoftjet_addone_dropMET_SSem");
            histonames.push_back("SRsoftjet_addone_dropMETMjj_SSem");
            histonames.push_back("SRsoftjet_addone_dropMETMll_SSem");
            histonames.push_back("SRsoftjet_addone_dropMjj_SSem");
            histonames.push_back("SRsoftjet_addone_dropMll_SSem");
            histonames.push_back("SRsoftjet_addone_dropDeta_SSem");
            histonames.push_back("SRsoftjet_addone_dropDetaMjj_SSem");
            histonames.push_back("SRsoftjet_addone_dropMTmax_SSem");
            histonames.push_back("SRsoftjet_addone_dropMTmaxMll_SSem");
            histonames.push_back("SRsoftjet_addone_dropMTmaxMjj_SSem");
            histonames.push_back("SRsoftjet_addone_SSmm");
            histonames.push_back("SRsoftjet_addone_dropMET_SSmm");
            histonames.push_back("SRsoftjet_addone_dropMETMjj_SSmm");
            histonames.push_back("SRsoftjet_addone_dropMETMll_SSmm");
            histonames.push_back("SRsoftjet_addone_dropMjj_SSmm");
            histonames.push_back("SRsoftjet_addone_dropMll_SSmm");
            histonames.push_back("SRsoftjet_addone_dropDeta_SSmm");
            histonames.push_back("SRsoftjet_addone_dropDetaMjj_SSmm");
        } if(Mjjside){
            histonames.push_back("SRsideMjj_addone_SSee");
            histonames.push_back("SRsideMjj_addone_dropMET_SSee");
            histonames.push_back("SRsideMjj_addone_dropMETMll_SSee");
            histonames.push_back("SRsideMjj_addone_dropMll_SSee");
            histonames.push_back("SRsideMjj_addone_dropDeta_SSee");
            histonames.push_back("SRsideMjj_addone_SSem");
            histonames.push_back("SRsideMjj_addone_dropMET_SSem");
            histonames.push_back("SRsideMjj_addone_dropMETMll_SSem");
            histonames.push_back("SRsideMjj_addone_dropMll_SSem");
            histonames.push_back("SRsideMjj_addone_dropDeta_SSem");
            histonames.push_back("SRsideMjj_addone_dropMTmax_SSem");
            histonames.push_back("SRsideMjj_addone_dropMTmaxMll_SSem");
            histonames.push_back("SRsideMjj_addone_SSmm");
            histonames.push_back("SRsideMjj_addone_dropMET_SSmm");
            histonames.push_back("SRsideMjj_addone_dropMETMll_SSmm");
            histonames.push_back("SRsideMjj_addone_dropMll_SSmm");
            histonames.push_back("SRsideMjj_addone_dropDeta_SSmm");
        }
    }
    if(SFOS){
        if(mainSR){
            histonames.push_back("SR_addone_0SFOS");
            histonames.push_back("SR_addone_dropMET_0SFOS");
            histonames.push_back("SR_addone_dropDPhilllMET_0SFOS");
            histonames.push_back("SR_addone_1SFOS");
            histonames.push_back("SR_addone_dropMET_1SFOS");
            histonames.push_back("SR_addone_droppTlll_1SFOS");
            histonames.push_back("SR_addone_dropDPhilllMET_1SFOS");
            histonames.push_back("SR_addone_droppTlllDPhilllMET_1SFOS");
            histonames.push_back("SR_addone_dropMT3rd_1SFOS");
            histonames.push_back("SR_addone_2SFOS");
            histonames.push_back("SR_addone_dropMET_2SFOS");
            histonames.push_back("SR_addone_droppTlll_2SFOS");
            histonames.push_back("SR_addone_dropDPhilllMET_2SFOS");
            histonames.push_back("SR_addone_droppTlllDPhilllMET_2SFOS");
        } if(twojet){
            histonames.push_back("SRsoftjet_addone_0SFOS");
            histonames.push_back("SRsoftjet_addone_dropMET_0SFOS");
            histonames.push_back("SRsoftjet_addone_dropDPhilllMET_0SFOS");
            histonames.push_back("SRsoftjet_addone_1SFOS");
            histonames.push_back("SRsoftjet_addone_dropMET_1SFOS");
            histonames.push_back("SRsoftjet_addone_droppTlll_1SFOS");
            histonames.push_back("SRsoftjet_addone_dropDPhilllMET_1SFOS");
            histonames.push_back("SRsoftjet_addone_droppTlllDPhilllMET_1SFOS");
            histonames.push_back("SRsoftjet_addone_dropMT3rd_1SFOS");
            histonames.push_back("SRsoftjet_addone_2SFOS");
            histonames.push_back("SRsoftjet_addone_dropMET_2SFOS");
            histonames.push_back("SRsoftjet_addone_droppTlll_2SFOS");
            histonames.push_back("SRsoftjet_addone_dropDPhilllMET_2SFOS");
            histonames.push_back("SRsoftjet_addone_droppTlllDPhilllMET_2SFOS");
        }
    }



    vector<int> is3lSS;
    samples.push_back("WWW");
    if(twosig) { samples.push_back("WHtoWWW");      }
    samples.push_back("photonfakes");
    samples.push_back("chargeflips");
    samples.push_back("fakes");
    if(!addlostlepto3l) { samples.push_back("3lLL");  }
    samples.push_back("SSLL");
    samples.push_back("true3L");
    samples.push_back("trueWWW");
    samples.push_back("trueSS");
    samples.push_back("bg");
    int nsig = 1;
    if(twosig) ++nsig;
    for(unsigned int n = 0; n<histonames.size();++n){
      string bgname = histonames[n] + "_bg";
      string stackname = histonames[n];
      string mapname = "";
      if(data) {
	mapname = histonames[n] + "_Data";
	string  mapnameX = histonames[n] + "_X";
	h[mapnameX ]=(TH1D*)f->Get(mapname.c_str());
	h[mapnameX]->SetName(mapnameX.c_str());
	h[mapname] = new TH1D(mapname.c_str(),"",h[mapnameX]->GetNbinsX(), h[mapnameX]->GetXaxis()->GetBinLowEdge(1), h[mapnameX]->GetXaxis()->GetBinLowEdge(h[mapnameX]->GetNbinsX()+1));
	h[mapname]->SetBinErrorOption(TH1::kPoisson);
	for(int i = 1; i<=h[mapname]->GetNbinsX(); ++i){
	  //cout << h[mapnameX]->GetBinContent(i) << endl;
	  for(int n = 1; n<=h[mapnameX]->GetBinContent(i); ++n){
	    h[mapname]->Fill(h[mapnameX]->GetXaxis()->GetBinCenter(i),1);
	  }
	}
	h[mapname ]=(TH1D*)f->Get(mapname.c_str());
      }
      for(unsigned int s = 0; s<samples.size()-1; ++s){
	mapname = histonames[n] + "_" + samples[s];
	cout << mapname << endl;
	h[mapname ]=(TH1D*)f->Get(mapname.c_str());
	if(addsig&&s==0){
	  string addname = histonames[n] + "_WHtoWWW";
	  h[addname ]=(TH1D*)f->Get(addname.c_str());
	  h[mapname ]->Add(h[addname ],1.);
	}
    if(addlostlepto3l&&samples[s].find("true3L")!=string::npos){
              string mapname2 =histonames[n] + "_3lLL";
              h[mapname2 ]=(TH1D*)f->Get(mapname2.c_str());
              h[mapname]->Add(h[mapname2],1.);
          }
	for(int b = 0; b<h[mapname]->GetNbinsX();++b){ if(h[mapname]->GetBinContent(b)<0){h[mapname]->SetBinContent(b,0); h[mapname]->SetBinError(b,0); } }
	if(s==nsig) h[bgname] = (TH1D*)h[mapname]->Clone(bgname.c_str());
	else if(s>nsig) {
	  h[bgname]->Add(h[mapname],1.);
	}

      }//samples

    }//histograms
    
    for(unsigned int n = 0; n<histonames.size();++n){
      string stackname = histonames[n];
      string bgname = stackname + "_bg";
        int SR = -1;
        //            ee    em     mm      0      1      2
        //LL         45%    35%    25%    50%    20%    25%
        //fake       80%    70%    50%    60%    80%    80%
        //irr        35%    25%    20%    35%    50%    100%
        //Qflip     100%   100%   100%   100%   100%   100%
        //photon    100%   100%   100%   100%   100%   100%
        //signal     45%    25%    15%    35%    25%    40%
        //stat       45%    15%    10%    30%  17.5%    35%
        
        if(histonames[n].find("SSee" )!=string::npos) SR = 1;
        if(histonames[n].find("SSem" )!=string::npos) SR = 2;
        if(histonames[n].find("SSmm" )!=string::npos) SR = 3;
        if(histonames[n].find("0SFOS")!=string::npos) SR = 4;
        if(histonames[n].find("1SFOS")!=string::npos) SR = 5;
        if(histonames[n].find("2SFOS")!=string::npos) SR = 6;
        if(SR<0) continue;
        float sigsyst(0), sigstat(0), LLsyst(0), fakesyst(0), Qflsyst(0), photsyst(0), irrsyst(0);
        
        if(SR==1){ sigsyst = 0.45; sigstat = 0.45; LLsyst = 0.45; fakesyst = 0.80; irrsyst = 0.35; Qflsyst = 1.0; photsyst = 1.0; }
        if(SR==2){ sigsyst = 0.25; sigstat = 0.15; LLsyst = 0.35; fakesyst = 0.70; irrsyst = 0.35; Qflsyst = 1.0; photsyst = 1.0; }
        if(SR==3){ sigsyst = 0.15; sigstat = 0.10; LLsyst = 0.25; fakesyst = 0.50; irrsyst = 0.30; Qflsyst = 1.0; photsyst = 1.0; }
        if(SR==4){ sigsyst = 0.35; sigstat = 0.30; LLsyst = 0.50; fakesyst = 0.60; irrsyst = 0.35; Qflsyst = 1.0; photsyst = 1.0; }
        if(SR==5){ sigsyst = 0.25; sigstat = 0.18; LLsyst = 0.20; fakesyst = 0.80; irrsyst = 0.50; Qflsyst = 1.0; photsyst = 1.0; }
        if(SR==6){ sigsyst = 0.40; sigstat = 0.35; LLsyst = 0.25; fakesyst = 0.80; irrsyst = 1.00; Qflsyst = 1.0; photsyst = 1.0; }
        
        for(int i = 1; i<=h[bgname]->GetNbinsX(); ++i){
            if(i==1){
                tZbi[stackname] = -99.;
                tSB[ stackname] = -99.;
            }
            if(h[stackname+"_WWW"]->GetBinContent(i) == 0) continue;
            float sig = h[stackname+"_WWW"]->GetBinContent(i);
            float sigunc = sqrt(pow(h[stackname+"_WWW"]->GetBinError(i),2)+pow(sigsyst*sig,2)-pow(sigstat*sig,2));
            float LL = h[stackname+"_SSLL"]->GetBinContent(i)+h[stackname+"_true3L"]->GetBinContent(i);
            float LLunc = std::max(LLsyst*LL,(float)sqrt(pow(h[stackname+"_SSLL"]->GetBinError(i),2)+pow(h[stackname+"_true3L"]->GetBinError(i),2)));
            float fake = h[stackname+"_fakes"]->GetBinContent(i);
            float fakeunc = std::max(fakesyst*fake,(float)h[stackname+"_fakes"]->GetBinError(i));
            //float fake = std::min(h[stackname+"_fakes"]->GetBinContent(i),3.);
            //float fakeunc = std::max(fakesyst*fake,(float)std::min(h[stackname+"_fakes"]->GetBinError(i),3.));
            float irr = h[stackname+"_trueSS"]->GetBinContent(i)+h[stackname+"_trueWWW"]->GetBinContent(i);
            float irrunc = std::max(irrsyst*irr,(float)sqrt(pow(h[stackname+"_trueSS"]->GetBinError(i),2)+pow(h[stackname+"_trueWWW"]->GetBinError(i),2)));
            float Qfl = h[stackname+"_chargeflips"]->GetBinContent(i);
            float Qflunc = std::max(Qflsyst*Qfl,(float)h[stackname+"_chargeflips"]->GetBinError(i));
            float phot = h[stackname+"_photonfakes"]->GetBinContent(i);
            float photunc = std::max(photsyst*phot,(float)h[stackname+"_photonfakes"]->GetBinError(i));
            float bg = LL + fake + irr + Qfl + phot;
            float bgunc = sqrt(pow(LLunc,2)+pow(fakeunc,2)+pow(irrunc,2)+pow(Qflunc,2)+pow(photunc,2));
            if(i==1){
                tZbi[stackname] = Zbi_(sig,bg,bgunc/bg);
                tSB[ stackname] = sig/sqrt(bg+bgunc*bgunc+sigunc*sigunc);
            }
            string s1 = "";
            if(stackname.find("SR_")       !=string::npos) s1 = "SR_addone_";
            if(stackname.find("SRonejet_") !=string::npos) s1 = "SRonejet_addone_";
            if(stackname.find("SRsoftjet_")!=string::npos) s1 = "SRsoftjet_addone_";
            if(stackname.find("SRsideMjj_")!=string::npos) s1 = "SRsideMjj_addone_";
            string s2 = "";
            if(stackname.find("_SSee") !=string::npos) s2 = "SSee";
            if(stackname.find("_SSem") !=string::npos) s2 = "SSem";
            if(stackname.find("_SSmm") !=string::npos) s2 = "SSmm";
            if(stackname.find("_0SFOS")!=string::npos) s2 = "0SFOS";
            if(stackname.find("_1SFOS")!=string::npos) s2 = "1SFOS";
            if(stackname.find("_2SFOS")!=string::npos) s2 = "2SFOS";
            string ss = s1+s2;
            if((ss==stackname&&i==1) || Zbi_(sig,bg,bgunc/bg)>tZbi[stackname] || sig/sqrt(bg+bgunc*bgunc+sigunc*sigunc)>tSB[stackname]){
                //cout << ss << " " << stackname << " " << i << " " << (ss==stackname) << " " << Zbi_(sig,bg,bgunc/bg) << " " <<tZbi[stackname] << " " << (Zbi_(sig,bg,bgunc/bg)>tZbi[stackname]) << " " << sig/sqrt(bg+bgunc*bgunc+sigunc*sigunc) << " " << tSB[stackname] << " " << (sig/sqrt(bg+bgunc*bgunc+sigunc*sigunc)>SB[stackname]) << endl;
                //cout << stackname << " bin " << int(i-1) << ": signal " << fixed << setprecision(2) << sig << " +/- " << sigunc << ", bg " << bg << " +/- " << bgunc << " --> Z_bi " << Zbi_(sig,bg,bgunc/bg) << ", S/sqrt(B+dS+dB) " << sig/sqrt(bg+bgunc+sigunc) << endl;
                cout << stackname << " bin " << int(i-1) << ": signal " << fixed << setprecision(2) << sig << " +/- " << sigunc << ", bg " << bg << " +/- " << bgunc << " --> Z_bi " << Zbi_(sig,bg,bgunc/bg) << ", S/sqrt(B+dS+dB) " << sig/sqrt(bg+bgunc*bgunc+sigunc*sigunc) << endl;
                cout << "     " << stackname << " bin " << int(i-1)<< ": irr " << fixed << setprecision(2)  << irr << " +/- " << irrunc << ", LL " << LL << " +/- " << LLunc << ", fake " << fake << " +/- " << fakeunc << ", phot " << phot << " +/- " << photunc << ", Qflip " << Qfl << " +/- " << Qflunc << endl;
            }
        }
        cout << endl;
    }
}
