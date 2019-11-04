void GetSRNumbersQuick(TString filename="www_private.root", TString filedir="hists/combineyearsLoose_v5.3.1/master/"){

  TString filetotalname = filedir + TString("/") + filename;
  TFile *f = TFile::Open(filetotalname);
  TH1F *heein  = (TH1F*)f->Get("SRSSeeMjjInFull__yield");
  TH1F *heeout = (TH1F*)f->Get("SRSSeeMjjOutFull__yield");
  TH1F *hemin  = (TH1F*)f->Get("SRSSemMjjInFull__yield");
  TH1F *hemout = (TH1F*)f->Get("SRSSemMjjOutFull__yield");
  TH1F *hmmin  = (TH1F*)f->Get("SRSSmmMjjInFull__yield");
  TH1F *hmmout = (TH1F*)f->Get("SRSSmmMjjOutFull__yield");
  TH1F *hee1J  = (TH1F*)f->Get("SRSS1JeeFull__yield");
  TH1F *hem1J  = (TH1F*)f->Get("SRSS1JemFull__yield");
  TH1F *hmm1J  = (TH1F*)f->Get("SRSS1JmmFull__yield");
  TH1F *h0SFOS = (TH1F*)f->Get("SR0SFOSFull__yield");
  TH1F *h1SFOS = (TH1F*)f->Get("SR1SFOSFull__yield");
  TH1F *h2SFOS = (TH1F*)f->Get("SR2SFOSFull__yield");

  cout << "┬──────────────┬" << endl;
  cout << "│    signal    │" << endl;
  cout << "┼──────────────┼" << endl;
  cout << "│  0.0 ± 0.0   │" << endl;
  cout << "│ " << fixed << setprecision(2) << heein ->GetBinContent(1) << " ± " << heein ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << hemin ->GetBinContent(1) << " ± " << hemin ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << hmmin ->GetBinContent(1) << " ± " << hmmin ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << heeout->GetBinContent(1) << " ± " << heeout->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << hemout->GetBinContent(1) << " ± " << hemout->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << hmmout->GetBinContent(1) << " ± " << hmmout->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << hee1J ->GetBinContent(1) << " ± " << hee1J ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << hem1J ->GetBinContent(1) << " ± " << hem1J ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << hmm1J ->GetBinContent(1) << " ± " << hmm1J ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h0SFOS->GetBinContent(1) << " ± " << h0SFOS->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h1SFOS->GetBinContent(1) << " ± " << h1SFOS->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h2SFOS->GetBinContent(1) << " ± " << h2SFOS->GetBinError(1) << "  │" << endl;
  cout << "│  0.0 ± 0.0   │" << endl;
  cout << "┴──────────────┴" << endl;

}
