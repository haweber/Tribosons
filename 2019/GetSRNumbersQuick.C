void GetSRNumbersQuick(TString filename="www_private.root", TString filedir="hists/combineyearsLoose_v5.3.2/master/"){

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

  TH1F *h2ee    = (TH1F*)f->Get("WZCRSSeeFull__yield");
  TH1F *h2em    = (TH1F*)f->Get("WZCRSSemFull__yield");
  TH1F *h2mm    = (TH1F*)f->Get("WZCRSSmmFull__yield");
  TH1F *h2eein  = (TH1F*)f->Get("WZCRSSeeMjjInFull__yield");
  TH1F *h2eeout = (TH1F*)f->Get("WZCRSSeeMjjOutFull__yield");
  TH1F *h2emin  = (TH1F*)f->Get("WZCRSSemMjjInFull__yield");
  TH1F *h2emout = (TH1F*)f->Get("WZCRSSemMjjOutFull__yield");
  TH1F *h2mmin  = (TH1F*)f->Get("WZCRSSmmMjjInFull__yield");
  TH1F *h2mmout = (TH1F*)f->Get("WZCRSSmmMjjOutFull__yield");
  TH1F *h2ee1J  = (TH1F*)f->Get("WZCRSS1JeeFull__yield");
  TH1F *h2em1J  = (TH1F*)f->Get("WZCRSS1JemFull__yield");
  TH1F *h2mm1J  = (TH1F*)f->Get("WZCRSS1JmmFull__yield");
  TH1F *h20SFOS = (TH1F*)f->Get("WZCR0SFOSFull__yield");
  TH1F *h21SFOS = (TH1F*)f->Get("WZCR1SFOSFull__yield");
  TH1F *h22SFOS = (TH1F*)f->Get("WZCR2SFOSFull__yield");

  TH1F *h3eein  = (TH1F*)f->Get("ARSSeeMjjInFull__yield");
  TH1F *h3eeout = (TH1F*)f->Get("ARSSeeMjjOutFull__yield");
  TH1F *h3emin  = (TH1F*)f->Get("ARSSemMjjInFull__yield");
  TH1F *h3emout = (TH1F*)f->Get("ARSSemMjjOutFull__yield");
  TH1F *h3mmin  = (TH1F*)f->Get("ARSSmmMjjInFull__yield");
  TH1F *h3mmout = (TH1F*)f->Get("ARSSmmMjjOutFull__yield");
  TH1F *h3ee1J  = (TH1F*)f->Get("ARSS1JeeFull__yield");
  TH1F *h3em1J  = (TH1F*)f->Get("ARSS1JemFull__yield");
  TH1F *h3mm1J  = (TH1F*)f->Get("ARSS1JmmFull__yield");
  TH1F *h30SFOS = (TH1F*)f->Get("AR0SFOSFull__yield");
  TH1F *h31SFOS = (TH1F*)f->Get("AR1SFOSFull__yield");
  TH1F *h32SFOS = (TH1F*)f->Get("AR2SFOSFull__yield");

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

  cout << "┬──────────────┬" << endl;
  cout << "│ signal WZ CR │" << endl;
  cout << "┼──────────────┼" << endl;
  cout << "│  0.0 ± 0.0   │" << endl;
  cout << "│ " << fixed << setprecision(2) << h2ee   ->GetBinContent(1) << " ± " << h2ee   ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h2em   ->GetBinContent(1) << " ± " << h2em   ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h2mm   ->GetBinContent(1) << " ± " << h2mm   ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h2ee1J ->GetBinContent(1) << " ± " << h2ee1J ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h2em1J ->GetBinContent(1) << " ± " << h2em1J ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h2mm1J ->GetBinContent(1) << " ± " << h2mm1J ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h20SFOS->GetBinContent(1) << " ± " << h20SFOS->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h21SFOS->GetBinContent(1) << " ± " << h21SFOS->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h22SFOS->GetBinContent(1) << " ± " << h22SFOS->GetBinError(1) << "  │" << endl;
  cout << "│  0.0 ± 0.0   │" << endl;
  cout << "┴──────────────┴" << endl;

  cout << "┬──────────────┬" << endl;
  cout << "│ signal WZ CR* │" << endl;
  cout << "┼──────────────┼" << endl;
  cout << "│  0.0 ± 0.0   │" << endl;
  cout << "│ " << fixed << setprecision(2) << h2eein ->GetBinContent(1) << " ± " << h2eein ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h2emin ->GetBinContent(1) << " ± " << h2emin ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h2mmin ->GetBinContent(1) << " ± " << h2mmin ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h2eeout->GetBinContent(1) << " ± " << h2eeout->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h2emout->GetBinContent(1) << " ± " << h2emout->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h2mmout->GetBinContent(1) << " ± " << h2mmout->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h2ee1J ->GetBinContent(1) << " ± " << h2ee1J ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h2em1J ->GetBinContent(1) << " ± " << h2em1J ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h2mm1J ->GetBinContent(1) << " ± " << h2mm1J ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h20SFOS->GetBinContent(1) << " ± " << h20SFOS->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h21SFOS->GetBinContent(1) << " ± " << h21SFOS->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h22SFOS->GetBinContent(1) << " ± " << h22SFOS->GetBinError(1) << "  │" << endl;
  cout << "│  0.0 ± 0.0   │" << endl;
  cout << "┴──────────────┴" << endl;

  cout << "┬──────────────┬" << endl;
  cout << "│ signal AR    │" << endl;
  cout << "┼──────────────┼" << endl;
  cout << "│  0.0 ± 0.0   │" << endl;
  cout << "│ " << fixed << setprecision(2) << h3eein ->GetBinContent(1) << " ± " << h3eein ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h3emin ->GetBinContent(1) << " ± " << h3emin ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h3mmin ->GetBinContent(1) << " ± " << h3mmin ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h3eeout->GetBinContent(1) << " ± " << h3eeout->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h3emout->GetBinContent(1) << " ± " << h3emout->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h3mmout->GetBinContent(1) << " ± " << h3mmout->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h3ee1J ->GetBinContent(1) << " ± " << h3ee1J ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h3em1J ->GetBinContent(1) << " ± " << h3em1J ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h3mm1J ->GetBinContent(1) << " ± " << h3mm1J ->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h30SFOS->GetBinContent(1) << " ± " << h30SFOS->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h31SFOS->GetBinContent(1) << " ± " << h31SFOS->GetBinError(1) << "  │" << endl;
  cout << "│ " << fixed << setprecision(2) << h32SFOS->GetBinContent(1) << " ± " << h32SFOS->GetBinError(1) << "  │" << endl;
  cout << "│  0.0 ± 0.0   │" << endl;
  cout << "┴──────────────┴" << endl;



}
