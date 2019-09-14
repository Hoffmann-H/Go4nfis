TH1F* GetToF(TFile *f, Int_t ch, string FC = "PuFC")
{
    char name[64] = "";
    sprintf(name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", ch+1);
    TH1I* pH1I= (TH1I*) f->Get(name);
    Int_t nbins = pH1I->GetNbinsX();
    Int_t xmin = pH1I->GetBinLowEdge(1);
    Int_t xmax = pH1I->GetBinLowEdge(nbins + 1);
    TH1F* pH1F = new TH1F(pH1I->GetName(), "ToF; t / ns; counts", nbins, xmin, xmax);
    pH1F->Add(pH1I, 1);
    pH1F->Sumw2();
    return pH1F;
}


void ToF()
{
    TFile *f = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/NIF.root");
    TH1F *pH = GetToF(f, 0);
    TCanvas *c1 = new TCanvas();
    pH->Draw("hist");
    c1->Draw();
}
