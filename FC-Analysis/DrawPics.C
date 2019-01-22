using namespace std;
TGraph* HistToGraph(TH1I* pH)
{
    int n = pH->GetNbinsX();
    Double_t x[n], y[n];
    for (int i = 1; i <= n; i++)
    {
        x[i] = pH->GetBinCenter(i);
        y[i] = pH->GetBinContent(i);
    }
    TGraph* pG = new TGraph(n, x, y);
    return pG;
}

void DrawQDCfit(TFile* f)
{
    /// Draw QDC fits
    TCanvas* pCRawQDC[8];
    TCanvas* pCAnaQDC[8];
    char hname[64] = "";
    for (int i = 0; i < 8; i++)
    {
        /// Draw fits for raw QDC spectrum
        // open histogram
        sprintf(hname, "/Histograms/Raw/QDC/low/H1RawQDCl_%i", i + 1);
        TH1I *pHRawQDC = (TH1I*)f->Get(hname);
        // open fits
        sprintf(hname, "/Analysis/QDC/Raw/Fit/fRawPed_%i", i+1);
        TF1* fRawPed = (TF1*)f->Get(hname);
        sprintf(hname, "/Analysis/QDC/Raw/Fit/fRawCut_%i", i+1);
        TF1* fRawCut = (TF1*)f->Get(hname);
        sprintf(hname, "/Analysis/QDC/Raw/Fit/fRawMax_%i", i+1);
        TF1* fRawMax = (TF1*)f->Get(hname);
        // draw together
        pCRawQDC[i] = new TCanvas();
        pHRawQDC->Draw();
        fRawPed->Draw("same");
        fRawCut->Draw("same");
        fRawMax->Draw("same");

        /// Draw fits for self-triggered ToF-gated QDC spectrum
        // open histogram
        sprintf(hname, "/Histograms/Analysis/PuFC/QDC/low/trig/H1AnaQDCl_trig_%i", i + 1);
        TH1I *pHAnaQDC = (TH1I*)f->Get(hname);
        // open fits
        sprintf(hname, "/Analysis/QDC/Trig/Fit/fTrigCut_%i", i+1);
        TF1* fAnaCut = (TF1*)f->Get(hname);
        sprintf(hname, "/Analysis/QDC/Trig/Fit/fTrigMax_%i", i+1);
        TF1* fAnaMax = (TF1*)f->Get(hname);
        // draw together
        pCAnaQDC[i] = new TCanvas();
        pHAnaQDC->Draw();
        fAnaCut->Draw("same");
        fAnaMax->Draw("same");
    }
}

void DrawQDCres(TFile* f)
{
    /// QDC Fit results
    TCanvas* c1 = new TCanvas("c1", "Multigraph fit results", 200, 10, 700, 500);
    TMultiGraph* mg1 = new TMultiGraph();
    mg1->SetTitle("QDC Fit results; Deposit; QDC channel (low gain)");
    TLegend* legend = new TLegend(0.6, 0.2, 0.85, 0.40, "Legend");

    TGraph* g1 = (TGraph*)f->Get("/Analysis/QDC/Trig/gTrigCut");
    g1->SetMarkerColor(kBlue);
    mg1->Add(g1);
    legend->AddEntry(g1, "self-triggered", "lp");

    TGraph* g0 = (TGraph*)f->Get("/Analysis/QDC/Raw/gRawCut");
    g0->SetMarkerColor(kRed);
    mg1->Add(g0);
    legend->AddEntry(g0, "Raw", "lp");

    TGraph* g3 = (TGraph*)f->Get("/Analysis/QDC/Trig/gTrigMax");
    g3->SetMarkerColor(kBlue);
    mg1->Add(g3, "P");

    TGraph* g2 = (TGraph*)f->Get("/Analysis/QDC/Raw/gRawMax");
    g2->SetMarkerColor(kRed);
    mg1->Add(g2, "P");

    mg1->Draw("AP");
    mg1->GetYaxis()->SetRangeUser(0, 2000);
    legend->Draw();
    c1->Modified();
    c1->Update();

    TCanvas* c2 = new TCanvas("c2", "relative Cut positions", 200, 10, 700, 500);
    TMultiGraph* mg2 = new TMultiGraph();
    mg2->SetTitle("Relative minimum position; Deposit; Ratio");
    TLegend* l2 = new TLegend(0.6, 0.6, 0.85, 0.85, "Legend");

    TGraph* g5 = (TGraph*)f->Get("/Analysis/QDC/Trig/gTrigRelCut");
    g5->SetMarkerColor(kBlue);
    mg2->Add(g5);
    l2->AddEntry(g5, "self-triggered", "lp");

    TGraph* g4 = (TGraph*)f->Get("/Analysis/QDC/Raw/gRawRelCut");
    g4->SetMarkerColor(kRed);
    mg2->Add(g4);
    l2->AddEntry(g4, "Raw", "lp");

    mg2->Draw("AP");
    mg2->GetYaxis()->SetRangeUser(0, 2);
    l2->Draw();
    c2->Modified();
    c2->Update();
}

void DrawDtInt(TFile* f)
{
    /// Draw counting results
    //TODO: change this for common pics
    TCanvas* c3 = new TCanvas("c3", "TimeDiff Counting results", 200, 10, 700, 500);
    TMultiGraph* mg3 = new TMultiGraph();
    mg3->SetTitle("SF detection rate; Deposit; events per second");

    TGraph* g6 = (TGraph*)f->Get("/Analysis/TimeDiff/SF_Rate");
    g6->SetMarkerColor(kBlue);
    mg3->Add(g6);

    mg3->Draw("AP");
    c3->Modified();
    c3->Update();

    TCanvas* c4 = new TCanvas("c4", "TimeDiff Counting results", 200, 10, 700, 500);
    TMultiGraph* mg4 = new TMultiGraph();
    mg4->SetTitle("NIF detection rate; Deposit; events per second");

    TGraph* g7 = (TGraph*)f->Get("/Analysis/TimeDiff/NIF_Rate");
    g7->SetMarkerColor(kBlue);
    mg4->Add(g7);

    mg4->Draw("AP");
    c4->Modified();
    c4->Update();

    TCanvas* c5 = new TCanvas("c5", "TimeDiff Counting results", 200, 10, 700, 500);
    TMultiGraph* mg5 = new TMultiGraph();
    mg5->SetTitle("NIF to SF detection ratio; Deposit; Ratio");

    TGraph* g8 = (TGraph*)f->Get("/Analysis/TimeDiff/NIF_SF_Ratio");
    g8->SetMarkerColor(kBlue);
    mg5->Add(g8);

    mg5->Draw("AP");
    c5->Modified();
    c5->Update();
}

void DrawDtUg(TFile* f)
{ // Draw TimeDiff integration
    TCanvas* pCAnaDt[8];
    char hname[64] = "";
    for (int i = 0; i < 8; i++)
    {
        // open Dt histogram nr i
        sprintf(hname, "/Histograms/Analysis/PuFC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i+1);
        TH1I* pH = (TH1I*)f->Get(hname);
        // open underground fit
        sprintf(hname, "/Analysis/TimeDiff/H1AnaDtFit_%i", i+1);
        TGraph* fUg = (TGraph*)f->Get(hname);
        // draw histogram and graph together
        pCAnaDt[i] = new TCanvas();//"c1", "TimeDiff Underground", 200, 10, 700, 500);
        pH->SetAxisRange(62000, 80000, "X");
        pH->Draw();
        fUg->SetLineColor(kRed);
        fUg->SetLineWidth(2);
        fUg->Draw("same");
        // draw number
        Double_t ChPerBin = pH->GetBinWidth(0);
        Double_t BinOffset = pH->GetBinCenter(0);
        Int_t low = (fUg->GetX()[0] - BinOffset) / ChPerBin;
        Int_t up = (fUg->GetX()[1] - BinOffset) / ChPerBin;
        Double_t NIF = pH->Integral(low, up) - (up - low + 1) * fUg->GetY()[0];
        char message[32] = "";
        sprintf(message, "NIF events: %i", (int)(NIF+0.5));
        TText *t = new TText();
        t->SetNDC();
        t->DrawText(0.4, 0.8, message);
    }
}

int TreatFile(TFile* f)
{
//    DrawQDCfit(f);
//    DrawQDCres(f);
//    DrawDtInt(f);
    DrawDtUg(f);



    return 1;
}

int DrawPics()
{
    TFile* fNIF = TFile::Open("/home/hoffma93/Go4nfis/offline/results/NIF.root");
    TFile* fSB = TFile::Open("/home/hoffma93/Go4nfis/offline/results/SB.root");
    /// single pics
    TreatFile(fNIF);
    TreatFile(fSB);

    /// common pics
    // SF, NIF
    return 1;
}
