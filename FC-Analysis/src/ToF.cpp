// Analysen der ToF-Spektren
#include "ToF.h"
#include "TFile.h"

using namespace std;

ToF::ToF(string file_name, string fc, string setup, string name)
{
    FilePath = file_name;
    FC = fc;
    Setup = setup;
//    std::stringstream s;
//    s << FC << "_" << Setup;
//    Name = s.str();
    Name = name;
    PuFC = strcmp(FC.c_str(), "UFC");
    Save = kTRUE;
    if (!strcmp(setup.c_str(), "SF"))
    {
        SpontaneousFission();
        return;
    }
    OpenToF(file_name);
//    MakeLimits(12, 37); //12.0821, 36.7100);//

//    Print();
//    FitBackground();
//    SubtractBackground();
//    InducedFission();
//    DrawDt(0);
}


ToF::~ToF()
{

}


string ToF::GetFCname()
{
    return FC;
}


Bool_t ToF::IsPuFC()
{
    return PuFC;
}


TH1F* ToF::TH1ItoTH1F(TH1I* pH)
{
    // copy manually...
    Int_t nbins = pH->GetNbinsX();
    Int_t xmin = pH->GetBinLowEdge(1);
    Int_t xmax = pH->GetBinLowEdge(nbins + 1);
//    cout << nbins << "   " << xmin << "   " << xmax << endl;
//    sprintf(name, "%s", pH->GetName());
    TH1F* pH1F = new TH1F(pH->GetName(), "ToF; t / ns; counts", nbins, xmin, xmax);
//    for (Int_t bin = 0; bin < nbins + 2; bin++)
//        pH1Dt[i]->SetBinContent(bin, pH->GetBinContent(bin));
    pH1F->Add(pH, 1);
    return pH1F;
}


void ToF::OpenToF(string file_name)
{
    cout << endl << "Opening histograms..." << endl;
    TFile *f = TFile::Open(file_name.c_str());

    TH1D* pHt = (TH1D*) f->Get("/Histograms/Raw/Scaler/Rates/H1RawRate_47")->Clone();
    t_live = pHt->Integral();

    char name[64] = "";
    char title[128] = "";
//    Double_t ch_ns = 40.97; // TDC channels per nanosecond
//    Double_t mPu[] = {5720, 5240, 5150, 5300, 5300, 5260, 5210, 5210}; // TDC channels from start of spectrum
//    Double_t mU[] = {11190, 11010, 10660, 10850, 10850, 10820, 10780, 10770};
    for (Int_t i = 0; i < NumCh; i++)
    {
        //// Open TDC channel spectrum
        sprintf(name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i + 1);
        TH1I* pH = (TH1I*) f->Get(name);
        pH1Dt[i] = TH1ItoTH1F(pH);
//        pH1Dt[i] = (TH1F*) f->Get(name)->Clone();
        pH1Dt[i]->Sumw2();
//        cout << "Opened " << pH1Dt[i]->GetName() << endl;
//        cout << " Integral: " << pH1DtCh[i]->Integral() << endl;

    }
//    f->Close();
    cout << "Done: Open histograms" << endl;
}


Float_t ToF_BG_min[]   = {25., 850.};               //ns
Float_t ToF_BG_max[]   = {95.,2450.};               //ns
Bool_t reject;
Double_t fline(Double_t *x, Double_t *par)
{   //compute the function only in selected sub-ranges
    if (reject)
    {   //only include points in sub-range into the fit
        if ((x[0] >= ToF_BG_min[0] && x[0] <= ToF_BG_max[0]) ||
            (x[0] >= ToF_BG_min[1] && x[0] <= ToF_BG_max[1]))
            return par[0];
        else
        {   TF1::RejectPoint();
            return 0;
        }
    }
    //compute the function values over the whole range
    else
        return par[0];
}


void ToF::MakeLimits(Double_t left, Double_t right)
{ // left, right: widths in ns. left = right = 0: Spontaneaus fission.
    cout << endl << "Making integration limits..." << endl;
    l0 = 42;
    l3 = 402;
    ToF_BG_min[0] = pH1Dt[0]->GetXaxis()->GetBinLowEdge(l0);
    ToF_BG_max[1] = pH1Dt[0]->GetXaxis()->GetBinLowEdge(l3);
    if (left == 0 && right == 0)
        return; // Spontaneous fission

    Double_t mPu[] = {140, 128, 126, 129, 129, 128, 127, 127}; // ns from start of spectrum
    Double_t mU[] = {273, 269, 260, 265, 265, 264, 263, 263};
    Double_t ns_bin = pH1Dt[0]->GetBinWidth(1);
    cout << " ch   bin-nr" << endl;
    for (Int_t i = 0; i < NumCh; i++)
    {
        if (this->IsPuFC())
            m[i] = mPu[i] / ns_bin;
        else
            m[i] = mU[i] / ns_bin;
        l1[i] = m[i] - left / ns_bin;
        l2[i] = m[i] + right / ns_bin;
        cout << " " << i+1 << "   " << l0 << "   " << l1[i] << "   " << l2[i] << "   " << l3 << endl;
    }
    cout << "Done: limits" << endl;
    FitBackground();
    SubtractBackground();
    InducedFission();
}


void ToF::SetLimits(Int_t *pl0, Int_t *pl1, Int_t *pl2, Int_t *pl3)
{
    cout << endl << Name << " setting integration limits..." << endl;
//    cout << " ch   bin-nr" << endl;
    l0 = *pl0;
    l3 = *pl3;
    ToF_BG_min[0] = pH1Dt[0]->GetXaxis()->GetBinLowEdge(l0);
    ToF_BG_max[1] = pH1Dt[0]->GetXaxis()->GetBinLowEdge(l3);
    for (Int_t i = 0; i < NumCh; i++)
    {
        l1[i] = pl1[i];
        l2[i] = pl2[i];
//        cout << " " << i+1 << "   " << l0 << "   " << l1[i] << "   " << l2[i] << "   " << l3 << endl;
    }
    cout << "Done: limits" << endl;
    FitBackground();
    SubtractBackground();
    InducedFission();
}


void ToF::SpontaneousFission()
{
    cout << endl << "Spontaneaus fission..." << endl;
    TFile *f = TFile::Open(FilePath.c_str());
    char name[64] = "";
//    char title[128] = "";
    for (Int_t i = 0; i < NumCh; i++)
    {
        //// Open TDC channel spectrum
        sprintf(name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i + 1);
        pH1Dt[i] = (TH1F*) f->Get(name)->Clone();
//        pH1Dt[i]->Sumw2();
        cout << "Opened " << pH1Dt[i]->GetName() << endl;
    }
    MakeLimits(0, 0);

    TH1D* pHt = (TH1D*) f->Get("/Histograms/Raw/Scaler/Rates/H1RawRate_47")->Clone();
    t_live = pHt->Integral();

    TGraphErrors *sf = new TGraphErrors(NumCh);
    sf->SetName("effSponFis");
    sf->SetTitle("eff. Spontaneous Fission rate; #Delta#font[12]{t}; Rate [1/s]");

    reject = kFALSE;
    cout << " ch   Fit   IntAndErr   user" << endl;
    for (Int_t i = 0; i < NumCh; i++)
    {
        sf->SetPoint(i, i+1, pH1Dt[i]->Integral() / t_live);
        sf->SetPointError(i, 0, sqrt(pH1Dt[i]->Integral()) / t_live);

        Double_t D1bg;
        Double_t bg = pH1Dt[i]->IntegralAndError(l0, l3-1, D1bg) / (l3 - l0);
        D1bg /= (l3 - l0);
        Double_t D0bg = sqrt(pH1Dt[i]->Integral(l0, l3-1)) / (l3 - l0);

        sprintf(name, "%s_fT_%i", Name.c_str(), i+1);
//        cout << name << "  " << ToF_BG_min[0] << "-" << ToF_BG_max[1] << endl;
        fTotal[i] = new TF1(name, fline, ToF_BG_min[0], ToF_BG_max[1], 1);
        pH1Dt[i]->Fit(name, "LR0Q");
        SaveToFile("Analysis/ToF/BG/Total", fTotal[i]);
        cout << " " << i+1 << "   " << fTotal[i]->GetParameter(0) << "+-" << fTotal[i]->GetParError(0)
                           << "   " << bg << "+-" << D1bg
                           << "   " << pH1Dt[i]->Integral(l0, l3-1) / (l3 - l0) << "+-" << D0bg << " " << l3 - l0 << endl;
//        DrawDt(i);
    }
    SaveToFile("Analysis/results", sf);
}


void ToF::DrawDt(Int_t i)
{
    cout << "Draw " << pH1Dt[i]->GetName() << endl;
    char name[64] = "";
    sprintf(name, "c_%s_Dt_%i", Name.c_str(), i + 1);
    TCanvas* c = new TCanvas(name, name, 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    pH1Dt[i]->SetStats(0);
    pH1Dt[i]->Draw("HIST");
    reject = kFALSE;
    fTotal[i]->Draw("same");
    Double_t ymin = 0;
    Double_t ymax = pH1Dt[i]->GetMaximum();
    TLine* line0 = new TLine(ToF_BG_min[0], ymin, ToF_BG_min[0], ymax);
    TLine* line3 = new TLine(ToF_BG_max[1], ymin, ToF_BG_max[1], ymax);
    line0->Draw();
    line3->Draw();
    c->Modified();
    c->Update();
}



void ToF::DrawDt()
{
    cout << endl << "Drawing ToF spectra..." << endl;
    TCanvas* c = new TCanvas("name", "title", 200, 10, 700, 500);
    c->Divide(2, 4);
    for (Int_t i = 0; i < NumCh; i++)
    {
        c->cd(i + 1);
        pH1Dt[i]->Draw("HIST");
        cout << " ch " << i+1 << endl;
    }
    c->Modified();
    c->Update();
    cout << "Done: Draw ToF spectra" << endl;
}


void ToF::FitBackground()
{
    cout << endl << "Fit constant ToF background..." << endl;
    cout << " ch   left   right   total   chi^2/dof   user" << endl;
    for (Int_t i = 0; i < NumCh; i++)
        FitBackground(i);
    cout << "Done: Fit constant ToF background" << endl;
}


void ToF::FitBackground(Int_t  i)
{
    char name[64] = "";
    Double_t xmin = pH1Dt[i]->GetBinLowEdge(1);
    Double_t xmax = pH1Dt[i]->GetBinLowEdge(pH1Dt[i]->GetNbinsX() + 2);
//    cout << "DEBUG: " << pH1Dt[i]->GetName() << ", xmin " << xmin << ", xmax " << xmax << ", " << pH1Dt[i]->GetNbinsX() << endl;
//    ToF_BG_min[0] = pH1Dt[i]->GetXaxis()->GetBinLowEdge(l0);
    ToF_BG_max[0] = pH1Dt[i]->GetXaxis()->GetBinLowEdge(l1[i]);
    ToF_BG_min[1] = pH1Dt[i]->GetXaxis()->GetBinLowEdge(l2[i]);
//    ToF_BG_max[1] = pH1Dt[i]->GetXaxis()->GetBinLowEdge(l3);
    Int_t npar = 1;

    sprintf(name, "%s_fL_%i", Name.c_str(), i+1);
    fLeft[i] = new TF1(name, fline, ToF_BG_min[0], ToF_BG_max[0], npar);
    reject = kTRUE;
    pH1Dt[i]->Fit(name, "LR0Q");
    reject = kFALSE;
    SaveToFile("Analysis/ToF/BG/Left", fLeft[i]);

    sprintf(name, "%s_fR_%i", Name.c_str(), i+1);
    fRight[i] = new TF1(name, fline, ToF_BG_min[1], ToF_BG_max[1], npar);
    reject = kTRUE;
    pH1Dt[i]->Fit(name, "LR0Q");
    reject = kFALSE;
    SaveToFile("Analysis/ToF/BG/Right", fRight[i]);

    sprintf(name, "%s_fT_%i", Name.c_str(), i+1);
    fTotal[i] = new TF1(name, fline, xmin, xmax, npar);
    reject = kTRUE;
    pH1Dt[i]->Fit(name, "LR0Q");
    reject = kFALSE;
    SaveToFile("Analysis/ToF/BG/Total", fTotal[i]);

    ug[i] = (pH1Dt[i]->Integral(l0, l1[i]-1) + pH1Dt[i]->Integral(l2[i], l3-1)) / (l1[i] - l0 + l3 - l2[i]);
    Dug[i] = sqrt(pH1Dt[i]->Integral(l0, l1[i]-1) + pH1Dt[i]->Integral(l2[i], l3-1)) / (l1[i] - l0 + l3 - l2[i]);

    cout << " " << i+1 << "   " << fLeft[i]->GetParameter(0) << "+-" << fLeft[i]->GetParError(0)
                       << "   " << fRight[i]->GetParameter(0) << "+-" << fRight[i]->GetParError(0)
                       << "   " << fTotal[i]->GetParameter(0) << "+-" << fTotal[i]->GetParError(0)
                       << "   " << fTotal[i]->GetChisquare() << "/" << fTotal[i]->GetNDF() << "=" << fTotal[i]->GetChisquare() / fTotal[i]->GetNDF()
                       << "   " << ug[i] << "+-" << Dug[i] << endl;
}


void ToF::SubtractBackground()
{
    cout << endl << "Subtract constant background" << endl;
//    cout << pH1Dt[0]->GetBinContent(100) << endl;
//    cout << " ch   bg-bg=0?" << endl;
    for (Int_t i = 0; i < NumCh; i++)
        SubtractBackground(i);
//    cout << pH1DtSub[0]->GetBinContent(100) << endl;
    cout << "Done: Subtract constant background" << endl;
}


void ToF::SubtractBackground(Int_t i)
{
    char name[64] = "";
    reject = kFALSE;
    pH1DtSub[i] = (TH1F*)pH1Dt[i]->Clone();
    pH1DtSub[i]->Add(fTotal[i], -1);
    SaveToFile("Analysis/ToF", pH1DtSub[i]);

    // Test: bg=0?
    Double_t ug = (pH1DtSub[i]->Integral(l0, l1[i]-1) + pH1DtSub[i]->Integral(l2[i], l3-1)) / (l1[i] - l0 + l3 - l2[i]);
    Double_t Dug = sqrt(pH1DtSub[i]->Integral(l0, l1[i]-1) + pH1DtSub[i]->Integral(l2[i], l3-1)) / (l1[i] - l0 + l3 - l2[i]);
//    cout << " " << i+1 << "   " << ug << "+-" << Dug << endl;

    // Create TGraphErrors
    Int_t nBins = pH1DtSub[i]->GetNbinsX();
    Double_t X[nBins];
    Double_t Xerr[nBins];
    Double_t Y[nBins];
    Double_t Yerr[nBins];
    pH1DtSub[i]->GetXaxis()->GetCenter(X);
    for (Int_t j = 0; j < nBins; j++)
    {
        Xerr[j] = 0;
        Y[j] = pH1DtSub[i]->GetBinContent(j+1);
        Yerr[j] = pH1DtSub[i]->GetBinError(j+1);
    }
    TGraphErrors *ge = new TGraphErrors(nBins, X, Y, Xerr, Yerr);

    // Fit zero
    Double_t xmin = pH1DtSub[i]->GetBinLowEdge(1);
    Double_t xmax = pH1DtSub[i]->GetBinLowEdge(nBins + 2);
    sprintf(name, "%s_fZ_%i", Name.c_str(), i+1);
    fZero[i] = new TF1(name, fline, xmin, xmax, 1);
    reject = kTRUE;
    ge->Fit(name, "R0Q");
    fZero[i]->SetParameter(0, fZero[i]->GetParameter(0) + 1.0);
    reject = kFALSE;
    SaveToFile("Analysis/ToF/BG/Zero", fZero[i]);
//    cout << " " << i+1 << "   " << fZero[i]->GetParameter(0) << "+-" << fZero[i]->GetParError(0) << endl;
}


void ToF::SaveToFile(string path, TObject *pObj)
{   //saves a TObject into the selected file fname
    if (!Save)
        return;
//    cout << "Saving " << pObj->GetName() << " to " << path << endl;
    TFile* f = TFile::Open(FilePath.c_str(), "UPDATE");
    TDirectory *EvalDir;
    TObject *pGraph;
    pGraph = (TObject*) pObj;
    string GraphName = pGraph->GetName();
    //check if folder "Analysis" already exists, otherwise create it
    if (f->Get(path.c_str())!=0)
        EvalDir = (TDirectory*) f->Get(path.c_str());
    else
    {   cout << "Creating " << path << endl;
        f->mkdir(path.c_str(), "Folder containing offline Analysis objects");
        EvalDir = f->GetDirectory(path.c_str());
    }
    EvalDir->cd();
    if (EvalDir->Get(GraphName.c_str())!=0)
        EvalDir->Delete((GraphName+";*").c_str());
    pGraph->Clone()->Write();
    f->Save(); f->Close();
}


void ToF::InducedFission()
{
    cout << endl << "Induced fission events..." << endl;

    TGraphErrors *gSF = new TGraphErrors(NumCh);
    gSF->SetName("effSponFis");
    gSF->SetTitle("eff. Spontaneous Fission rate; #Delta#font[12]{t}; Rate [1/s]");

    TGraphErrors *gNF = new TGraphErrors(NumCh);
    gNF->SetName("effIndFis");
    gNF->SetTitle("eff. n-Induced Fission rate; #Delta#font[12]{t}; Rate [1/s]");

    cout << " ch   (n,f)   user   sf" << endl;
    for (Int_t i = 0; i < NumCh; i++)
    {
        InducedFission(i);
        gSF->SetPoint(i, i+1, sf[i] / t_live);
        gSF->SetPointError(i, 0, Dsf[i] / t_live);
        gNF->SetPoint(i, i+1, nf[i] / t_live);
        gNF->SetPointError(i, 0, Dnf[i] / t_live);
    }
    SaveToFile("Analysis/results", gSF);
    SaveToFile("Analysis/results", gNF);
    cout << "Done: Induced fission" << endl;
}


void ToF::InducedFission(Int_t i)
{
    // first: Sum bg-subtracted spectrum
    Double_t DnNIF;
    Double_t nNIF = pH1DtSub[i]->IntegralAndError(l1[i], l2[i], DnNIF);
    nf[i] = nNIF;
    Dnf[i] = DnNIF;
    // second: Subtract bg from integral
    nNIF = pH1Dt[i]->Integral(l1[i], l2[i]) - (l2[i] - l1[i] + 1) * ug[i];
    DnNIF = sqrt( pH1Dt[i]->Integral(l1[i], l2[i]) + pow((l2[i] - l1[i] + 1) * Dug[i], 2) );
    unf[i] = nNIF;
    Dunf[i] = DnNIF;
    // Find spontaneous fission rate
    sf[i] = pH1Dt[i]->Integral() - nf[i];
    Dsf[i] = pH1Dt[i]->Integral() / sqrt(sf[i]); // = Poisson counting statistics excluding peak interval, scaled up.

    cout << " " << i+1 << "   " << nf[i] << "+-" << Dnf[i]
                       << "   " << unf[i] << "+-" << Dunf[i]
                       << "   " << sf[i] << "+-" << Dsf[i] << endl;
}


void ToF::Print()
{
    cout << endl;
    for (Int_t i = 0; i < NumCh; i++)
    {
        cout << l1[i] - l0 + l3 - l2[i] << " " << pH1Dt[i]->Integral(l0, l1[i]-1) + pH1Dt[i]->Integral(l2[i], l3-1) << " " << l2[i] - l1[i] << " " << pH1Dt[i]->Integral(l1[i], l2[i]-1) << " " << t_live << endl;
    }
}


void ToF::DrawPeakLim(Int_t l, Int_t r_start, Int_t r_stop, Int_t r_step)
{
    char name[64] = "";
    Int_t n = (r_stop - r_start) / r_step + 1;
    SetSave(kFALSE);
    TGraphErrors* ge[NumCh];
    sprintf(name, "%s Deposit", FC.c_str());
    TLegend *legend = new TLegend(0.9, 0.5, 1.0, 1.0, name);
    for (Int_t i = 0; i < NumCh; i++)
    {
        ge[i] = new TGraphErrors(n);
        sprintf(name, "gPeak_%i", i+1);
        ge[i]->SetName(name);
    }
    TGraphErrors* ge_av = new TGraphErrors(n);
    ge_av->SetName("gPeak_av");

    TMultiGraph* mg = new TMultiGraph("mgPeak", "ToF Peak Integration; Gate [ns]; Peak");
    for (Int_t i = 0; i < NumCh; i++)
    {
        mg->Add(ge[i]);
        sprintf(name, "%i", i+1);
        legend->AddEntry(ge[i], name);
    }
    mg->Add(ge_av);
    legend->AddEntry(ge_av, "av");

    for (Int_t j = 0; j < n; j++)
    {
        Int_t r = r_start + j * r_step;
        cout << endl << "DrawPeakLim starting analysis with r = " << r << "ns" << endl;
        MakeLimits(l, r);
        Double_t av = 0;
        Double_t D2av = 0;
        for (Int_t i = 0; i < NumCh; i++)
        {
            ge[i]->SetPoint(j, l + r, nf[i]);
            ge[i]->SetPointError(j, 0, Dnf[i]);
            av += nf[i] / NumCh;
            D2av += pow(Dnf[i] / NumCh, 2);
        }
        ge_av->SetPoint(j, l + r, av);
        ge_av->SetPointError(j, 0, sqrt(D2av));
    }
    SetSave(kTRUE);
    for (Int_t i = 0; i < NumCh; i++)
    {
        SaveToFile("Analysis/ToF/Gate", ge[i]);
        ge[i]->SetLineColorAlpha(i+2, 0.5);
        ge[i]->SetLineWidth(1);
        ge[i]->SetMarkerColorAlpha(i+2, 0.5);
        ge[i]->SetMarkerStyle(20);
        ge[i]->SetMarkerSize(2);
    }
    SaveToFile("Analysis/ToF/Gate", ge_av);
    ge_av->SetLineColorAlpha(1, 0.5);
    ge_av->SetLineWidth(2);
    ge_av->SetMarkerColorAlpha(1, 0.5);
    ge_av->SetMarkerStyle(20);
    ge_av->SetMarkerSize(2);
    TCanvas* c1 = new TCanvas("cPeak", "Test Integration Limits", 200, 10, 700, 500);
    c1->SetTicks(1, 1);
    mg->Draw("AP");
    legend->Draw();
    c1->Modified();
    c1->Update();
}


//void ToF::DrawDtSub(Int_t i)
//{
//    cout << "Draw " << pH1Dt[i]->GetName() << endl;
//    cout << pH1Dt[i]->GetNbinsX() << " bins" << endl;
//    char name[64] = "";
//    sprintf(name, "c_%s_Dt_%i", Name.c_str(), i + 1);
//    TCanvas* c = new TCanvas(name, name, 200, 10, 700, 500);
//    gPad->SetTicks(1, 1);
//    pH1Dt[i]->SetStats(0);
//    pH1Dt[i]->Draw("HIST");

////    reject = kFALSE;
////    fTotal[i]->Draw("same");
////    fLeft[i]->Draw("same");
////    fRight[i]->Draw("same");

//    Double_t x0 = pH1Dt[i]->GetXaxis()->GetBinLowEdge(l0);
//    Double_t x1 = pH1Dt[i]->GetXaxis()->GetBinLowEdge(l1[i]);
//    Double_t x2 = pH1Dt[i]->GetXaxis()->GetBinLowEdge(l2[i]+1);
//    Double_t x3 = pH1Dt[i]->GetXaxis()->GetBinLowEdge(l3+1);
//    Double_t ymin = 0;
//    Double_t ymax = pH1Dt[i]->GetMaximum();
//    cout << "x limits in ns: " << x0 << " " << x1 << " " << x2 << " " << x3 << endl
//         << " y limits: " << ymin << " " << ymax << endl;
//    TLine* line0 = new TLine(x0, ymin, x0, ymax);
//    TLine* line1 = new TLine(x1, ymin, x1, ymax);
//    TLine* line2 = new TLine(x2, ymin, x2, ymax);
//    TLine* line3 = new TLine(x3, ymin, x3, ymax);
//    line0->Draw();
//    line1->Draw();
//    line2->Draw();
//    line3->Draw();

//    c->Modified();
//    c->Update();
//}


//void ToF::DrawDtSub()
//{
//    cout << endl << "Drawing ToF spectra..." << endl;
//    TCanvas* c = new TCanvas("name", "title", 200, 10, 700, 500);
//    c->Divide(2, 4);
//    for (Int_t i = 0; i < NumCh; i++)
//    {
//        c->cd(i + 1);
//        pH1Dt[i]->Draw();
//        cout << " ch " << i+1 << endl;
//    }
//    c->Modified();
//    c->Update();
//    cout << "Done: Draw ToF spectra" << endl;
//}
