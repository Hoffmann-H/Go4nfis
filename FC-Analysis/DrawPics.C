#define lw 2
//#define a 0.25 // alpha
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


TGraph* Rebin(TF1* func, Int_t nPoints, Int_t yFactor)
{
    Double_t x[nPoints];
    Double_t y[nPoints];
    Double_t xmin, xmax;
    func->GetRange(xmin, xmax);
    for (int i = 0; i < nPoints; i++)
    {
        x[i] = xmin + i * (xmax - xmin) / (nPoints - 1);
        y[i] = yFactor * func->Eval(x[i]);
    }
    TGraph* TheGraph = new TGraph(nPoints, x, y);
    return TheGraph;
}


Double_t func_peak(Double_t *x, Double_t *par)
{
    return par[0] * exp( - pow((x[0] - par[1]) / par[2], 2)) + par[3];
}


void BiasX(TH1I* pH, Double_t dx)
{
    Int_t offset = dx / pH->GetBinWidth(0);
    Int_t underflow = 0;
    for (int bin = 0; bin < offset+1; bin++)
        underflow += pH->GetBinContent(bin);
    pH->SetBinContent(0, underflow);
    for (int bin = offset + 1; bin < pH->GetNbinsX() + 1; bin++)
        pH->SetBinContent(bin - offset, pH->GetBinContent(bin));
    for (int bin = pH->GetNbinsX() - offset + 1; bin < pH->GetNbinsX() + 1; bin++)
        pH->SetBinContent(bin, 0);
}


TH1F* CopyRange(TH1I* pH, char* name, Double_t x0, Double_t x1, Double_t yoffset)
{
    Double_t offset = pH->GetBinLowEdge(0);
    Double_t ChPerBin = pH->GetBinWidth(0);
    int b0 = (x0 - offset) / ChPerBin;
    int b1 = (x1 - offset) / ChPerBin + 1;
    TH1F* pH2 = new TH1F(name, name, b1 - b0, pH->GetBinLowEdge(b0), pH->GetBinLowEdge(b1));
    for (int bin = 0; bin < pH2->GetNbinsX(); bin++)
        pH2->SetBinContent(bin + 1, pH->GetBinContent(bin + b0) + yoffset);
    return pH2;
}


void DrawQDCfit(TFile* f, string FC, string Setup)
{
    /// Draw QDC fits
    TCanvas* pCRawQDC[8];
    TCanvas* pCAnaQDC[8];
    char hname[64] = "";
    char name[32] = "";
    for (int i = 0; i < 8; i++)
    {
        /// Draw fits for raw QDC spectrum
        sprintf(name, "%s_%s_QDCfit_%i", FC.c_str(), Setup.c_str(), i+1);
        pCRawQDC[i] = new TCanvas(name, name, 200, 10, 700, 500);
        gPad->SetTicks(1, 1);
        gPad->SetLogy(1);
        TLegend* legend = new TLegend(0.6, 0.2, 0.85, 0.70, "Legend");
        // open histogram
        sprintf(hname, "/Histograms/Raw/QDC/low/H1RawQDCl_%i", i + 1);
        TH1I *pHRawQDC = (TH1I*)f->Get(hname);
        pHRawQDC->SetStats(0);
        pHRawQDC->GetXaxis()->SetTitleSize(0.05);
        pHRawQDC->GetYaxis()->SetTitleSize(0.05);
        pHRawQDC->SetLineColor(kBlue);
        legend->AddEntry(pHRawQDC, "Data");
        pHRawQDC->Draw();
        // open fits
        sprintf(hname, "/Analysis/QDC/Fit/f%s%sCut_%i", FC.c_str(), Setup.c_str(), i+1);
        TF1* fCut = (TF1*)f->Get(hname);
        fCut->SetLineWidth(lw);
        fCut->SetLineColor(kRed);
        legend->AddEntry(fCut, "pol4 Minimum fit");
        fCut->Draw("same"); // draw together

//        sprintf(hname, "/Analysis/QDC/Fit/f%s%sMax_%i", FC.c_str(), Setup.c_str(), i+1);
//        TF1* fMax = (TF1*)f->Get(hname);
//        fMax->SetLineWidth(lw);//CanvasPreferGL
//        fMax->SetLineColor(kGreen);
//        legend->AddEntry(fMax, "pol4 Maximum fit");
//        fMax->Draw("same");

//        if(strcmp(FC.c_str(), "PuFC") == 0)
//        { // if drawing PuFC results, draw Pedestal fit.
//            sprintf(hname, "/Analysis/QDC/Fit/f%s%sPed_%i", FC.c_str(), Setup.c_str(), i+1);
//            TF1* fPed = (TF1*)f->Get(hname);
//            fPed->SetLineColor(kYellow);
//            legend->AddEntry(fPed, "pol2 Pedestal fit");
//            fPed->Draw("same");
//        }
        legend->Draw("same");
    }
}

void DrawQDCeff(TFile* f, string FC, string Setup)
{ // Plot internal efficiency on self-triggered QDC spectrum
    TCanvas* pC[8];
    char hname[67] = "";
    char name[67] = "";
    TGraphErrors* eInt = (TGraphErrors*) f->Get("/Analysis/QDC/eInt");
    Double_t x, eff, Deff;
    for (int i = 0; i < 8; i++)
    {
        // open Dt histogram nr i
        sprintf(hname, "/Histograms/Analysis/FC/QDC/low/trig/H1AnaQDCl_trig_%i", i + 1);
        TH1I *pH = (TH1I*)f->Get(hname);
        pH->SetLineColor(kBlue);
        pH->SetStats(0);
        pH->GetXaxis()->SetLabelSize(0.06);
        pH->GetXaxis()->SetTitleSize(0.08);
        pH->GetYaxis()->SetLabelSize(0.06);
        pH->GetYaxis()->SetTitleSize(0.08);
        pH->GetYaxis()->SetTitleOffset(0.9);
        sprintf(name, "; #font[12]{Q} [Kanal]; #frac{Ereignisse}{Kanal}");
        pH->SetTitle(name);

        // open underground line
        sprintf(hname, "/Analysis/QDC/Ug/f%s%sUg_%i", FC.c_str(), Setup.c_str(), i+1);
        TGraph* fUg = (TGraph*)f->Get(hname);
        Double_t level = 0, cut = 0;
        fUg->GetPoint(3, cut, level);
        TH1I *pH1 = (TH1I*)pH->Clone();
        sprintf(hname, "QDC signal integration, channel %i", i+1);
        TH1F* pH2 = (TH1F*)CopyRange(pH, hname, cut, 4096, 0);
        pH2->SetLineWidth(0);

        // Make legend
//        TLegend *l = new TLegend(0.6, 0.85, 0.95, 0.95, "");
//        sprintf(name, "%s Kanal %i", FC.c_str(), i+1);
//        l->AddEntry(pH1, name);

        // draw histogram and graph together
        sprintf(name, "%s_%s_QDCeff_%i", FC.c_str(), Setup.c_str(), i+1);
        pC[i] = new TCanvas(name, name, 200, 10, 700, 500);
        gPad->SetTicks(1, 1);
        gPad->SetLogy(1);
//        pH1->SetAxisRange(0, 4096, "X");
        pH1->Draw();
        pH2->SetFillColorAlpha(kBlue, 0.25);
        pH2->Draw("same");
        fUg->SetLineColor(kRed);
        fUg->SetLineWidth(lw);
        fUg->SetFillColorAlpha(kRed, 0.25);
//        legend->AddEntry(fUg, "FF constant extrapolation");
        fUg->Draw("fsame");
//        l->Draw();

        // Get eff number
        eInt->GetPoint(i, x, eff);
        Deff = eInt->GetErrorY(i);
        sprintf(name, "#varepsilon #approx %.3f", eff);
        TLatex* tNIF = new TLatex();
        tNIF->SetNDC();
        tNIF->DrawLatex(0.25, 0.25, name);
//        TLatex* thr = new TLatex();
//        thr->SetNDC();
//        thr->DrawLatex(0.25, 0.5, "thr");
        TLatex* FF = new TLatex();
        FF->SetNDC();
        FF->DrawLatex(0.6, 0.65, "Spaltprodukte");
        TLatex* alpha = new TLatex();
        alpha->SetNDC();
        alpha->DrawLatex(0.25, 0.65, "#alpha");

        Double_t x[] = {cut, cut};
        Double_t y[] = {level, pH->GetMaximum()};
//        TArrow *arrow = new TArrow(0.5, 0.6, 0.9, 0.6);
//        arrow->Draw();

        pC[i]->Modified();
        pC[i]->Update();
    }
}

void DrawQDCres(TFile* f, string FC, string Setup)
{ // Plot QDC fit results
    char name[32] = "";
    sprintf(name, "%s_%s_QDCfit", FC.c_str(), Setup.c_str());
    /// QDC Fit results
    TCanvas* c1 = new TCanvas(name, name, 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg1 = new TMultiGraph();
    mg1->SetTitle("QDC Fit results; Deposit; QDC channel (low gain)");
    TLegend* legend = new TLegend(0.6, 0.2, 0.85, 0.70, "Legend");

    TGraph* g1 = (TGraph*)f->Get("/Analysis/QDC/gPed");
    g1->SetMarkerStyle(20);
    g1->SetMarkerColor(kBlue);
    mg1->Add(g1);
    legend->AddEntry(g1, "Pedestal position", "lp");

    TGraph* g0 = (TGraph*)f->Get("/Analysis/QDC/gCut");
    g0->SetMarkerColor(kRed);
    mg1->Add(g0);
    legend->AddEntry(g0, "Minimum position", "lp");

    TGraph* g2 = (TGraph*)f->Get("/Analysis/QDC/gMax");
    g2->SetMarkerColor(kGreen);
    mg1->Add(g2, "P");
    legend->AddEntry(g2, "Maximum position", "lp");

    mg1->Draw("AP");
    mg1->GetYaxis()->SetRangeUser(0, 2000);
    legend->Draw();
    c1->Modified();
    c1->Update();

//    sprintf(name, "%s_%s_QDCcut", FC.c_str(), Setup.c_str());
//    TCanvas* c2 = new TCanvas(name, name, 200, 10, 700, 500);
//    gPad->SetTicks(1, 1);
//    TMultiGraph* mg2 = new TMultiGraph();
//    mg2->SetTitle("Relative minimum position; Deposit; Ratio");

//    TGraph* g7 = (TGraph*)f->Get("/Analysis/QDC/gRelCut");
//    g7->SetMarkerColor(kRed);
//    mg2->Add(g7);

//    mg2->Draw("AP");
//    mg2->GetYaxis()->SetRangeUser(0, 2);
//    c2->Modified();
//    c2->Update();
}

void DrawDtInt(TFile* f, string FC, string Setup)
{ // Plot TimeDiff Integration results
    char name[64] = "";
    sprintf(name, "%s_%s_rSF", FC.c_str(), Setup.c_str());
    /// Draw counting results
    TCanvas* c3 = new TCanvas(name, name, 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg3 = new TMultiGraph();
    mg3->SetTitle("SF detection rate; Deposit; events per second");

    sprintf(name, "/%s/TimeDiff/Underground/SF_Rate", FC.c_str());
    TGraph* g6 = (TGraph*)f->Get(name);
    g6->SetLineWidth(lw);
    g6->SetMarkerColor(kBlue);
    mg3->Add(g6);

    mg3->Draw("AP");
    c3->Modified();
    c3->Update();

    if (strcmp(Setup.c_str(), "SF") != 0)
    {
        sprintf(name, "%s_%s_rNIF", FC.c_str(), Setup.c_str());
        TCanvas* c7 = new TCanvas(name, name, 200, 10, 700, 500);
        TMultiGraph* mg7 = new TMultiGraph();
        mg7->SetTitle("NIF detection rate; Deposit; events per second");

        sprintf(name, "/%s/TimeDiff/NIF_Rate", FC.c_str());
        TGraph* g7 = (TGraph*)f->Get(name);
        g7->SetLineWidth(lw);
        g7->SetMarkerColor(kBlue);
        mg7->Add(g7);

        mg7->Draw("AP");
        c7->Modified();
        c7->Update();

        sprintf(name, "%s_%s_NIFtoSF", FC.c_str(), Setup.c_str());
        TCanvas* c5 = new TCanvas(name, name, 200, 10, 700, 500);
        gPad->SetTicks(1, 1);
        TMultiGraph* mg5 = new TMultiGraph();
        mg5->SetTitle("NIF to SF detection ratio; Deposit; Ratio");

        sprintf(name, "/%s/TimeDiff/NIF_SF_Ratio", FC.c_str());
        TGraph* g8 = (TGraph*)f->Get(name);
        g8->SetLineWidth(lw);
        g8->SetMarkerColor(kBlue);
        mg5->Add(g8);

        mg5->Draw("AP");
        c5->Modified();
        c5->Update();

        sprintf(name, "%s_%s_NIFtoN", FC.c_str(), Setup.c_str());
        TCanvas* c6 = new TCanvas(name, name, 200, 10, 700, 500);
        gPad->SetTicks(1, 1);
        TMultiGraph* mg6 = new TMultiGraph();
        mg6->SetTitle("NIF rate normed to neutron flux; Deposit; Ratio / mm^2 s");

        sprintf(name, "/%s/TimeDiff/NIF_Flux_Ratio", FC.c_str());
        TGraph* g9 = (TGraph*)f->Get(name);
        g9->SetLineWidth(lw);
        g9->SetMarkerColor(kBlue);
        mg6->Add(g9);

        mg6->Draw("AP");
        c6->Modified();
        c6->Update();
    }
}


void DrawUGDt(TFile* f, TFile* fCom, string FC, string Setup)
{// Draw Dt spectrum for beam off.
    TCanvas* pCDtUg[8];
    char name[64] = "";
    for (int i = 0; i < 8; i++)
    {
        // open Dt histogram nr i
        sprintf(name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i+1);
        TH1I* pH = (TH1I*)f->Get(name);

        // open underground fit nr i
        sprintf(name, "/%s/TimeDiff/Underground/UG/f%sUGUg_%i", FC.c_str(), FC.c_str(), i+1);
        TGraph* fUg = (TGraph*)fCom->Get(name);
        fUg->SetLineColor(kRed);
        fUg->SetLineWidth(lw);

        // Draw together
        sprintf(name, "%s_%s_DtUg_%i", FC.c_str(), Setup.c_str(), i+1);
        pCDtUg[i] = new TCanvas(name, name, 200, 10, 700, 500);
        gPad->SetTicks(1,1);
        pH->GetXaxis()->SetRangeUser(62000, 80000);
        pH->Draw();
        fUg->Draw("same");
    }
}


TH1I* GetDt(TFile* f, string FC, string Setup, Int_t channel)
{
    char name[64] = "";
    char title[64] = "";
    Double_t chMin, chMax, tMin, tMax;
    Double_t TimePerCh = 25.E-03; // 25 ps in ns

    // Open PH-gated TimeDiff spectrum
    sprintf(name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", channel+1);
    TH1I *pHDtG = (TH1I*) f->Get(name);
    Int_t nbins = pHDtG->GetXaxis()->GetNbins();
    chMin = pHDtG->GetXaxis()->GetBinLowEdge(1);
    chMax = pHDtG->GetXaxis()->GetBinUpEdge(nbins);
//        ChPerBin = pHDtG->GetXaxis()->GetBinWidth(0);

    // Create counts over time histogram
    sprintf(name, "%s_%s_DtG_%i", FC.c_str(), Setup.c_str(), channel+1);
    sprintf(title, "; #font[12]{t} [ns]; Counts");
    tMin = TimePerCh * chMin;
    tMax = TimePerCh * chMax;
    TH1I *pHDtRate = new TH1I(name, title, nbins, tMin, tMax);
    pHDtRate->GetXaxis()->SetRangeUser(1550, 2000);
    pHDtRate->GetXaxis()->SetTitleSize(0.05);
    pHDtRate->GetYaxis()->SetTitleSize(0.05);
    pHDtRate->GetYaxis()->SetTitleOffset(1.0);
    pHDtRate->SetStats(0);
    pHDtRate->SetLineWidth(lw);

    // Fill histogram
    Int_t counts;
    for (int j = 0; j < nbins + 2; j++)
    {
        counts = pHDtG->GetBinContent(j);
        pHDtRate->SetBinContent(j, counts);
    }
    return pHDtRate;
}


void DrawDtUg(TFile* f, TFile* fCom, string FC, string Setup)
{ // Draw TimeDiff integration
    if(!strcmp(Setup.c_str(), "SF"))
    { // if beam off, call the corresponding function.
        DrawUGDt(f, fCom, FC, Setup);
        return;
    }
    char name[64] = "";
    TCanvas* pCDtPeak[8];
    TCanvas* pCDtUg[8];
    TH1I* pHsum = new TH1I();
    Double_t UgX[2];
    Double_t UgY[2];
    char hname[64] = "";
    char message[64] = "";
    for (int i = 0; i < 8; i++)
    {
        // open Dt histogram nr i
        TH1I* pH = (TH1I*)GetDt(f, FC, Setup, i);

        // open underground fit nr i
        sprintf(hname, "/%s/TimeDiff/Underground/%s/f%s%sUg_%i", FC.c_str(), Setup.c_str(), FC.c_str(), Setup.c_str(), i+1);
        TGraph* fUg = (TGraph*)fCom->Get(hname);
        sprintf(hname, "/%s/TimeDiff/Underground/%s/f%s%sUgPeak_%i", FC.c_str(), Setup.c_str(), FC.c_str(), Setup.c_str(), i+1);
        TGraph* fUgP = (TGraph*)fCom->Get(hname);

        //// Draw peak integration
        // extract peak integration limits from underground graph
        Double_t x0 = 0, x1 = 0, level = 0; // prepare ug level loading
        fUgP->GetPoint(0, x0, level);
        fUgP->GetPoint(1, x1, level);

        // create full TH1F histogram
        sprintf(name, "%s_%s_Dt_%i", FC.c_str(), Setup.c_str(), i+1);
        TH1F* pH1 = CopyRange(pH, name, 1550, 2000, - level);
        pH1->SetTitle("Untergrund-bereinigtes Zeitdifferenz-Spektrum; #font[12]{t} [ns]; Counts");

        // create peak region TH1F histogram
        sprintf(name, "%s_%s_DtInt_%i", FC.c_str(), Setup.c_str(), i+1);
        TH1F* pH2 = CopyRange(pH, name, .025*x0, .025*x1, - level);

        // manipulate histograms
        TH1I* pH3 = (TH1I*)pH->Clone();
        BiasX(pH3, 0.5*0.025*(x0 + x1));
        pHsum->Add(pH3);
        pH2->SetFillColorAlpha(kRed, 0.25);
        pH2->SetLineColor(kWhite);

        // draw histograms together
        sprintf(name, "%s_%s_DtPeak_%i", FC.c_str(), Setup.c_str(), i+1);
        pCDtPeak[i] = new TCanvas(name, name, 200, 10, 700, 500);
        gPad->SetTicks(1, 1);
        pH1->SetStats(0);
        pH1->GetXaxis()->SetLabelSize(0.05);
        pH1->GetXaxis()->SetTitleSize(0.07);
        pH1->GetXaxis()->SetTitleOffset(0.7);
        pH1->GetYaxis()->SetLabelSize(0.05);
        pH1->GetYaxis()->SetTitleSize(0.07);
        pH1->GetYaxis()->SetTitleOffset(0.7);
        pH1->Draw();
        pH2->Draw("same");
        pH1->Draw("same");

        // draw number
        Double_t x = 0, SF = 0, NIF = 0, DSF = 0, DNIF = 0; // prepare number drawing
        sprintf(name, "/%s/TimeDiff/Signal/%s/nNIF", FC.c_str(), Setup.c_str());
        TGraphErrors* gNIF = (TGraphErrors*)fCom->Get(name); // get NIF counts
        gNIF->GetPoint(i, x, NIF);
        DNIF = gNIF->GetErrorY(i);
        sprintf(message, "(n, f): %i #pm %i", (int)(NIF+0.5), (int)(DNIF+0.5));
        TLatex* tNIF = new TLatex();
        tNIF->SetNDC();
        tNIF->DrawLatex(0.6, 0.8, message);

        //// Draw original Dt spectrum with underground lines
        sprintf(name, "%s_%s_DtUg_%i", FC.c_str(), Setup.c_str(), i+1);
        pCDtUg[i] = new TCanvas(name, name, 200, 10, 700, 500);
        gPad->SetTicks(1,1);
        pH->SetStats(0);
        pH->GetXaxis()->SetRangeUser(62000, 80000);
        pH->Draw();

//        fUg->GetPoint(0, UgX[0], UgY[0]);
//        fUgP->GetPoint(0, UgX[1], UgY[1]);
//        TGraph* g0 = new TGraph(2, UgX, UgY);
//        g0->SetLineColor(kRed);
//        g0->SetLineWidth(lw);
//        g0->Draw("same");
//        fUgP->GetPoint(1, UgX[0], UgY[0]);
//        fUg->GetPoint(1, UgX[1], UgY[1]);
//        TGraph* g1 = new TGraph(2, UgX, UgY);
//        g1->SetLineColor(kRed);
//        g1->SetLineWidth(lw);
//        g1->Draw("same");

        SF = pH->Integral() - NIF;
        DSF = sqrt(pH->Integral() + DNIF*DNIF);
        sprintf(message, "Untergrund: %i #pm %i", (int)(SF+0.5), (int)(DSF+0.5));
        TLatex* tSF = new TLatex();
        tSF->SetNDC();
        tSF->DrawLatex(0.2, 0.2, message);
    }

    sprintf(name, "%s_%s_DtUg_all", FC.c_str(), Setup.c_str());
    TCanvas* pCsum = new TCanvas(name, name, 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    pHsum->SetNameTitle(name, "Time difference Sum over all deposits; TDC channels from maximum; counts");
    pHsum->SetStats(0);
    pHsum->GetXaxis()->SetTitleSize(0.05);
    pHsum->GetYaxis()->SetTitleSize(0.05);
    pHsum->Rebin(20);
    pHsum->SetAxisRange(-5000, 13000, "X");
    pHsum->Draw();
}


TH1F* GetDtRate(TFile* f, string FC, string Setup, Int_t channel)
{
    char name[64] = "";
    char title[64] = "";
    Double_t chMin, chMax, tMin, tMax;
    TH1D *pHtLive = (TH1D*) f->Get("/Histograms/Raw/Scaler/Rates/H1RawRate_47");
    Double_t t_live = pHtLive->Integral();
//    Double_t ChPerBin;
    Double_t TimePerCh = 25.E-03; // 25 ps in ns

    // Open PH-gated TimeDiff spectrum
    sprintf(name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", channel+1);
    TH1I *pHDtG = (TH1I*) f->Get(name);
    Int_t nbins = pHDtG->GetXaxis()->GetNbins();
    chMin = pHDtG->GetXaxis()->GetBinLowEdge(1);
    chMax = pHDtG->GetXaxis()->GetBinUpEdge(nbins);
//        ChPerBin = pHDtG->GetXaxis()->GetBinWidth(0);

    // Create rate over time histogram
    sprintf(name, "%s_%s_DtG_%i", FC.c_str(), Setup.c_str(), channel+1);
    sprintf(title, "; #font[12]{t} / ns; fission rate / 1/s");
    tMin = TimePerCh * chMin;
    tMax = TimePerCh * chMax;
    TH1F *pHDtRate = new TH1F(name, title, nbins, tMin, tMax);
    pHDtRate->GetXaxis()->SetRangeUser(1550, 2000);
    pHDtRate->GetXaxis()->SetTitleSize(0.05);
    pHDtRate->GetYaxis()->SetTitleSize(0.05);
    pHDtRate->GetYaxis()->SetTitleOffset(1.0);
    pHDtRate->SetStats(0);
    pHDtRate->SetLineWidth(lw);

    // Fill histogram
    Int_t counts;
    Double_t rate;
    for (int j = 0; j < nbins + 2; j++)
    {
        counts = pHDtG->GetBinContent(j);
        rate = 1.0 * counts / t_live;
        pHDtRate->SetBinContent(j, rate);
    }
    return pHDtRate;
}


void DrawDtRate(TFile* f, string FC, string Setup)
{
    char name[64] = "";
    TCanvas *pC[8];

    for (int i = 0; i < 8; i++)
    {
        // Create
        TH1F* pHDtRate = GetDtRate(f, FC, Setup, i);

        // Draw
        sprintf(name, "%s_%s_DtG_%i", FC.c_str(), Setup.c_str(), i+1);
        pC[i] = new TCanvas(name, name, 200, 10, 700, 500);
        gPad->SetTicks(1, 1);
        TLegend *l = new TLegend(0.65, 0.75, 0.85, 0.85, "");
        sprintf(name, "%s channel %i", FC.c_str(), i+1);
        l->AddEntry(pHDtRate, name);
        pHDtRate->Draw();
        l->Draw();
    }
}


void DrawDtFC(TFile *fNIF, TFile *fSB, string FC, Int_t RebinFactor)
{
    char name[64] = "";
    TCanvas *pC[8];

    for (int i = 0; i < 8; i++)
    {
        // Create
        TH1F* pHDtNIF = GetDtRate(fNIF, FC, "NIF", i);
        pHDtNIF->Rebin(RebinFactor);
        pHDtNIF->GetXaxis()->SetRangeUser(1550, 2000);
        pHDtNIF->SetLineColor(kGray);
        TH1F* pHDtSB = GetDtRate(fSB, FC, "SB", i);
        pHDtSB->Rebin(RebinFactor);
        pHDtSB->SetLineColor(kBlue);

        // Draw
        sprintf(name, "%s_DtComp_%i", FC.c_str(), i+1);
        pC[i] = new TCanvas(name, name, 200, 10, 700, 500);
        gPad->SetTicks(1, 1);
        TLegend *l = new TLegend(0.7, 0.8, 0.95, 0.95, FC.c_str());
        l->AddEntry(pHDtNIF, "Open geometry");
        l->AddEntry(pHDtSB, "Shadow Cone");
        pHDtNIF->Draw();
//        pHDtNIF->Draw("same");
        pHDtSB->Draw("same");
        l->Draw();
    }
}


void DrawDtComp(TFile *fPu, TFile *fU, string Setup, Int_t RebinFactor)
{
    char name[128] = "";
    TCanvas *pC[8];

    string SFrates[] = {"4.3268#pm0.0028",
                        "3.9011+-0.0027",
                        "3.4881+-0.0025",
                        "3.3631+-0.0025",
                        "3.5231+-0.0025",
                        "3.5272+-0.0025",
                        "3.2923+-0.0024",
                        "4.2637#pm0.0028"};

    for (int i = 0; i < 8; i++)
    {
        // Create
        TH1F* pHDtPu = GetDtRate(fPu, "PuFC", Setup, i);
        pHDtPu->Rebin(RebinFactor);
        pHDtPu->Scale(1000);
        sprintf(name, "Kanal %i; #font[12]{ToF} [ns]; #frac{#font[12]{SR}}{ns} [10^{-3} s^{-1}]", i+1);
        pHDtPu->SetTitle(name);
        pHDtPu->SetTitleSize(0.08f, "t");
        pHDtPu->GetXaxis()->SetRangeUser(1550, 2000);
        pHDtPu->GetXaxis()->SetLabelSize(0.06);
        pHDtPu->GetXaxis()->SetTitleSize(0.08);
        pHDtPu->GetXaxis()->SetTitleOffset(0.8);
        pHDtPu->GetYaxis()->SetLabelSize(0.06);
        pHDtPu->GetYaxis()->SetTitleSize(0.08);
        pHDtPu->GetYaxis()->SetTitleOffset(1);
        pHDtPu->SetLineColor(kMagenta);

        TH1F* pHDtU = GetDtRate(fU, "UFC", Setup, i);
        pHDtU->Rebin(RebinFactor);
        pHDtU->Scale(1000);
        pHDtU->SetLineColor(kBlue);

        // Draw
        sprintf(name, "%s_DtComp_%i", Setup.c_str(), i+1);
        pC[i] = new TCanvas(name, name, 200, 10, 700, 500);
        gPad->SetTicks(1, 1);

        TLegend *l = new TLegend(0.7, 0.8, 1.0, 1.0);
        l->AddEntry(pHDtPu, "^{242}Pu");
        l->AddEntry(pHDtU, "^{235}U");
        pHDtPu->Draw();
        pHDtU->Draw("same");
        l->Draw();

        TLatex *t = new TLatex();
        t->SetNDC();
//        sprintf(name, "Pu #font[12]{SR}: %s s^{-1}", SFrates[i].c_str());
        t->DrawLatex(0.2, 0.3, "Pu #font[12]{SR}: (4.2637#pm0.0028) s^{-1}");
        pC[i]->Modified();
        pC[i]->Update();
    }
}


void DrawDtComp(TFile *fPuNIF, TFile *fPuSB, TFile *fUNIF, TFile *fUSB, Int_t RebinFactor)
{
    char name[128] = "";
    TCanvas *pC[8];
    Double_t PuFg = 0.0412;
    Double_t PuBg = 0.0064;
    Double_t UFg = 0.1054;
    Double_t UBg = 0.0057;

    for (int i = 7; i < 8; i++)
    {
        // Create
        TH1F* pHDtPu = GetDtRate(fPuNIF, "PuFC", "NIF", i);
        pHDtPu->Rebin(RebinFactor);
        pHDtPu->Scale(1000);
        sprintf(name, "Kanal %i; #font[12]{ToF} [ns]; #frac{#font[12]{SR}}{ns} [10^{-3} s^{-1}]", i+1);
        pHDtPu->SetTitle(name);
        pHDtPu->SetTitleSize(0.08f, "t");
        pHDtPu->GetXaxis()->SetRangeUser(1550, 2000);
        pHDtPu->GetXaxis()->SetLabelSize(0.06);
        pHDtPu->GetXaxis()->SetTitleSize(0.08);
        pHDtPu->GetXaxis()->SetTitleOffset(0.8);
        pHDtPu->GetYaxis()->SetLabelSize(0.06);
        pHDtPu->GetYaxis()->SetTitleSize(0.08);
        pHDtPu->GetYaxis()->SetTitleOffset(0.6);
        pHDtPu->SetLineColor(kGray);

        TH1F* pHDtU = GetDtRate(fUNIF, "UFC", "NIF", i);
        pHDtU->Rebin(RebinFactor);
        pHDtU->Scale(1000);
        pHDtU->SetLineColor(kGray);

        TH1F* pHDtPuSB = GetDtRate(fPuSB, "PuFC", "SB", i);
        pHDtPuSB->Rebin(RebinFactor);
        pHDtPuSB->Scale(1000);
        pHDtPuSB->SetLineColor(kMagenta);

        TH1F* pHDtUSB = GetDtRate(fUSB, "UFC", "SB", i);
        pHDtUSB->Rebin(RebinFactor);
        pHDtUSB->Scale(1000);
        pHDtUSB->SetLineColor(kBlue);

        // Draw
        sprintf(name, "DtComp_%i", i+1);
        pC[i] = new TCanvas(name, name, 200, 10, 700, 500);
        gPad->SetTicks(1, 1);

        TLegend *l = new TLegend(0.7, 0.8, 1.0, 1.0);
        l->AddEntry(pHDtPuSB, "^{242}Pu");
        l->AddEntry(pHDtUSB, "^{235}U");
        pHDtPu->Draw();
        pHDtU->Draw("same");
        pHDtPuSB->Draw("same");
        pHDtUSB->Draw("same");
        l->Draw();

        TLatex *t1 = new TLatex();
        t1->SetNDC();
        sprintf(name, "%.4f s^{-1}", PuFg);
        t1->SetTextColor(14);
        t1->DrawLatex(0.4, 0.8, name);

        TLatex *t2 = new TLatex();
        t2->SetNDC();
        sprintf(name, "%.4f s^{-1}", PuBg);
        t2->SetTextColor(kMagenta);
        t2->DrawLatex(0.4, 0.6, name);

        TLatex *t3 = new TLatex();
        t3->SetNDC();
        sprintf(name, "%.4f s^{-1}", UFg);
        t3->SetTextColor(14);
        t3->DrawLatex(0.6, 0.4, name);

        TLatex *t4 = new TLatex();
        t4->SetNDC();
        sprintf(name, "%.4f s^{-1}", UBg);
        t4->SetTextColor(kBlue);
        t4->DrawLatex(0.6, 0.2, name);

        pC[i]->Modified();
        pC[i]->Update();
    }
}


void DrawDtComp4(TFile *fPuNIF, TFile *fPuSB, TFile *fUNIF, TFile *fUSB, Int_t RebinFactor)
{
    char name[64] = "";
    TCanvas *pC[8];
//    Int_t PuFg[] = {0.03483, 0.02902, 0.0258, 0.0265, 0.0251, 0.0306, 0.0289, 0.0411};
//    Int_t PuBg[] = {61, -18, -54, 23, 8, 98, 24, 147};
//    Int_t UFg[] = {5099, 5659, 5925, 5793, 5781, 6180, 6241, 6226};
//    Int_t UBg[] = {91, 41, 74, 75, 57, 74, 62, 62};

    Double_t PuFg = 0.0412;
    Double_t PuBg = 0.0064;
    Double_t UFg = 0.1054;
    Double_t UBg = 0.0057;
    TH1F* pHDtPuNIF[8];
    TH1F* pHDtPuSB[8];
    TH1F* pHDtUNIF[8];
    TH1F* pHDtUSB[8];

    for (int i = 0; i < 8; i++)
    {
        // Create
        pHDtPuNIF[i] = GetDtRate(fPuNIF, "PuFC", "NIF", i);
        pHDtPuNIF[i]->Rebin(RebinFactor);
        pHDtPuNIF[i]->Scale(1000);
        sprintf(name, "Kanal %i; #font[12]{ToF} [ns]; #frac{#font[12]{SR}}{ns} [10^{-3} s^{-1}]", i+1);
        pHDtPuNIF[i]->SetTitle(name);
        pHDtPuNIF[i]->SetTitleSize(0.08f, "t");
        pHDtPuNIF[i]->GetXaxis()->SetRangeUser(1550, 2000);
        pHDtPuNIF[i]->GetXaxis()->SetLabelSize(0.06);
        pHDtPuNIF[i]->GetXaxis()->SetTitleSize(0.08);
        pHDtPuNIF[i]->GetXaxis()->SetTitleOffset(0.8);
        pHDtPuNIF[i]->GetYaxis()->SetLabelSize(0.06);
        pHDtPuNIF[i]->GetYaxis()->SetTitleSize(0.08);
        pHDtPuNIF[i]->GetYaxis()->SetTitleOffset(0.8);
        pHDtPuNIF[i]->SetLineColor(14);

        pHDtPuSB[i] = GetDtRate(fPuSB, "PuFC", "SB", i);
        pHDtPuSB[i]->SetName("Dia2");
        pHDtPuSB[i]->Rebin(RebinFactor);
        pHDtPuSB[i]->Scale(1000);
        pHDtPuSB[i]->SetLineColor(kMagenta);

        pHDtUNIF[i] = GetDtRate(fUNIF, "UFC", "NIF", i);
        pHDtUNIF[i]->SetName("Dia3");
        pHDtUNIF[i]->Scale(1000);
        pHDtUNIF[i]->Rebin(RebinFactor);
        pHDtUNIF[i]->SetLineColor(14);

        pHDtUSB[i] = GetDtRate(fUSB, "UFC", "SB", i);
        pHDtUSB[i]->SetName("Dia4");
        pHDtUSB[i]->Scale(1000);
        pHDtUSB[i]->Rebin(RebinFactor);
        pHDtUSB[i]->SetLineColor(kBlue);

        // Draw
        sprintf(name, "DtComp_%i", i+1);
        pC[i] = new TCanvas(name, name, 200, 10, 700, 500);
        gPad->SetTicks(1, 1);
        TLegend *l = new TLegend(0.7, 0.75, 1.0, 1.0);
        l->AddEntry(pHDtPuNIF[i], "Offen");
        l->AddEntry(pHDtPuSB[i], "^{242}Pu, Schattenkegel");
        l->AddEntry(pHDtUSB[i], "^{235}U, Schattenkegel");
        pHDtPuNIF[i]->Draw();
        pHDtUNIF[i]->Draw("same");
        pHDtPuSB[i]->Draw("same");
        pHDtUSB[i]->Draw("same");
        l->Draw();

        TLatex *t1 = new TLatex();
        t1->SetNDC();
        sprintf(name, "%.4f s^{-1}", PuFg);
        t1->SetTextColor(14);
        t1->DrawLatex(0.4, 0.8, name);

        TLatex *t2 = new TLatex();
        t2->SetNDC();
        sprintf(name, "%.4f s^{-1}", PuBg);
        t2->SetTextColor(kMagenta);
        t2->DrawLatex(0.4, 0.6, name);

        TLatex *t3 = new TLatex();
        t3->SetNDC();
        sprintf(name, "%.4f s^{-1}", UFg);
        t3->SetTextColor(14);
        t3->DrawLatex(0.6, 0.4, name);

        TLatex *t4 = new TLatex();
        t4->SetNDC();
        sprintf(name, "%.4f s^{-1}", UBg);
        t4->SetTextColor(kBlue);
        t4->DrawLatex(0.6, 0.2, name);
    }
}


void DrawDtPeak(TFile* f, TFile* fCom, string FC, string Setup)
{   // Draw the Dt-peak fits into Dt spectrum
    if(!strcmp(Setup.c_str(), "SF")) // if Underground measurement
    {
        cout << "No TimeDiff peak for underground!" << endl;
        return;
    }
    char name[64] = "";
    char hname[64] = "";
    char message[64] = "";
    TCanvas* pCDtFit[8];
    // open integration results
//    sprintf(name, "/%s/TimeDiff/Signal/%s/nNIFfit", FC.c_str(), Setup.c_str());
//    TGraphErrors* gIntegral = (TGraphErrors*)fCom->Get(name);
    for (int i = 0; i < 8; i++)
    {
        // Create canvas
        sprintf(name, "%s_%s_DtPeakFit_%i", FC.c_str(), Setup.c_str(), i+1);
        pCDtFit[i] = new TCanvas(name, name, 200, 10, 700, 500);
        gPad->SetTicks(1, 1);
        // open Dt histogram nr i
        sprintf(hname, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i+1);
        TH1I* pH = (TH1I*)f->Get(hname);
        // open corresponding peak fit
        sprintf(name, "/%s/TimeDiff/Signal/%s/Fit/f%s%sDtPeak_%i", FC.c_str(), Setup.c_str(), FC.c_str(), Setup.c_str(), i+1);
        TF1* fP = (TF1*)fCom->Get(name);

        // manipulate graphs
        Int_t yFactor = 10;
        pH->Rebin(yFactor);
        TGraph* gr = Rebin(fP, 1000, yFactor);
        pH->SetAxisRange(62000, 80000, "X");
        pH->SetLineWidth(lw);
        gr->SetLineColor(kRed);
        gr->SetLineWidth(lw);
        // Draw
        pH->Draw();
        gr->Draw("same");
        // Draw integral number...
        Double_t x = 0, content = 0, Dcontent = 0; // prepare number drawing
        sprintf(name, "/%s/TimeDiff/Signal/%s/nNIFfit", FC.c_str(), Setup.c_str());
        TGraphErrors* gNIF = (TGraphErrors*)fCom->Get(name); // get NIF counts
        gNIF->GetPoint(i, x, content);
        Dcontent = gNIF->GetErrorY(i);
        sprintf(message, "peak content: %i +- %i", (int)(content+0.5), (int)(Dcontent+0.5));
        TText* tNIF = new TText();
        tNIF->SetNDC();
        tNIF->DrawText(0.3, 0.8, message);
//        cout << "checkpoint " << i+1 << endl;
    }
}


void BiasX(TGraphErrors* pG, Double_t dx)
{
    Double_t x = 0, y = 0;
    for (int i = 0; i < 8; i++)
    {
        pG->GetPoint(i, x, y);
        pG->SetPoint(i, x + dx, y);
    }
}


void DrawSFRate(TFile* fNIF, TFile* fSB, TFile* fSF, TFile* fCom)
{
    /// SF rate
    TCanvas* c7 = new TCanvas("SFRate", "TimeDiff Counting results", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg1 = new TMultiGraph();
    mg1->SetTitle("PuFC SF detection rate; Deposit; events per second");
    TLegend* l1 = new TLegend(0.6, 0.2, 0.85, 0.70, "Measurement");

    TGraphErrors* g1 = (TGraphErrors*)fNIF->Get("/Analysis/TimeDiff/SF_Rate");
    BiasX(g1, -0.15);
    g1->SetLineWidth(lw);
    g1->SetLineColor(kRed);
    mg1->Add(g1);
    l1->AddEntry(g1, "Open geometry", "lp");

    TGraphErrors* g2 = (TGraphErrors*)fSB->Get("/Analysis/TimeDiff/SF_Rate");
    BiasX(g2, -0.05);
    g2->SetLineWidth(lw);
    g2->SetLineColor(kGreen);
    mg1->Add(g2);
    l1->AddEntry(g2, "Shadow cone", "lp");

    TGraphErrors* g3 = (TGraphErrors*)fSF->Get("/Analysis/TimeDiff/SF_Rate");
    BiasX(g3, 0.05);
    g3->SetLineWidth(lw);
    g3->SetLineColor(kBlue);
    mg1->Add(g3);
    l1->AddEntry(g3, "Beam off", "lp");

    TGraphErrors* g7 = (TGraphErrors*)fCom->Get("/Analysis/PuFC/SF/SF_Rate");
    BiasX(g7, 0.15);
    g7->SetLineWidth(lw);
    g7->SetLineColor(kBlack);
    mg1->Add(g7);
    l1->AddEntry(g7, "Combined", "lp");

    mg1->Draw("AP");
    mg1->GetYaxis()->SetRangeUser(0, 7.5);
    l1->Draw();
    c7->Modified();
    c7->Update();
}


void DrawNPu(TFile* f)
{
    TCanvas* c2 = new TCanvas("NPuEff", "NPuEff", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg2 = new TMultiGraph();
    mg2->SetTitle("Effective number of Pu-272 atoms; Deposit; N");
    TLegend* l2 = new TLegend(0.6, 0.2, 0.85, 0.70, "Evaluated");

    TGraphErrors* g5 = (TGraphErrors*)f->Get("/Analysis/PuFC/SF/NPuEff");
    g5->SetLineColor(kBlue);
    g5->SetLineWidth(lw);
    mg2->Add(g5);
    l2->AddEntry(g5, "Spontaneaus fission", "lp");

    TGraphErrors* g6 = (TGraphErrors*)f->Get("/Analysis/Evaluation/NPuLit");
    g5->SetLineColor(kRed);
    g6->SetLineWidth(lw);
    mg2->Add(g6);
    l2->AddEntry(g6, "Neutron-induced fission", "lp");

    mg2->Draw("AP");
    mg2->GetYaxis()->SetRangeUser(0, 2.E19);
    l2->Draw();
    c2->Modified();
    c2->Update();

    TCanvas* c3 = new TCanvas("NPu", "NPu", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg3 = new TMultiGraph();
    mg3->SetTitle("Number of Pu-272 atoms using internal efficiency; Deposit; N");

    TGraphErrors* g7 = (TGraphErrors*)f->Get("/Analysis/PuFC/SF/NPu");
    g7->SetLineWidth(lw);
    mg3->Add(g7);

    mg3->Draw("AP");
    mg3->GetYaxis()->SetRangeUser(0, 2.E19);
    c3->Modified();
    c3->Update();
}


void DrawQDCmin(TFile* fNIF, TFile* fSB, TFile* fSF, TFile* fUNIF, TFile* fUSB, TFile* fCom)
{
    /// QDC extrema fit results
    TCanvas* c8 = new TCanvas("PuFC_QDCfit", "Common minimum fit results for PuFC", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg2 = new TMultiGraph();
    mg2->SetTitle("Common minimum fit results for PuFC; Deposit; QDC channel (low gain)");
    TLegend* l2 = new TLegend(0.6, 0.2, 0.85, 0.70, "Setup");

    TGraph* g6 = (TGraph*)fNIF->Get("/Analysis/QDC/gCut");
    g6->SetMarkerColor(kRed);
    mg2->Add(g6);
    l2->AddEntry(g6, "Open geometry", "lp");

    TGraph* g8 = (TGraph*)fSB->Get("/Analysis/QDC/gCut");
    g8->SetMarkerColor(kGreen);
    mg2->Add(g8);
    l2->AddEntry(g8, "Shadow bar", "lp");

    TGraph* g7 = (TGraph*)fSF->Get("/Analysis/QDC/gCut");
    g7->SetMarkerColor(kBlue);
    mg2->Add(g7);
    l2->AddEntry(g7, "Beam off", "lp");

    TGraphErrors* g5 = (TGraphErrors*)fCom->Get("/Analysis/PuFC/QDC/gCut");
    BiasX(g5, +0.1);
    g5->SetLineWidth(lw);
    g5->SetLineColor(kBlack);
    mg2->Add(g5);
    l2->AddEntry(g5, "Combined", "lp");

    mg2->Draw("AP");
    l2->Draw();
    c8->Modified();
    c8->Update();

    TCanvas* c9 = new TCanvas("PuFC_RelCut", "Relative minimum positions", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg3 = new TMultiGraph();
    mg3->SetTitle("Relative minimum positions; Deposit; Ratio");
    TLegend* l3 = new TLegend(0.6, 0.2, 0.85, 0.70, "Setup");

    TGraph* g9 = (TGraph*)fNIF->Get("/Analysis/QDC/gRelCut");
    g9->SetMarkerColor(kRed);
    mg3->Add(g9);
    l3->AddEntry(g9, "PuFC NIF", "lp");

    TGraph* g0 = (TGraph*)fSB->Get("/Analysis/QDC/gRelCut");
    g0->SetMarkerColor(kGreen);
    mg3->Add(g0);
    l3->AddEntry(g0, "PuFC SB", "lp");

    TGraph* g1 = (TGraph*)fSF->Get("/Analysis/QDC/gRelCut");
    g1->SetMarkerColor(kBlue);
    mg3->Add(g1);
    l3->AddEntry(g1, "PuFC SF", "lp");

    TGraph* g2 = (TGraph*)fUNIF->Get("/Analysis/QDC/gRelCut");
    g2->SetMarkerColor(kMagenta);
    mg3->Add(g2);
    l3->AddEntry(g2, "UFC NIF", "lp");

    TGraph* g3 = (TGraph*)fUSB->Get("/Analysis/QDC/gRelCut");
    g3->SetMarkerColor(kOrange);
    mg3->Add(g3);
    l3->AddEntry(g3, "UFC SB", "lp");

    mg3->Draw("AP");
    mg3->GetYaxis()->SetRangeUser(0.0, 1.0);
    l3->Draw();
    c9->Modified();
    c9->Update();
}


void DrawIntEff(TFile* fNIF, TFile* fSB, TFile* fSF, TFile* fCom)
{
    /// Internal efficiencies
    TCanvas* c1 = new TCanvas("PuFC_eInt", "Internal efficiencies", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg1 = new TMultiGraph();
    mg1->SetTitle("Efficiencies of PuFC Deposits (internal); Deposit; Ratio");
    TLegend* l1 = new TLegend(0.6, 0.2, 0.85, 0.70, "Setup");

    TGraphErrors* g1 = (TGraphErrors*)fNIF->Get("/Analysis/QDC/eInt");
    BiasX(g1, -0.15);
    g1->SetLineWidth(lw);
    g1->SetLineColor(kRed);
    mg1->Add(g1);
    l1->AddEntry(g1, "Open geometry", "lp");

    TGraphErrors* g2 = (TGraphErrors*)fSB->Get("/Analysis/QDC/eInt");
    BiasX(g2, -0.05);
    g2->SetLineWidth(lw);
    g2->SetLineColor(kGreen);
    mg1->Add(g2);
    l1->AddEntry(g2, "Shadow bar", "lp");

    TGraphErrors* g3 = (TGraphErrors*)fSF->Get("/Analysis/QDC/eInt");
    BiasX(g3, 0.05);
    g3->SetLineWidth(lw);
    g3->SetLineColor(kBlue);
    mg1->Add(g3);
    l1->AddEntry(g3, "Beam off", "lp");

    mg1->Draw("AP");
    mg1->GetYaxis()->SetRangeUser(0.0, 1.0);
    l1->Draw();
    c1->Modified();
    c1->Update();
}


void DrawInScat(TFile* fCom)
{
    //// percents
    TCanvas* c1 = new TCanvas("InScat", "In-scattering part", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg1 = new TMultiGraph();
    mg1->SetTitle("In-scattering part for fission chambers; Deposit; In-scattering / %");
    TLegend* l1 = new TLegend(0.6, 0.2, 0.85, 0.70, "Chamber");

    TGraphErrors* g0 = (TGraphErrors*)fCom->Get("/PuFC/TimeDiff/gScat");
    BiasX(g0, -0.1);
    g0->SetLineWidth(lw);
    g0->SetLineColor(kGreen);
    mg1->Add(g0);
    l1->AddEntry(g0, "PuFC", "lp");

    TGraphErrors* g1 = (TGraphErrors*)fCom->Get("/UFC/TimeDiff/gScat");
    BiasX(g1, +0.1);
    g1->SetLineWidth(lw);
    g1->SetLineColor(kBlue);
    mg1->Add(g1);
    l1->AddEntry(g1, "UFC", "lp");

    mg1->Draw("AP");
    mg1->GetYaxis()->SetRangeUser(0.0, 10.0);
    l1->Draw();
    c1->Modified();
    c1->Update();

    //// Rates normed to neutron monitor
    TCanvas* c2 = new TCanvas("PuFC_NIFtoMon", "PuFC_NIFtoMon", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg2 = new TMultiGraph();
    mg2->SetTitle("NIF detection rates normed to monitor; Deposit; NIF events / monitor counts");
    TLegend* l2 = new TLegend(0.6, 0.2, 0.85, 0.70, "Setup");

    TGraphErrors* g2 = (TGraphErrors*)fCom->Get("/PuFC/TimeDiff/Signal/NIF/NIFtoMon");
    BiasX(g2, -0.1);
    g2->SetLineWidth(lw);
    g2->SetLineColor(kBlue);
    mg2->Add(g2);
    l2->AddEntry(g2, "Open", "lp");

    TGraphErrors* g3 = (TGraphErrors*)fCom->Get("/PuFC/TimeDiff/Signal/SB/NIFtoMon");
    BiasX(g3, +0.1);
    g3->SetLineWidth(lw);
    g3->SetLineColor(kGray);
    mg2->Add(g3);
    l2->AddEntry(g3, "Shadow bar", "lp");

    mg2->Draw("AP");
    l2->Draw();
    c2->Modified();
    c2->Update();

    TCanvas* c3 = new TCanvas("UFC_NIFtoMon", "UFC_NIFtoMon", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg3 = new TMultiGraph();
    mg3->SetTitle("NIF detection rates normed to monitor; Deposit; NIF events / monitor counts");
    TLegend* l3 = new TLegend(0.6, 0.2, 0.85, 0.70, "Setup");

    TGraphErrors* g4 = (TGraphErrors*)fCom->Get("/UFC/TimeDiff/Signal/NIF/NIFtoMon");
    BiasX(g4, -0.1);
    g4->SetLineWidth(lw);
    g4->SetLineColor(kBlue);
    mg3->Add(g4);
    l3->AddEntry(g4, "Open", "lp");

    TGraphErrors* g5 = (TGraphErrors*)fCom->Get("/UFC/TimeDiff/Signal/SB/NIFtoMon");
    BiasX(g5, +0.1);
    g5->SetLineWidth(lw);
    g5->SetLineColor(kGray);
    mg3->Add(g5);
    l3->AddEntry(g5, "Shadow bar", "lp");

    mg3->Draw("AP");
    l3->Draw();
    c3->Modified();
    c3->Update();
}


void DrawEff(TFile* fCom)
{ // efficiency conclusion
    TCanvas* c2 = new TCanvas("EffPuFC", "Efficiencies of PuFC deposits", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg2 = new TMultiGraph();
    mg2->SetTitle("Efficiencies of PuFC Deposits; Deposit; Ratio");
    TLegend* l2 = new TLegend(0.6, 0.2, 0.85, 0.70, "Source");

    TGraphErrors* g4 = (TGraphErrors*)fCom->Get("/PuFC/Efficiency/eSF");
    BiasX(g4, -0.2);
    g4->SetLineWidth(lw);
    g4->SetLineColor(kRed);
    mg2->Add(g4);
    l2->AddEntry(g4, "Spontaneaus fission", "lp");

//    TGraphErrors* g5 = (TGraphErrors*)fCom->Get("/Analysis/PuFC/Efficiency/eRel");
//    BiasX(g5, -0.1);
//    g5->SetLineWidth(lw);
//    g5->SetLineColor(kGreen);
//    mg2->Add(g5);
//    l2->AddEntry(g5, "NIF:SF", "lp");

    TGraphErrors* g6 = (TGraphErrors*)fCom->Get("/Analysis/PuFC/Efficiency/eNIF");
    BiasX(g6, 0.0);
    g6->SetLineWidth(lw);
    g6->SetLineColor(kBlue);
    mg2->Add(g6);
    l2->AddEntry(g6, "Neutron-induced fission", "lp");

//    TGraphErrors* g7 = (TGraphErrors*)fCom->Get("/Analysis/PuFC/Efficiency/eSimG");
//    BiasX(g7, +0.1);
//    g7->SetLineWidth(lw);
//    g7->SetLineColor(kMagenta);
//    mg2->Add(g7);
//    l2->AddEntry(g7, "Simulation (Gayther)", "lp");

    TGraphErrors* g8 = (TGraphErrors*)fCom->Get("/Analysis/PuFC/Efficiency/eSimM");
    BiasX(g8, +0.2);
    g8->SetLineWidth(lw);
    g8->SetLineColor(kOrange);
    mg2->Add(g8);
    l2->AddEntry(g8, "Simulation (Minimum)", "lp");

    mg2->Draw("AP");
    mg2->GetYaxis()->SetRangeUser(0.0, 1.2);
    l2->Draw();
    c2->Modified();
    c2->Update();
}


void DrawEval(TFile* fCom)
{
    TCanvas* c1 = new TCanvas("XS_SF", "XS_SF", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg1 = new TMultiGraph();
    mg1->SetTitle("Cross section; Deposit; sigma in mm^2");

    TGraphErrors* g1 = (TGraphErrors*)fCom->Get("/Analysis/Evaluation/XS_SF");
    g1->SetLineWidth(lw);
    g1->SetLineColor(kBlue);
    mg1->Add(g1);

    mg1->Draw("AP");
    c1->Modified();
    c1->Update();

    TCanvas* c3 = new TCanvas("eps", "eps", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg3 = new TMultiGraph();
    mg3->SetTitle("Efficiencies; Deposit; Ratio");
    TLegend* l3 = new TLegend(0.6, 0.2, 0.85, 0.70, "Fission fragments");

    TGraphErrors* g7 = (TGraphErrors*)fCom->Get("/Analysis/PuFC/Efficiency/eNIF");
    BiasX(g7, -0.1);
    g7->SetLineColor(kBlue);
    g7->SetLineWidth(lw);
    mg3->Add(g7);
    l3->AddEntry(g7, "Neutron-induced");

    TGraphErrors* g8 = (TGraphErrors*)fCom->Get("/Analysis/PuFC/Efficiency/eSF");
    g8->SetLineColor(kGreen);
    g8->SetLineWidth(lw);
    mg3->Add(g8);
    l3->AddEntry(g8, "Spontaneaus (average)");

    TGraphErrors* g9 = (TGraphErrors*)fCom->Get("/Analysis/PuFC/Efficiency/eRel");
    BiasX(g9, +0.1);
    g9->SetLineColor(kRed);
    g9->SetLineWidth(lw);
    mg3->Add(g9);
    l3->AddEntry(g9, "Ratio NIF:SF");

    mg3->Draw("AP");
    mg3->GetYaxis()->SetRangeUser(0, 1.5);
    l3->Draw();
    c3->Modified();
    c3->Update();
}


void DrawScaledDt(TFile* fNIF, TFile* fSB, TFile* fUG)
{
    char name[64] = "";
    TCanvas* c[8];
    TH1F* hTlNIF = (TH1F*)fNIF->Get("/Histograms/Raw/Scaler/Rates/H1RawRate_47");
    Double_t tNIF = hTlNIF->Integral();
    TH1F* hTlSB = (TH1F*)fSB->Get("/Histograms/Raw/Scaler/Rates/H1RawRate_47");
    Double_t tSB = hTlSB->Integral();
    TH1F* hTlUG = (TH1F*)fUG->Get("/Histograms/Raw/Scaler/Rates/H1RawRate_47");
    Double_t tUG = hTlUG->Integral();
    cout << "Times: " << tNIF << ", " << tSB << ", " << tUG << endl;
    for (int i = 0; i < 1; i++)
    {
        sprintf(name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i+1);
        TH1I* h;
        h = (TH1I*)fNIF->Get(name);
        TH1F* hNIF = (TH1F*)(h->Clone());
        h = (TH1I*)fSB->Get(name);
        TH1F* hSB = (TH1F*)(h->Clone());
        h = (TH1I*)fUG->Get(name);
        TH1F* hUG = (TH1F*)(h->Clone());
        TLegend* l = new TLegend(0.6, 0.2, 0.85, 0.70, "Setup");
        hNIF->Scale(sqrt(tUG*tSB/tNIF));
        hSB->Scale(sqrt(tNIF*tUG/tSB));
        hUG->Scale(sqrt(tNIF*tSB/tUG));
        sprintf(name, "DtSc_%i", i+1);
        c[i] = new TCanvas(name, "TimeDiff scaled by live time", 200, 10, 700, 500);
        hNIF->SetLineColor(kRed);
        l->AddEntry(hNIF, "Open beam", "lp");
        hSB->SetLineColor(kGreen);
        l->AddEntry(hSB, "Shadow bar", "lp");
        hUG->SetLineColor(kBlue);
        l->AddEntry(hUG, "Beam off", "lp");
        hSB->SetAxisRange(62000, 80000, "X");
        hSB->SetAxisRange(0, 130000, "Y");
        hSB->Draw();
        hNIF->Draw("same");
        hUG->Draw("same");
        l->Draw();
        c[i]->Modified();
        c[i]->Update();
    }
}


void DrawScaledDt(TFile* fNIF, TFile* fSB)
{
    char name[64] = "";
    TCanvas* c[8];
    TH1F* hTlNIF = (TH1F*)fNIF->Get("/Histograms/Raw/Scaler/Rates/H1RawRate_47");
    Double_t tNIF = hTlNIF->Integral();
    TH1F* hTlSB = (TH1F*)fSB->Get("/Histograms/Raw/Scaler/Rates/H1RawRate_47");
    Double_t tSB = hTlSB->Integral();
    cout << "Times: " << tNIF << ", " << tSB << endl;
    for (int i = 0; i < 1; i++)
    {
        sprintf(name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i+1);
        TH1I* h;
        h = (TH1I*)fNIF->Get(name);
        TH1F* hNIF = (TH1F*)(h->Clone());
        h = (TH1I*)fSB->Get(name);
        TH1F* hSB = (TH1F*)(h->Clone());
        TLegend* l = new TLegend(0.6, 0.2, 0.85, 0.70, "Setup");
        hNIF->Scale(sqrt(tSB/tNIF));
        hSB->Scale(sqrt(tNIF/tSB));
        sprintf(name, "DtSc_%i", i+1);
        c[i] = new TCanvas(name, "TimeDiff scaled by live time", 200, 10, 700, 500);
        hNIF->SetLineColor(kRed);
        l->AddEntry(hNIF, "Open beam", "lp");
        hSB->SetLineColor(kGreen);
        l->AddEntry(hSB, "Shadow bar", "lp");
        hSB->SetAxisRange(62000, 80000, "X");
        hSB->SetAxisRange(0, 130000, "Y");
        hSB->Draw();
        hNIF->Draw("same");
        l->Draw();
        c[i]->Modified();
        c[i]->Update();
    }
}


void DrawCompareDt(TFile* fPu, TFile* fU, string Setup)
{
    char name[64] = "";
    TCanvas* pC[8];
    Double_t tlPu, tlU; // live times
    sprintf(name, "/Histograms/Raw/Scaler/Rates/H1RawRate_47");
    TH1D* pHtlPu = (TH1D*)fPu->Get(name);
    tlPu = pHtlPu->Integral();
    TH1D* pHtlU = (TH1D*)fU->Get(name);
    tlU = pHtlU->Integral();
    cout << "Times " << tlPu << " " << tlU << endl;
    for (int i = 0; i < 8; i++)
    {
        // open hists
        sprintf(name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i+1);
        TH1I* pHIPu = (TH1I*)fPu->Get(name);
        TH1I* pHIU = (TH1I*)fU->Get(name);
        // scale uranium plot
        sprintf(name, "ToF spectra, %s setup, channel %i", Setup.c_str(), i+1);
        Int_t n = pHIPu->GetNbinsX();
        Double_t xmin = pHIPu->GetBinLowEdge(1);
        Double_t xmax = pHIPu->GetBinLowEdge(n) + pHIPu->GetBinWidth(n);
        TH1F* pHFU = new TH1F(name, name, n, xmin, xmax);
        for (int bin = 0; bin <= n + 1; bin++)
            pHFU->SetBinContent(bin, pHIU->GetBinContent(bin) * tlPu / tlU);
        // manipulate hists
        pHIPu->GetXaxis()->SetRangeUser(62000, 80000);
        sprintf(name, "ToF spectra, %s setup, channel %i", Setup.c_str(), i+1);
        pHIPu->SetTitle(name);
        pHIPu->SetStats(0);
        pHIPu->GetXaxis()->SetTitleSize(0.05);
        pHIPu->GetYaxis()->SetTitleSize(0.05);
        pHIPu->SetLineColor(kGreen);
        pHFU->SetLineColor(kBlue);
        // create legend
        TLegend* l = new TLegend(0.6, 0.2, 0.85, 0.70, "Fission Chamber");
        l->AddEntry(pHIPu, "PuFC", "lp");
        l->AddEntry(pHFU, "UFC (scaled by live time)", "lp");
        // draw
        sprintf(name, "%s_DtComp_%i", Setup.c_str(), i+1);
        pC[i] = new TCanvas(name, name, 200, 10, 700, 500);
        gPad->SetTicks(1, 1);
        pHIPu->Draw();
        pHFU->Draw("same");
        l->Draw();
    }
}


//void TestRebin(TFile* f, string FC, string Setup, int i)
//{
//    char name[64] = "";
//    sprintf(name, "/Histograms/Raw/QDC/low/H1RawQDCl_%i", i+1);
//    TH1I* pH = (TH1I*)f->Get("");
//    sprintf(name, "/Analysis/QDC/Fit/f%s%sCut_%i", FC.c_str(), Setup.c_str(), i+1);
//    TF1* f = (TF1*)f->Get(name);
//}


void GetStability(string file, string FC, TMultiGraph* mg1, TMultiGraph* mg2, Int_t chStart = 0, Int_t chStop = 7, TLegend* l = 0, Int_t comp = 0)
{
    char name[64] = "";
    Int_t l0 = 159;
    Int_t l1[8], l2[8];
    Int_t l3 = 1649;
    if (!strcmp(FC.c_str(), "PuFC"))
    {
        Int_t L1[] = {522, 474, 465, 480, 480, 476, 471, 471};
        Int_t L2[] = {722, 674, 665, 680, 680, 676, 671, 671};
        for (Int_t i = 0; i < 8; i++)
        {
            l1[i] = L1[i];
            l2[i] = L2[i];
        }
    } else {
        Int_t L1[] = {1069, 1051, 1016, 1035, 1035, 1032, 1028, 1027};
        Int_t L2[] = {1269, 1251, 1216, 1235, 1235, 1232, 1228, 1227};
        for (Int_t i = 0; i < 8; i++)
        {
            l1[i] = L1[i];
            l2[i] = L2[i];
        }
    }
    sprintf(name, "/home/hoffma93/Programme/Go4nfis/offline/results/%s.root", file.c_str());
    TFile* f1 = TFile::Open(name);
    sprintf(name, "/Histograms/Raw/Scaler/Rates/H1RawRate_47");
    TH1D* h2 = (TH1D*)f1->Get(name);

    for (Int_t i = chStart; i <= chStop; i++)
    {
        sprintf(name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H2DtGvsTime_%i", i+1);
        TH2I* h1 = (TH2I*)f1->Get(name);
        Int_t tBins = h1->GetNbinsY();
        Int_t tStart = h1->GetYaxis()->GetBinLowEdge(1);
        Int_t tEnd = h1->GetYaxis()->GetBinLowEdge(tBins + 1);
//        cout << tStart << "-" << tEnd << endl;
        Int_t tStep = 900;
        Int_t tPoints = (tEnd - tStart) / tStep;
//        cout << tPoints << " << " << tBins << endl;
        Double_t X[tPoints];
        Double_t Xerr[tPoints];
        Double_t rNIF[tPoints];
        Double_t DrNIF[tPoints];
        Double_t rSF[tPoints];
        Double_t DrSF[tPoints];
        Double_t nNIF, DnNIF, nSF, DnSF;

        for (Int_t t = 0; t < tPoints; t++)
        {
            Int_t start = t * h2->GetNbinsX() / tPoints;
            Int_t stop = start + h2->GetNbinsX() / tPoints - 1;
            Double_t tLive = h2->Integral(start, stop);

            start = t * tBins / tPoints;
            stop = start + tBins / tPoints - 1;
            X[t] = tStart + t * tStep + tStep / 2 + 100*i;
            Xerr[t] = 0;
            if (tLive > 0)
            {
//                cout << t << " " << strcmp(file.c_str(), "SF") << endl;
                if (strcmp(file.c_str(), "SF"))
                { // if beam on
                    Double_t PeakInt = h1->Integral(l1[i], l2[i], start, stop);
                    Double_t UgInt = h1->Integral(l0, l1[i] - 1, start, stop) + h1->Integral(l2[i] + 1, l3, start, stop);
                    nNIF = PeakInt - (l2[i]-l1[i]+1.0) / (l3-l2[i]+l1[i]-l0) * UgInt;
                    DnNIF = sqrt( PeakInt + pow((l2[i]-l1[i]+1) / (l3-l2[i]+l1[i]-l0), 2) * UgInt );
                } else { // SF measurement - no neutrons - no (n,f)
                    nNIF = 0;
                    DnNIF = 0;
                }
                Double_t Int = h1->Integral(0, h1->GetNbinsX()+1, start, stop);
                nSF = Int - nNIF;//Int - h1->Integral(l1[i], l2[i], start, stop);//
                DnSF = sqrt( nSF );
                rNIF[t] = nNIF / tLive;
                DrNIF[t] = DnNIF / tLive;
                rSF[t] = nSF / tLive;
                DrSF[t] = DnSF / tLive;
//                cout << t << ", " << start << "-" << stop << ", " << X[t] - 1400000000 << "+-" << Xerr[t] << ", " << rSF[t] << "+-" << DrSF[t] << endl;
            } else {
                rNIF[t] = 0;
                DrNIF[t] = 0;
                rSF[t] = 0;
                DrSF[t] = 0;
            }
        }
        TGraphErrors* ge1 = new TGraphErrors(tPoints, X, rNIF, Xerr, DrNIF);
        sprintf(name, "%s_%s_%i", FC.c_str(), file.c_str(), i+1);
        ge1->SetName(name);
        ge1->SetLineColorAlpha(comp ? 2 : i+1, 0.5);
        ge1->SetMarkerColorAlpha(0, 0.0);
        ge1->SetLineWidth(1);
        mg1->Add(ge1);
        TGraphErrors* ge2 = new TGraphErrors(tPoints, X, rSF, Xerr, DrSF);
        sprintf(name, "%s_%s_%i", FC.c_str(), file.c_str(), i+1);
        ge2->SetName(name);
        ge2->SetLineColorAlpha(comp ? 3 : i+1, 0.5);
        ge2->SetMarkerColorAlpha(0, 0.0);
        ge2->SetLineWidth(1);
        mg2->Add(ge2);
        if (l != 0)
        {
            if (comp == 1)
            {
                l->AddEntry(ge1, "(n,f)");
                l->AddEntry(ge2, "(s,f)");
            } else {
                sprintf(name, "%i", i+1);
                l->AddEntry(ge1, name);
            }
        }
    }
}


void DrawStability(Int_t i)
{
    char name[64] = "";
    TMultiGraph *mg1 = new TMultiGraph("mg_sr", "Spaltrate; #font[12]{t}; Rate [1/s]");
    sprintf(name, "Deposit %i", i+1);
    TLegend *l = new TLegend(0.8, 0.7, 1.0, 0.9, name);

    GetStability("UFC_NIF", "UFC", mg1, mg1, i, i, l, 1);
    GetStability("UFC_SB", "UFC", mg1, mg1, i, i, 0, 1);
    GetStability("NIF", "PuFC", mg1, mg1, i, i, 0, 1);
    GetStability("SB", "PuFC", mg1, mg1, i, i, 0, 1);
    GetStability("SF", "PuFC", mg1, mg1, i, i, 0, 1);

    sprintf(name, "Stab_%i", i+1);
    TCanvas *c = new TCanvas(name, "(n,f) Stability", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    mg1->Draw("AP");
    mg1->GetXaxis()->SetTimeDisplay(1);
    mg1->GetXaxis()->SetTimeFormat("%d.%m. %H:%M%F1970-01-01");
    l->Draw();
    c->Modified();
    c->Update();
}


void DrawStability()
{
    TMultiGraph *mg1 = new TMultiGraph("mg_nf", "(n,f)-Rate; #font[12]{t}; Rate [1/s]");
    TMultiGraph *mg2 = new TMultiGraph("mg_sf", "UG-Rate; #font[12]{t}; Rate [1/s]");
    TLegend *l = new TLegend(0.9, 0.5, 1.0, 1.0, "Deposit");

    Int_t chStart = 0, chStop = 7;
    GetStability("UFC_NIF", "UFC", mg1, mg2, chStart, chStop, l, 0);
    GetStability("UFC_SB", "UFC", mg1, mg2, chStart, chStop);
//    GetStability("NIF", "PuFC", mg1, mg2, chStart, chStop);
//    GetStability("SB", "PuFC", mg1, mg2, chStart, chStop);
//    GetStability("SF", "PuFC", mg1, mg2, chStart, chStop);
//    GetStability("UFC_FG_MS20_2", "UFC", mg1, mg2, chStart, chStop, l, 0);
//    GetStability("UFC_FG_MS20_3", "UFC", mg1, mg2, chStart, chStop);
//    GetStability("UFC_FG_MS20_4", "UFC", mg1, mg2, chStart, chStop);
//    GetStability("UFC_BG_MS20_5", "UFC", mg1, mg2, chStart, chStop);
//    GetStability("UFC_FG_MS21_2", "UFC", mg1, mg2, chStart, chStop);
//    GetStability("UFC_FG_MS21_3", "UFC", mg1, mg2, chStart, chStop);
//    GetStability("UFC_BG_MS21_4", "UFC", mg1, mg2, chStart, chStop);
//    GetStability("PuFC_FG_MS4", "PuFC", mg1, mg2, chStart, chStop);
//    GetStability("PuFC_FG_MS5", "PuFC", mg1, mg2, chStart, chStop);
//    GetStability("PuFC_FG_MS6", "PuFC", mg1, mg2, chStart, chStop);
//    GetStability("PuFC_FG_MS6", "PuFC", mg1, mg2, chStart, chStop);
////    GetStability("PuFC_BG_MS8", "PuFC", mg1, mg2, chStart, chStop);
//    GetStability("PuFC_BG_MS9", "PuFC", mg1, mg2, chStart, chStop);
//    GetStability("PuFC_BG_MS10", "PuFC", mg1, mg2, chStart, chStop);
//    GetStability("PuFC_BG_MS11", "PuFC", mg1, mg2, chStart, chStop);

    TCanvas *c1 = new TCanvas("Stab1", "(n,f) Stability", 200, 10, 700, 500);
//    l->Draw();
    TPad *p1 = new TPad("p1", "(n,f)", 0.0, 0.5, 1.0, 1.0);
    gPad->SetTicks(1, 1);
    p1->Draw();
    p1->cd();
    mg1->Draw("AP");
    mg1->GetXaxis()->SetTimeDisplay(1);
    mg1->GetXaxis()->SetTimeFormat("%d.%m. %H:%M%F1970-01-01");
    c1->cd();
    TPad *p2 = new TPad("p2", "n field", 0.0, 0.0, 1.0, 0.5);
    p2->Draw();
    p2->cd();
    mg2->Draw("AP");
    mg2->GetXaxis()->SetTimeDisplay(1);
    mg2->GetXaxis()->SetTimeFormat("%d.%m. %H:%M%F1970-01-01 00:00:00s0");

    c1->Modified();
    c1->Update();

//    TCanvas *c2 = new TCanvas("Stab2", "(s,f) Stability", 200, 10, 700, 500);
//    gPad->SetTicks(1, 1);
//    l->Draw();
//    c2->Modified();
//    c2->Update();
}


void GetPeak(Int_t i, TH1I* pH, string FC, Int_t left, Int_t right, Double_t* nif, Double_t* Dnif)
{
    char name[64] = "";
    Int_t l0 = 159;
    Int_t l1;
    Int_t l2;
    Int_t l3 = 1649;
    if (!strcmp(FC.c_str(), "PuFC")) // if PuFC
    {
        Int_t m[] = {572, 524, 515, 530, 530, 526, 521, 521};
        l1 = m[i] - left;
        l2 = m[i] + right;
    } else {
        Int_t m[] = {1119, 1101, 1066, 1085, 1085, 1082, 1078, 1077};
        l1 = m[i] - left;
        l2 = m[i] + right;
    }
    Double_t PeakInt = pH->Integral(l1, l2);
    Double_t UgInt = pH->Integral(l0, l3) - PeakInt;
    *nif = PeakInt - (l2-l1+1.0) / (l3-l2+l1-l0) * UgInt;
    *Dnif = sqrt( PeakInt + pow((l2-l1+1) / (l3-l2+l1-l0), 2) * UgInt );
}


void DrawPeak(TFile* f, string FC)
{
    char name[64] = "";
    TMultiGraph* mg = new TMultiGraph("mgPeak", "ToF Peak Integration; Breite [ns]; Peak");
    sprintf(name, "%s Deposit", FC.c_str());
    TLegend *l = new TLegend(0.9, 0.5, 1.0, 1.0, name);
    TGraphErrors* ge[8];
    Int_t W[] = {20, 40, 60, 80, 100, 120, 140, 160, 180, 200,
                 220, 240, 260, 280, 300, 320, 340, 360, 380, 400,
                 420, 440, 460, 480, 500, 520, 540, 560, 580, 600
//                 620, 640, 660, 680, 700, 720, 740, 760, 780, 800,
//                 820, 840, 860, 880, 900, 920, 940, 960, 980, 1000
                }; // width samples
    Int_t n = sizeof(W) / sizeof(W[0]);
    Double_t X[n];
    Double_t Xerr[n];
    Double_t avNIF[n];
    Double_t DavNIF[n];
    for (Int_t j = 0; j < n; j++)
    {
        X[j] = W[j] / 4.0;
        Xerr[j] = 0;
        avNIF[j] = 0;
        DavNIF[j] = 0;
    }

    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i+1);
        TH1I* pH = (TH1I*)f->Get(name);
        Double_t NIF[n];
        Double_t DNIF[n];
        for (Int_t j = 0; j < n; j++)
        {
            GetPeak(i, pH, FC, 60, W[j], &NIF[j], &DNIF[j]);
            avNIF[j] += NIF[j] / 8;
            DavNIF[j] = sqrt( pow(DavNIF[j], 2) + pow(DNIF[j] / 8, 2) );
//            cout << W[j] << ", " << NIF[j] << "+-" << DNIF[j] << endl;
        }
        ge[i] = new TGraphErrors(n, X, NIF, Xerr, DNIF);
        sprintf(name, "gPeak_%i", i+1);
        ge[i]->SetName(name);
        ge[i]->SetLineColorAlpha(i+2, 0.5);
        ge[i]->SetLineWidth(1);
        ge[i]->SetMarkerColorAlpha(i+2, 0.5);
        ge[i]->SetMarkerStyle(20);
        ge[i]->SetMarkerSize(2);
        sprintf(name, "%i", i+1);
        l->AddEntry(ge[i], name);
        mg->Add(ge[i]);
    }
    TGraphErrors* av = new TGraphErrors(n, X, avNIF, Xerr, DavNIF);
    av->SetName("gPeak_av");
    av->SetLineWidth(2);
    av->SetLineColor(1);
    l->AddEntry(av, "av");
    mg->Add(av);

    TCanvas* c = new TCanvas("cPeak", "Test Integration Limits", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    mg->Draw("AP");
    l->Draw();
    c->Modified();
    c->Update();
}


void GetLRUg(Int_t i, TH1I* pH, string FC, Int_t left, Int_t right, Double_t* lrUg, Double_t* DlrUg)
{
    char name[64] = "";
    Int_t l0 = 159;
    Int_t l1;
    Int_t l2;
    Int_t l3 = 1649;
    Double_t lInt, rInt, lBins, rBins; // left/right integrals and interval length
    if (!strcmp(FC.c_str(), "PuFC")) // if PuFC
    {
        Int_t m[] = {572, 524, 515, 530, 530, 526, 521, 521};
        l1 = m[i] - left;
        l2 = m[i] + right;
        lInt = pH->Integral(l0, l1 - 1) + pH->Integral(l3 - 499, l3);
        rInt = pH->Integral(l2 + 1, l3 - 500);
        lBins = l1 - l0 + 500;
        rBins = l3 - 500 - l2;
    } else {
        Int_t m[] = {1119, 1101, 1066, 1085, 1085, 1082, 1078, 1077};
        l1 = m[i] - left;
        l2 = m[i] + right;
        lInt = pH->Integral(l0, l1 - 1);
        rInt = pH->Integral(l2 + 1, l3);
        lBins = l1 - l0;
        rBins = l3 - l2;
    }
//    if ()
    Double_t lUg = lInt / lBins;
    Double_t rUg = rInt / rBins;
    Double_t DlUg = sqrt(lInt) / lBins;
    Double_t DrUg = sqrt(rInt) / rBins;
//    cout << pH->GetBinLowEdge(l0) << " " << pH->GetBinLowEdge(l1) << endl;
//    TF1* fl = new TF1("lUg", "[0]", pH->GetBinLowEdge(l0), pH->GetBinLowEdge(l1));
//    pH->Fit("lUg", "LRQ0");
//    TF1* fr = new TF1("rUg", "[0]", pH->GetBinLowEdge(l2+1), pH->GetBinLowEdge(l3+1));
//    pH->Fit("rUg", "LRQ0");
//    cout << " l " << left  << ", " << lUg << "+-" << DlUg << " / " << fl->GetParameter(0) << "+-" << fl->GetParError(0) << endl
//         << " r " << right << ", " << rUg << "+-" << DrUg <<  " / " << fr->GetParameter(0) << "+-" << fr->GetParError(0) << endl;
    *lrUg = rUg - lUg;
    *DlrUg = sqrt( pow(DlUg, 2) + pow(DrUg, 2) );
}


void DrawLRUg(TFile* f, string FC)
{
    char name[64] = "";
    TMultiGraph* mg = new TMultiGraph("mgLR", "R-L Background Difference; Breite [ns]; R-L [1/s]");
    sprintf(name, "%s Deposit", FC.c_str());
    TLegend *l = new TLegend(0.9, 0.5, 1.0, 1.0, name);
    TGraphErrors* ge[8];
    Int_t W[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400}; // width samples
    Int_t n = sizeof(W) / sizeof(W[0]);
    Double_t X[n];
    Double_t Xerr[n];
    Double_t avLRUg[n];
    Double_t DavLRUg[n];
    for (Int_t j = 0; j < n; j++)
    {
        X[j] = W[j] / 4.0;
        Xerr[j] = 0;
        avLRUg[j] = 0;
        DavLRUg[j] = 0;
    }

    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i+1);
        TH1I* pH = (TH1I*)f->Get(name);
        Double_t LRUg[n];
        Double_t DLRUg[n];
        for (Int_t j = 0; j < n; j++)
        {
            GetLRUg(i, pH, FC, 60, W[j], &LRUg[j], &DLRUg[j]);
            avLRUg[j] += LRUg[j] / n;
            DavLRUg[j] = sqrt(pow(DavLRUg[j], 2) + pow(DLRUg[j] / n, 2));
//            cout << W[j] << ", " << NIF[j] << "+-" << DNIF[j] << endl;
        }
        ge[i] = new TGraphErrors(n, X, LRUg, Xerr, DLRUg);
        sprintf(name, "gLR_%i", i+1);
        ge[i]->SetName(name);
        ge[i]->SetLineColorAlpha(i+2, 0.5);
        ge[i]->SetLineWidth(1);
        ge[i]->SetMarkerColorAlpha(i+2, 0.5);
        ge[i]->SetMarkerStyle(20);
        ge[i]->SetMarkerSize(2);
        sprintf(name, "%i", i+1);
        l->AddEntry(ge[i], name);
        mg->Add(ge[i]);
    }
    TGraphErrors* av = new TGraphErrors(n, X, avLRUg, Xerr, DavLRUg);
    av->SetName("gLR_av");
    av->SetLineColorAlpha(1, 1.0);
    av->SetLineWidth(2);
    av->SetMarkerColorAlpha(1, 1.0);
    l->AddEntry(av, "av");
    mg->Add(av);

    TCanvas* c = new TCanvas("cLR", "Test Backgrounds", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    mg->Draw("AP");
    l->Draw();
    c->Modified();
    c->Update();
}


int DrawSingle(TFile* f, TFile* fCom, string FC, string Setup)
{ // Draw all pictures concerning one single file
    DrawLRUg(f, FC);
//    DrawQDCfit(f, FC, Setup);
//    DrawQDCeff(f, FC, Setup);
//    DrawQDCres(f, FC, Setup);
//    DrawDtInt(f, FC, Setup);
//    DrawDtUg(f, fCom, FC, Setup);
//    DrawDtRate(f, FC, Setup);
//    DrawDtPeak(f, fCom, FC, Setup);

    return 1;
}


int DrawPics()
{
//    gStyle->SetCanvasPreferGL();

    TFile* fNIF = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/NIF.root");
    TFile* fSB = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/SB.root");
//    TFile* fSF = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/SF.root");
    TFile* fUNIF = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/UFC_NIF.root");
    TFile* fUSB = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/UFC_SB.root");
    TFile* fCom = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/Evaluation.root");
    /// single pics
//    DrawSingle(fNIF, fCom, "PuFC", "NIF");
//    DrawSingle(fSB, fCom, "PuFC", "SB");
//    DrawSingle(fSF, fCom, "PuFC", "SF");
//    DrawSingle(fUNIF, fCom, "UFC", "NIF");
//    DrawSingle(fUSB, fCom, "UFC", "SB");

//    DrawScaledDt(fNIF, fSB, fSF);
//    DrawScaledDt(fUNIF, fUSB);
    /// common pics
//    DrawSFRate(fNIF, fSB, fSF, fCom);
//    DrawNPu(fCom);
//    DrawQDCmin(fNIF, fSB, fSF, fUNIF, fUSB, fCom);
//    DrawIntEff(fNIF, fSB, fSF, fCom);
//    DrawEff(fCom);
//    DrawInScat(fCom);
//    DrawEval(fCom);
//    DrawDtComp(fNIF, fUNIF, "NIF", 4);
//    DrawDtComp(fSB, fUSB, "SB", 4);
//    DrawDtComp(fNIF, fSB, fUNIF, fUSB, 4);
//    DrawDtFC(fNIF, fSB, "PuFC", 10);
//    DrawDtFC(fUNIF, fUSB, "UFC", 10);

    DrawStability();
//    for (Int_t i = 0; i < 8; i++)
//        DrawStability(i);
//    DrawPeak(fNIF, "PuFC");
//    DrawPeak(fUNIF, "UFC");

    return 1;
}
