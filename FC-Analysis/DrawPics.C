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

        sprintf(hname, "/Analysis/QDC/Fit/f%s%sMax_%i", FC.c_str(), Setup.c_str(), i+1);
        TF1* fMax = (TF1*)f->Get(hname);
        fMax->SetLineWidth(lw);//CanvasPreferGL
        fMax->SetLineColor(kGreen);
        legend->AddEntry(fMax, "pol4 Maximum fit");
        fMax->Draw("same");

        if(strcmp(FC.c_str(), "PuFC") == 0)
        { // if drawing PuFC results, draw Pedestal fit.
            sprintf(hname, "/Analysis/QDC/Fit/f%s%sPed_%i", FC.c_str(), Setup.c_str(), i+1);
            TF1* fPed = (TF1*)f->Get(hname);
            fPed->SetLineColor(kYellow);
            legend->AddEntry(fPed, "pol2 Pedestal fit");
            fPed->Draw("same");
        }
        legend->Draw("same");
    }
}

void DrawQDCeff(TFile* f, string FC, string Setup)
{ // Plot internal efficiency
    TCanvas* pC[8];
    char hname[67] = "";
    char name[67] = "";
    for (int i = 0; i < 8; i++)
    {
        // open Dt histogram nr i
        sprintf(hname, "/Histograms/Raw/QDC/low/H1RawQDCl_%i", i + 1);
        TH1I *pH = (TH1I*)f->Get(hname);
        pH->SetLineColor(kBlue);
//        legend->AddEntry(pHRawQDC, "Data");
        // open underground line
        sprintf(hname, "/Analysis/QDC/Ug/f%s%sUg_%i", FC.c_str(), Setup.c_str(), i+1);
        TGraph* fUg = (TGraph*)f->Get(hname);
        Double_t unused = 0, cut = 0;
        fUg->GetPoint(3, cut, unused);
        TH1I *pH1 = (TH1I*)pH->Clone();
        TH1I *pH2 = (TH1I*)pH->Clone();
        for(int bin = 0; bin < pH->GetNbinsX(); bin++)
        {
            if(pH->GetBinCenter(bin) < cut)
                pH2->SetBinContent(bin, 0);
            else
                pH1->SetBinContent(bin, 0);
        }
        // draw histogram and graph together
        sprintf(name, "%s_%s_QDCeff_%i", FC.c_str(), Setup.c_str(), i+1);
        pC[i] = new TCanvas(name, name, 200, 10, 700, 500);
        gPad->SetTicks(1, 1);
        gPad->SetLogy(1);
//        pH1->SetAxisRange(0, 7096, "X");
        pH1->Draw();
        pH2->SetFillColorAlpha(kBlue, 0.25);
//        pH2->SetFillStyle(3001);
        pH2->Draw("same");
        fUg->SetLineColor(kRed);
        fUg->SetLineWidth(lw);
        fUg->SetFillColorAlpha(kRed, 0.25);
//        fUg->SetFillStyle(3001);
//        legend->AddEntry(fUg, "FF constant extrapolation");
        fUg->Draw("fsame");
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


void DrawUGDt(TFile* f, TFile* fCom, string FC, string Setup)
{
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
        sprintf(hname, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i+1);
        TH1I* pH = (TH1I*)f->Get(hname);

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
        TH1F* pH1 = CopyRange(pH, name, 62000, 80000, - level);

        // create peak region TH1F histogram
        sprintf(name, "%s_%s_DtInt_%i", FC.c_str(), Setup.c_str(), i+1);
        TH1F* pH2 = CopyRange(pH, name, x0, x1, - level);

        // manipulate histograms
        TH1I* pH3 = (TH1I*)pH->Clone();
        BiasX(pH3, 0.5*(x0 + x1));
        pHsum->Add(pH3);
        pH2->SetLineColor(kRed);
        pH2->SetFillColorAlpha(kRed, 0.25);
        pH2->SetLineColor(kWhite);

        // draw histograms together
        sprintf(name, "%s_%s_DtPeak_%i", FC.c_str(), Setup.c_str(), i+1);
        pCDtPeak[i] = new TCanvas(name, name, 200, 10, 700, 500);
        gPad->SetTicks(1, 1);
        pH1->Draw();
        pH2->Draw("same");
        pH1->Draw("same");

        // draw number
        Double_t x = 0, SF = 0, NIF = 0; // prepare number drawing
        sprintf(name, "/%s/TimeDiff/Signal/%s/nNIF", FC.c_str(), Setup.c_str());
        TGraphErrors* gNIF = (TGraphErrors*)fCom->Get(name); // get NIF counts
        gNIF->GetPoint(i, x, NIF);
        sprintf(message, "NIF events: %i", (int)(NIF+0.5));
        TText* tNIF = new TText();
        tNIF->SetNDC();
        tNIF->DrawText(0.6, 0.8, message);

        //// Draw original Dt spectrum with underground lines
        sprintf(name, "%s_%s_DtUg_%i", FC.c_str(), Setup.c_str(), i+1);
        pCDtUg[i] = new TCanvas(name, name, 200, 10, 700, 500);
        gPad->SetTicks(1,1);
        pH->GetXaxis()->SetRangeUser(62000, 80000);
        pH->Draw();

        fUg->GetPoint(0, UgX[0], UgY[0]);
        fUgP->GetPoint(0, UgX[1], UgY[1]);
        TGraph* g0 = new TGraph(2, UgX, UgY);
        g0->SetLineColor(kRed);
        g0->SetLineWidth(lw);
        g0->Draw("same");
        fUgP->GetPoint(1, UgX[0], UgY[0]);
        fUg->GetPoint(1, UgX[1], UgY[1]);
        TGraph* g1 = new TGraph(2, UgX, UgY);
        g1->SetLineColor(kRed);
        g1->SetLineWidth(lw);
        g1->Draw("same");

        SF = pH->Integral() - NIF;
        sprintf(message, "Underground events: %i", (int)(SF+0.5));
        TText* tSF = new TText();
        tSF->SetNDC();
        tSF->DrawText(0.2, 0.2, message);
    }

    sprintf(name, "%s_%s_DtUg_all", FC.c_str(), Setup.c_str());
    TCanvas* pCsum = new TCanvas(name, name, 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    pHsum->SetNameTitle(name, "Time difference Sum over all deposits; #font[12]{t} / ch; counts");
    pHsum->Rebin(20);
    pHsum->SetAxisRange(-5000, 13000, "X");
    pHsum->Draw();
}


Double_t func_peak(Double_t *x, Double_t *par)
{
    return par[0] * exp( - pow((x[0] - par[1]) / par[2], 2)) + par[3];
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
        pH->SetAxisRange(62000, 80000, "X");
        fP->SetLineColor(kRed);
        fP->SetLineWidth(lw);
        // Draw
        pH->Draw();
        fP->Draw("same");
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


void DrawEff(TFile* fCom)
{ // efficiency conclusion
    TCanvas* c2 = new TCanvas("EffPuFC", "Efficiencies of PuFC deposits", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg2 = new TMultiGraph();
    mg2->SetTitle("Efficiencies of PuFC Deposits; Deposit; Ratio");
    TLegend* l2 = new TLegend(0.6, 0.2, 0.85, 0.70, "Source");

    TGraphErrors* g4 = (TGraphErrors*)fCom->Get("/Analysis/PuFC/Efficiency/eSF");
    BiasX(g4, -0.2);
    g4->SetLineWidth(lw);
    g4->SetLineColor(kRed);
    mg2->Add(g4);
    l2->AddEntry(g4, "Spontaneaus fission", "lp");

    TGraphErrors* g5 = (TGraphErrors*)fCom->Get("/Analysis/PuFC/Efficiency/eRel");
    BiasX(g5, -0.1);
    g5->SetLineWidth(lw);
    g5->SetLineColor(kGreen);
    mg2->Add(g5);
    l2->AddEntry(g5, "NIF:SF", "lp");

    TGraphErrors* g6 = (TGraphErrors*)fCom->Get("/Analysis/PuFC/Efficiency/eNIF");
    BiasX(g6, 0.0);
    g6->SetLineWidth(lw);
    g6->SetLineColor(kBlue);
    mg2->Add(g6);
    l2->AddEntry(g6, "Neutron-induced fission", "lp");

    TGraphErrors* g7 = (TGraphErrors*)fCom->Get("/Analysis/PuFC/Efficiency/eSimG");
    BiasX(g7, +0.1);
    g7->SetLineWidth(lw);
    g7->SetLineColor(kMagenta);
    mg2->Add(g7);
    l2->AddEntry(g7, "Simulation (Gayther)", "lp");

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


int DrawSingle(TFile* f, TFile* fCom, string FC, string Setup)
{ // Draw all pictures concerning one single file
//    DrawQDCfit(f, FC, Setup);
//    DrawQDCeff(f, FC, Setup);
//    DrawQDCres(f, FC, Setup);
//    DrawDtInt(f, FC, Setup);
//    DrawDtUg(f, fCom, FC, Setup);
    DrawDtPeak(f, fCom, FC, Setup);

    return 1;
}


int DrawPics()
{

//    gStyle->SetCanvasPreferGL();

    TFile* fNIF = TFile::Open("/home/hoffma93/Go4nfis/offline/results/NIF.root");
    TFile* fSB = TFile::Open("/home/hoffma93/Go4nfis/offline/results/SB.root");
    TFile* fSF = TFile::Open("/home/hoffma93/Go4nfis/offline/results/SF.root");
    TFile* fUNIF = TFile::Open("/home/hoffma93/Go4nfis/offline/results/UFC_NIF.root");
    TFile* fUSB = TFile::Open("/home/hoffma93/Go4nfis/offline/results/UFC_SB.root");
    TFile* fCom = TFile::Open("/home/hoffma93/Go4nfis/offline/results/Evaluation.root");
    /// single pics
//    DrawSingle(fNIF, fCom, "PuFC", "NIF");
//    DrawSingle(fSB, fCom, "PuFC", "SB");
    DrawSingle(fSF, fCom, "PuFC", "SF");
//    DrawSingle(fUNIF, fCom, "UFC", "NIF");
//    DrawSingle(fUSB, fCom, "UFC", "SB");

//    DrawScaledDt(fNIF, fSB, fSF);
    /// common pics
    // SF
//    DrawSFRate(fNIF, fSB, fSF, fCom);
//    DrawNPu(fCom);
//    DrawQDCmin(fNIF, fSB, fSF, fUNIF, fUSB, fCom);
//    DrawIntEff(fNIF, fSB, fSF, fCom);
//    DrawEff(fCom);
//    DrawEval(fCom);
    return 1;
}
