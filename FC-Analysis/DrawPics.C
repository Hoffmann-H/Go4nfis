#define lw 2
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
    char hname[67] = "";
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
        fMax->SetLineWidth(lw);
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
        pH2->SetFillColor(kBlue);
//        pH2->SetFillColorAlpha(kBlue, 0.25); // possible if flagOpenGL.CanvasPreferGL is set to 1 in $ROOTSYS/etc/system.rootrc
        pH2->SetFillStyle(3001);
        pH2->Draw("same");
        fUg->SetLineColor(kRed);
        fUg->SetLineWidth(lw);
        fUg->SetFillColor(kRed);
        fUg->SetFillStyle(3001);
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
    char name[32] = "";
    sprintf(name, "%s_%s_rSF", FC.c_str(), Setup.c_str());
    /// Draw counting results
    TCanvas* c3 = new TCanvas(name, name, 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg3 = new TMultiGraph();
    mg3->SetTitle("SF detection rate; Deposit; events per second");

    TGraph* g6 = (TGraph*)f->Get("/Analysis/TimeDiff/SF_Rate");
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

        TGraph* g7 = (TGraph*)f->Get("/Analysis/TimeDiff/NIF_Rate");
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

        TGraph* g8 = (TGraph*)f->Get("/Analysis/TimeDiff/NIF_SF_Ratio");
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

        TGraph* g9 = (TGraph*)f->Get("/Analysis/TimeDiff/NIF_Flux_Ratio");
        g9->SetLineWidth(lw);
        g9->SetMarkerColor(kBlue);
        mg6->Add(g9);

        mg6->Draw("AP");
        c6->Modified();
        c6->Update();
    }
}

void DrawDtUg(TFile* f, string FC, string Setup)
{ // Draw TimeDiff integration
    char name[32] = "";
    TCanvas* pCAnaDt[8];
    TH1I* pHsum = new TH1I();
    char hname[67] = "";
    for (int i = 0; i < 8; i++)
    {
        // open Dt histogram nr i
        sprintf(hname, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i+1);
        TH1I* pH = (TH1I*)f->Get(hname);
        // open underground fit
        sprintf(hname, "/Analysis/TimeDiff/Ug/f%s%sUg_%i", FC.c_str(), Setup.c_str(), i+1);
        TGraph* fUg = (TGraph*)f->Get(hname);
        sprintf(hname, "/Analysis/TimeDiff/Ug/f%s%sUgPeak_%i", FC.c_str(), Setup.c_str(), i+1);
        TGraph* fUgP = (TGraph*)f->Get(hname);

        // create canvas
        sprintf(name, "%s_%s_DtUg_%i", FC.c_str(), Setup.c_str(), i+1);
        pCAnaDt[i] = new TCanvas(name, name, 200, 10, 700, 500);
        gPad->SetTicks(1, 1);

        // manipulate hists
        pHsum->Add(pH);
        pH->SetLineColor(kBlue);
        pH->SetAxisRange(62000, 80000, "X");
        pH->SetAxisRange(0, pH->GetMaximum(), "Y");
        TH1I* pH1 = (TH1I*)pH->Clone();
        TH1I* pH2 = (TH1I*)pH->Clone();
        pH2->SetFillColor(kBlue);
        pH2->SetFillStyle(3001);
        Double_t x0 = 0, x1 = 0, level = 0; // extract peak integration limits from underground graph
        fUgP->GetPoint(0, x0, level);
        fUgP->GetPoint(1, x1, level);
        cout << x0 << " " << x1 << " " << level << endl;
        for(Int_t bin = 0; bin < pH->GetNbinsX(); bin++)
        {
            if(pH->GetBinCenter(bin) > x0 && pH->GetBinCenter(bin) < x1)
            {
//                pH2->SetBinContent(bin, pH->GetBinContent(bin) - level);
                pH1->SetBinContent(bin, 0);
            }
            else
                pH2->SetBinContent(bin, 0);
        }

        // draw histogram and graphs together
        pH2->Draw();
        pH1->Draw("same");
        fUg->SetLineColor(kRed);
        fUg->SetLineWidth(lw);
        fUg->Draw("same");
        fUgP->SetLineColor(kOrange);
        fUgP->SetLineWidth(lw);
        fUgP->Draw("same");
        TGraphErrors* gSF = (TGraphErrors*)f->Get("/Analysis/TimeDiff/SF"); // get SF counts
        Double_t x = 0, SF = 0, NIF = 0;
        gSF->GetPoint(i, x, SF);
        char message[32] = "";
        sprintf(message, "SF events: %i", (int)(SF+0.5));
        TText* tSF = new TText();
        tSF->SetNDC();
        tSF->DrawText(0.2, 0.2, message);
        if(strcmp(Setup.c_str(), "SF")) // if beam was on
        {
            // draw number
            TGraphErrors* gNIF = (TGraphErrors*)f->Get("/Analysis/TimeDiff/NIF"); // get NIF counts
            gNIF->GetPoint(i, x, NIF);
            sprintf(message, "NIF events: %i", (int)(NIF+0.5));
            TText* tNIF = new TText();
            tNIF->SetNDC();
            tNIF->DrawText(0.6, 0.8, message);
        }
    }

    sprintf(name, "%s_%s_DtUg_all", FC.c_str(), Setup.c_str());
    TCanvas* pCsum = new TCanvas(name, name, 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg1 = new TMultiGraph();
    pHsum->SetNameTitle(name, "Time difference Sum over all deposits; #font[12]{t} / ch; counts");
    pHsum->Rebin(10);
    pHsum->SetAxisRange(62000, 80000, "X");
    pHsum->Draw();
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

    TGraphErrors* g7 = (TGraphErrors*)fCom->Get("/Analysis/Evaluation/eNIF");
    BiasX(g7, -0.1);
    g7->SetLineColor(kBlue);
    g7->SetLineWidth(lw);
    mg3->Add(g7);
    l3->AddEntry(g7, "Neutron-induced");

    TGraphErrors* g8 = (TGraphErrors*)fCom->Get("/Analysis/Evaluation/eSF");
    g8->SetLineColor(kGreen);
    g8->SetLineWidth(lw);
    mg3->Add(g8);
    l3->AddEntry(g8, "Spontaneaus (average)");

    TGraphErrors* g9 = (TGraphErrors*)fCom->Get("/Analysis/Evaluation/eRel");
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

int DrawSingle(TFile* f, string FC, string Setup)
{ // Draw all pictures concerning one single file
    DrawQDCfit(f, FC, Setup);
    DrawQDCeff(f, FC, Setup);
//    DrawQDCres(f, FC, Setup);
//    DrawDtInt(f, FC, Setup);
//    DrawDtUg(f, FC, Setup);

    return 1;
}


int DrawPics()
{
    TFile* fNIF = TFile::Open("/home/hoffma93/Go4nfis/offline/results/NIF.root");
    TFile* fSB = TFile::Open("/home/hoffma93/Go4nfis/offline/results/SB.root");
    TFile* fSF = TFile::Open("/home/hoffma93/Go4nfis/offline/results/SF.root");
    TFile* fUNIF = TFile::Open("/home/hoffma93/Go4nfis/offline/results/UFC_NIF.root");
    TFile* fUSB = TFile::Open("/home/hoffma93/Go4nfis/offline/results/UFC_SB.root");
    /// single pics
//    DrawSingle(fNIF, "PuFC", "NIF");
//    DrawSingle(fSB, "PuFC", "SB");
//    DrawSingle(fSF, "PuFC", "SF");
    DrawSingle(fUNIF, "UFC", "NIF");
    DrawSingle(fUSB, "UFC", "SB");

    /// common pics
    TFile* fCom = TFile::Open("/home/hoffma93/Go4nfis/offline/results/Evaluation.root");
    // SF
//    DrawSFRate(fNIF, fSB, fSF, fCom);
//    DrawNPu(fCom);
//    DrawQDCmin(fNIF, fSB, fSF, fUNIF, fUSB, fCom);
//    DrawIntEff(fNIF, fSB, fSF, fCom);
//    DrawEff(fCom);
//    DrawEval(fCom);
    return 1;
}
