#include "Plot.h"
#define NumCh 8
#define lw 2

using namespace std;

Plot::Plot()
{
    for (Int_t i = 0; i < NumCh; i++)
    { // ch == {1,2,3,4,5,6,7,8}
        x[i] = i+1;
        xerr[i] = 0;
    }
}


Plot::Plot(string fc, string setup)
{
    FC = fc;
    Setup = setup;
    for (Int_t i = 0; i < NumCh; i++)
    { // ch == {1,2,3,4,5,6,7,8}
        x[i] = i+1;
        xerr[i] = 0;
    }
}


Plot::~Plot()
{

}


void Plot::simply(TH1F *pH, string Name, string xTitle, string yTitle)
{
    TCanvas *C = new TCanvas(Name.c_str(), Name.c_str(), 200, 10, 700, 500);
    gPad->SetTicks(1,1);
    pH->SetStats(0);
//    char name[64] = "";
//    sprintf(name, "%s; %s; %s", Name.c_str(), xTitle.c_str(), yTitle.c_str());
//    pH->SetTitle(name);
    pH->SetStats(0);
    pH->GetXaxis()->SetTitle(xTitle.c_str());
    pH->GetYaxis()->SetTitle(yTitle.c_str());
    pH->GetXaxis()->SetTitleSize(0.05);
    pH->GetYaxis()->SetTitleSize(0.05);
    pH->GetXaxis()->SetLabelSize(0.05);
    pH->GetYaxis()->SetLabelSize(0.05);
    pH->Draw("hist");
    C->Modified();
    C->Update();
}


void Plot::Sigma(TGraph *gPu242, TGraph *gU235, TGraph *gU238)
{
    TCanvas *pC1 = new TCanvas("Pu242", "Pu-242 cross sections", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    gPu242->SetTitle("Pu-242 (n,f); E / MeV; cross section / barn");
    gPu242->Draw();

    TCanvas *pC2 = new TCanvas("U235", "U-235 cross sections", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    gU235->SetTitle("U-235 (n,f); E / MeV; cross section / barn");
    gU235->Draw();

    TCanvas *pC3 = new TCanvas("U238", "U-238 cross sections", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    gU238->SetTitle("U-238 (n,f); E / MeV; cross section / barn");
    gU238->Draw();
}


void Plot::Source_E(TGraph *g, Double_t E, Double_t DE)
{
    char name[64] = "";
//    TGraph *g = new TGraph("/gpfs/home/hoffma93/Programme/ROOT/Data/Source_E.dat");
    TCanvas *c = new TCanvas("Source_E", "Source_E", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
//    sprintf(name, "E = %.2f +- %.2f MeV", E, DE);
    TLegend *l = new TLegend(0.3, 0.15, 0.7, 0.3, "E = <E> +- dE");

    g->SetNameTitle("Source_E", "Neutron Source Energy Spectrum");
    g->GetXaxis()->SetTitle("E / MeV");
    g->GetXaxis()->SetTitleSize(0.05);
    g->GetXaxis()->SetLabelSize(0.05);
    g->GetYaxis()->SetTitle("n(E)");
    g->GetYaxis()->SetTitleSize(0.05);
    g->GetYaxis()->SetLabelSize(0.05);
    g->SetLineWidth(lw);

    Double_t ymax = 0.816189;
    TLine *lLeft = new TLine(E-DE, 0, E-DE, ymax);
    lLeft->SetLineWidth(lw);
    lLeft->SetLineColor(kRed);
    TLine *lCenter = new TLine(E, 0, E, ymax);
    lCenter->SetLineWidth(lw);
    lCenter->SetLineColor(kBlue);
    sprintf(name, "<E> = %.2f MeV", E);
    l->AddEntry(lCenter, name);
    TLine *lRight = new TLine(E+DE, 0, E+DE, ymax);
    lRight->SetLineWidth(lw);
    lRight->SetLineColor(kRed);
    sprintf(name, "dE = %.2f MeV", DE);
    l->AddEntry(lRight, name);

    g->Draw();
    lLeft->Draw();
    lCenter->Draw();
    lRight->Draw();
    l->Draw();
    c->Modified();
    c->Update();
}


//void Plot::QDCfit(Int_t ch, TH1I *pH, Double_t cut, TF1 *fCut)
//{
//    char name[32] = "";
//    sprintf(name, "QDCfit_%i", ch+1);
//    TCanvas *pC = new TCanvas(name, name, 200, 10, 700, 500);
//    gPad->SetTicks(1, 1);
//    gPad->SetLogy(1);
//    TLegend* legend = new TLegend(0.6, 0.2, 0.85, 0.70, "Legend");

//    pH->SetStats(0);
//    pH->GetXaxis()->SetTitleSize(0.05);
//    pH->GetYaxis()->SetTitleSize(0.05);
//    pH->SetLineColor(kBlue);
//    legend->AddEntry(pH, "Data");
//    pH->Draw();

//    fCut->SetLineWidth(lw);
//    fCut->SetLineColor(kRed);
//    legend->AddEntry(fCut, "pol4 Minimum fit");
//    fCut->Draw("same"); // draw together

//    if (fMax != 0)
//    {
//        fMax->SetLineWidth(lw);//CanvasPreferGL
//        fMax->SetLineColor(kGreen);
//        legend->AddEntry(fMax, "pol4 Maximum fit");
//        fMax->Draw("same");
//    }
//        if(strcmp(FC.c_str(), "PuFC") == 0)
//        { // if drawing PuFC results, draw Pedestal fit.
//            sprintf(hname, "/Analysis/QDC/Fit/f%s%sPed_%i", FC.c_str(), Setup.c_str(), i+1);
//            TF1* fPed = (TF1*)f->Get(hname);
//            fPed->SetLineColor(kYellow);
//            legend->AddEntry(fPed, "pol2 Pedestal fit");
//            fPed->Draw("same");
//        }
//    legend->Draw("same");
//    pC->Modified();
//    pC->Update();
//}


void Plot::QDCeff(Int_t ch, TH1I *pH, Double_t pedestal, Double_t cut, Double_t level, Double_t eInt, Double_t DeInt)
{
    char name[64] = "";
    pH->SetLineColor(kBlue);
    pH->SetStats(0);
    pH->GetXaxis()->SetTitleSize(0.05);
    pH->GetYaxis()->SetTitleSize(0.05);
    sprintf(name, "Fission fragment energy deposition, channel %i", ch+1);
    pH->SetTitle(name);

    TH1I *pH1 = (TH1I*)pH->Clone();
//    pH1->SetAxisRange(0, 4096, "X");
    sprintf(name, "QDC signal integration, channel %i", ch+1);
    TH1F* pH2 = (TH1F*)CopyRange(pH, name, pH->GetBinLowEdge(cut), 4096, 0);
    pH2->SetLineWidth(0);
    pH2->SetFillColorAlpha(kBlue, 0.25);

    // make underground line
    Double_t xval[] = {pedestal, pedestal, pH->GetBinLowEdge(cut), pH->GetBinLowEdge(cut)};
    Double_t yval[] = {0, level, level, 0};
    TGraph* fUg = new TGraph(4, xval, yval);
    sprintf(name, "fUg_%i", ch+1);
    fUg->SetNameTitle(name, "FF underground");
    fUg->SetLineWidth(lw);
    fUg->SetLineColor(kRed);
    fUg->SetFillColorAlpha(kRed, 0.25);

    // Make legend
    TLegend *l = new TLegend(0.6, 0.85, 0.95, 0.95, "");
    sprintf(name, "Channel %i", ch+1);
    l->AddEntry(pH1, name);
//        legend->AddEntry(fUg, "FF constant extrapolation");

    // draw histogram and graph together
    sprintf(name, "QDCeff_%i", ch+1);
    TCanvas *pC = new TCanvas(name, name, 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    gPad->SetLogy(1);

    pH1->Draw("hist");
    pH2->Draw("same");
    fUg->Draw("fsame");
    l->Draw();
    // Show eff number
    sprintf(name, "eff =  %.3f", eInt);
    TText* tNIF = new TText();
    tNIF->SetNDC();
    tNIF->DrawText(0.2, 0.2, name);
//    TText* thr = new TText();
//    thr->SetNDC();
//    thr->DrawText(0.25, 0.5, "thr");
//    TText* FF = new TText();
//    FF->SetNDC();
//    FF->DrawText(0.5, 0.6, "Fission fragments");
//    TText* alpha = new TText();
//    alpha->SetNDC();
//    alpha->DrawText(0.15, 0.6, "Alphas");

    pC->Modified();
    pC->Update();
    cout << pC->GetName() << "  " << pH1->GetName() << "  " << pH2->GetName() << "  " << fUg->GetName() << endl;
}


void Plot::ExpT(Double_t* uT, Double_t* DuT, Double_t* T, Double_t* DT)
{
    cout << endl << "Plotting Exp. T" << endl;
//    for (int i = 0; i < NumCh; i++)
//        cout << uT[i] << " " << DuT[i] << " " << T[i] << " " << DT[i] << endl;
    Double_t X[] = {1,2,3,4,5,6,7,8};
    Double_t Xerr[] = {0,0,0,0,0,0,0,0};
    TCanvas * c1 = new TCanvas("uExpT", "(n,f) rate over SF rate", 200, 10, 700, 500);
    TGraphErrors *g1 = new TGraphErrors(NumCh, X, uT, Xerr, DuT);
    g1->SetNameTitle("uncExpT", "(n,f) rate over SF rate; Channel; N(n,f) / N(SF)");
    gPad->SetTicks(1, 1);
    g1->Draw();
    c1->Modified();
    c1->Update();

    TCanvas * c2 = new TCanvas("ExpT", "Experimental tramsmission", 200, 10, 700, 500);
    TGraphErrors *g2 = new TGraphErrors(NumCh, X, T, Xerr, DT);
    gPad->SetTicks(1, 1);
    g2->Draw();
    c2->Modified();
    c2->Update();
}


void Plot::SimT(Double_t* T, Double_t* DT, Double_t* S, Double_t* DS, Double_t* F, Double_t* DF)
{
    cout << endl << "Plotting simulated correction factors" << endl;
    for (int i = 0; i < NumCh; i++)
        cout << F[i] << " " << DF[i] << endl;
    Double_t X1[] = {0.9,1.9,2.9,3.9,4.9,5.9,6.9,7.9};
    Double_t X2[] = {1,2,3,4,5,6,7,8};
    Double_t X3[] = {1.1,2.1,3.1,4.1,5.1,6.1,7.1,8.1};
    Double_t Xerr[] = {0,0,0,0,0,0,0,0};

    TCanvas * c1 = new TCanvas("SimFac", "Simulation correction factors", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg = new TMultiGraph();
    mg->SetTitle("Correction factor simulation; Channel; Factor");
    TLegend* l = new TLegend(0.15, 0.4, 0.4, 0.6, "Legend");

    TGraphErrors *g1 = new TGraphErrors(NumCh, X1, T, Xerr, DT);
    g1->SetName("SimT");
    g1->SetLineColor(kRed);
    g1->SetLineWidth(lw);
    mg->Add(g1);
    l->AddEntry(g1, "Transmission", "lp");

    TGraphErrors *g2 = new TGraphErrors(NumCh, X2, S, Xerr, DS);
    g2->SetName("SimS");
    g2->SetLineColor(kGreen);
    g2->SetLineWidth(lw);
    mg->Add(g2);
    l->AddEntry(g2, "Scattering", "lp");

    TGraphErrors *g3 = new TGraphErrors(NumCh, X3, F, Xerr, DF);
    g3->SetName("SimF");
    g3->SetLineWidth(lw);
    mg->Add(g3);
    l->AddEntry(g3, "Correction factor", "lp");

    mg->Draw("AP");
    l->Draw();
    c1->Modified();
    c1->Update();
    cout << "Done: Draw correction factor simulation" << endl;
}


void Plot::CompSc(Double_t* pExp, Double_t* DpExp, Double_t* pSim1, Double_t* DpSim1, Double_t* pSim2, Double_t* DpSim2)
{
    cout << endl << "Plotting scattering comparison..." << endl;
    Double_t X[] = {1,2,3,4,5,6,7,8};
    Double_t Xerr[] = {0,0,0,0,0,0,0,0};
    char name[64] = "";

    TCanvas *c1 = new TCanvas("CompSc", "Compare scattering");
    gPad->SetTicks(1, 1);
    TMultiGraph* mg = new TMultiGraph();
    sprintf(name, "%s Scattering; Channel; unsc. portion", FC.c_str());
    mg->SetTitle(name);
    TLegend *l = new TLegend(0.2, 0.4, 0.5, 0.6, "Legend");

    TGraphErrors *g1 = new TGraphErrors(NumCh, X, pExp, Xerr, DpExp);
    g1->SetName("pExp");
    g1->SetLineColor(kRed);
    g1->SetLineWidth(lw);
    mg->Add(g1);
    l->AddEntry(g1, "FG+BG Measurement");

    TGraphErrors *g2 = new TGraphErrors(NumCh, X, pSim1, Xerr, DpSim1);
    g2->SetName("pSim1");
    g2->SetLineColor(kGreen);
    g2->SetLineWidth(lw);
    mg->Add(g2);
    l->AddEntry(g2, "FG Simulation");

    TGraphErrors *g3 = new TGraphErrors(NumCh, X, pSim2, Xerr, DpSim2);
    g3->SetName("pSim2");
    g3->SetLineColor(kBlue);
    g3->SetLineWidth(lw);
    mg->Add(g3);
    l->AddEntry(g3, "FG+BG Simulation");

    mg->Draw("AP");
    l->Draw();
    c1->Modified();
    c1->Update();
    cout << "Done: Scattering comparison" << endl;

}


void Plot::Dt(Int_t ch, TH1I *pH, Double_t nf, Double_t Dnf, Int_t l0, Int_t l1, Int_t l2, Int_t l3, Double_t level, string Setup)
{
//    cout << "Plotting Dt..." << endl << "(n,f) " << nf << endl;
    char name[64] = "";
    sprintf(name, "Dt_%i", ch+1);
    TCanvas *C = new TCanvas(name, "Time differecne spectrum", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    sprintf(name, "Ug_%s_%i", Setup.c_str(), ch+1);
    TH1F *pHbg = CopyRange(pH, name, 62000, 80000, -level);
    sprintf(name, "Peak_%s_%i", Setup.c_str(), ch+1);
    TH1F *pHpeak = CopyBins(pH, name, l1, l2, -level);
    Double_t BgX[] = {pH->GetBinLowEdge(l0), pH->GetBinLowEdge(l3+1)};
    Double_t BgY[] = {0, 0};
    TGraph *g = new TGraph(2, BgX, BgY);
    g->SetLineColor(kRed);
    g->SetLineWidth(lw);
    pHbg->SetStats(0);
    sprintf(name, "Bg-subtracted time spectrum, ch.%i; TDC channel; Fissions", ch+1);
    pHbg->SetTitle(name);
    pHbg->Draw();
    g->Draw("same");
    pHpeak->SetLineWidth(0);
    pHpeak->SetFillColorAlpha(kBlue, 0.25);
    pHpeak->Draw("same");
//    g->Draw("same");
//    pH->Draw("same");
    sprintf(name, "(n,f) events: %.0f +- %.0f", nf, Dnf);
    TText *tnf = new TText();
    tnf->SetNDC();
    tnf->DrawText(0.5, 0.8, name);
    sprintf(name, "SF background: %.2f events / bin", level);
    TText *tbg = new TText();
    tbg->SetNDC();
    tbg->DrawText(0.2, 0.2, name);
    C->Modified();
    C->Update();
}


void Plot::TvsE(Int_t ch, TH2F *pH2TvsE, Double_t binToFmin, Double_t binToFmax, Int_t directN)
{
//    cout << "Going to plot 2D E-t-hist for channel " << ch+1 << endl;
    Double_t ToFmin = pH2TvsE->GetXaxis()->GetBinLowEdge(binToFmin);
    Double_t ToFmax = pH2TvsE->GetXaxis()->GetBinLowEdge(binToFmax + 1);
    Int_t NbinsE = pH2TvsE->GetNbinsY();
    Int_t ScatteredN = pH2TvsE->Integral(binToFmin, binToFmax, 0, NbinsE+1) - directN;

    char name[64] = "";
    sprintf(name, "%s_%s_ToFvsEkin_Ch.%i", FC.c_str(), Setup.c_str(), ch+1);
    TCanvas *C = new TCanvas(name, name, 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    gPad->SetLogz(1);
//    gStyle->SetPalette(kRainBow);
    pH2TvsE->Rebin2D(1, 8);
    sprintf(name, "%s, %s, ToF vs Ekin, Ch.%i", FC.c_str(), Setup.c_str(), ch+1);
    pH2TvsE->SetTitle(name);
    pH2TvsE->SetStats(0);
    Double_t xmax = pH2TvsE->GetXaxis()->GetBinLowEdge(pH2TvsE->GetNbinsX()+1);
    Double_t ymax = pH2TvsE->GetYaxis()->GetBinLowEdge(pH2TvsE->GetNbinsY()+1);
    TLine *lLeft = new TLine(ToFmin, 0, ToFmin, ymax);
    lLeft->SetLineColor(kRed);
    TLine *lRight = new TLine(ToFmax, 0, ToFmax, ymax);
    lRight->SetLineColor(kRed);

    pH2TvsE->Draw("colz");
    lLeft->Draw();
    lRight->Draw();

    TText *tD = new TText();
    tD->SetNDC();
    sprintf(name, "%i direct neutrons", directN);
    tD->DrawText(0.3, 0.8, name);
    TText *tS = new TText();
    tS->SetNDC();
    sprintf(name, "%i scattered neutrons", ScatteredN);
    tS->DrawText(0.3, 0.5, name);

    C->Modified();
    C->Update();
}


void Plot::TvsE(Int_t ch, TH2F *pH2TvsE, Double_t binToFmin, Double_t binToFmax, Double_t binEmin, Double_t binEmax)
{
//    cout << "Going to plot 2D E-t-hist for channel " << ch+1 << endl;
    Double_t ToFmin = pH2TvsE->GetXaxis()->GetBinLowEdge(binToFmin);
    Double_t ToFmax = pH2TvsE->GetXaxis()->GetBinLowEdge(binToFmax + 1);
    Double_t Emin = pH2TvsE->GetYaxis()->GetBinLowEdge(binEmin);
    Double_t Emax = pH2TvsE->GetYaxis()->GetBinLowEdge(binEmax + 1);
//    cout << ToFmin << " " << ToFmax << " " << Emin << " " << Emax << endl;
    Int_t DirectN = pH2TvsE->Integral(binToFmin, binToFmax, binEmin, binEmax);
    Int_t ScatteredN = pH2TvsE->Integral(binToFmin, binToFmax, 0, binEmin-1);
    cout << DirectN << ";  " << ScatteredN << endl;

    char name[64] = "";
    sprintf(name, "ToFvsEkin_Ch.%i", ch+1);
    TCanvas *C = new TCanvas(name, name, 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    gPad->SetLogz(1);
//    gStyle->SetPalette(kRainBow);
    pH2TvsE->Rebin2D(1, 8);
    pH2TvsE->SetTitle("");
    pH2TvsE->SetStats(0);
    Double_t xmax = pH2TvsE->GetXaxis()->GetBinLowEdge(pH2TvsE->GetNbinsX()+1);
    Double_t ymax = pH2TvsE->GetYaxis()->GetBinLowEdge(pH2TvsE->GetNbinsY()+1);
    TLine *lLeft = new TLine(ToFmin, 0, ToFmin, ymax);
    lLeft->SetLineColor(kRed);
    TLine *lRight = new TLine(ToFmax, 0, ToFmax, ymax);
    lRight->SetLineColor(kRed);
    TLine *lBottom = new TLine(0, Emin, xmax, Emin);
    lBottom->SetLineColor(kRed);
    TLine *lTop = new TLine(0, Emax, xmax, Emax);
    lTop->SetLineColor(kRed);

    pH2TvsE->Draw("colz");
    lLeft->Draw();
    lRight->Draw();
    lBottom->Draw();
    lTop->Draw();

    TText *tD = new TText();
    tD->SetNDC();
    sprintf(name, "%i direct neutrons", DirectN);
    tD->DrawText(0.3, 0.8, name);
    TText *tS = new TText();
    tS->SetNDC();
    sprintf(name, "%i scattered neutrons", ScatteredN);
    tS->DrawText(0.3, 0.5, name);

    C->Modified();
    C->Update();
}


void Plot::DtPeakForm(Int_t ch, TH1F *pSc, TH1F *pUnsc, Double_t Emin, Int_t Projectiles)
{
//    cout << nSc << "  " << nUnsc << "  " << Projectiles << endl;
    char name[64] = "";
    sprintf(name, "Fission time spectrum channel %i", ch+1);
    TCanvas *C = new TCanvas(name, name, 200, 10, 700, 500);
    gPad->SetTicks(1,1);
    gPad->SetLogy(1);
//    pUnsc->Add(pSc, -1);

    pUnsc->SetStats(0);
    pUnsc->GetXaxis()->SetTitle("ToF / ns");
    pUnsc->GetYaxis()->SetTitle("N");
    pUnsc->GetXaxis()->SetTitleSize(0.05);
    pUnsc->GetYaxis()->SetTitleSize(0.05);
    pUnsc->GetXaxis()->SetLabelSize(0.05);
    pUnsc->GetYaxis()->SetLabelSize(0.05);
    pUnsc->SetLineColor(kBlue);
    pUnsc->SetLineWidth(lw);
    pSc->SetLineColor(kRed);
    pSc->SetLineWidth(lw);

    sprintf(name, "%i projectiles", Projectiles);
    TLegend *l = new TLegend(0.4, 0.6, 0.8, 0.8, name);
    sprintf(name, "E > %.2f MeV", Emin);
    l->AddEntry(pUnsc, "Unscattered");
    sprintf(name, "E < %.2f MeV", Emin);
    l->AddEntry(pSc, "Scattered");
    pUnsc->Draw();
    pSc->SetLineColor(kRed);
    pSc->Draw("same");
    l->Draw();
    C->Update();
}


TH1F* Plot::CopyRange(TH1I* pH, char* name, Double_t x0, Double_t x1, Double_t yoffset)
{
    Double_t offset = pH->GetBinLowEdge(0);
    Double_t ChPerBin = pH->GetBinWidth(0);
    int b0 = (x0 - offset) / ChPerBin;
    int b1 = (x1 - offset) / ChPerBin + 1;
    return CopyBins(pH, name, b0, b1, yoffset);
}


TH1F* Plot::CopyBins(TH1I *pH, char *name, Int_t b0, Int_t b1, Double_t yoffset)
{
    TH1F* pH2 = new TH1F(name, name, b1 - b0, pH->GetBinLowEdge(b0), pH->GetBinLowEdge(b1));
    for (int bin = 0; bin < pH2->GetNbinsX(); bin++)
        pH2->SetBinContent(bin + 1, pH->GetBinContent(bin + b0) + yoffset);
    return pH2;
}
