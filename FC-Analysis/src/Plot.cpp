#include "Plot.h"
#include "TGaxis.h"
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


Plot::Plot(string fc, string setup, TFile *file)
{
    FC = fc;
    Setup = setup;
    f = file;
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
    cout << endl << "Cross sections input @ 15MeV" << endl
         << " 242Pu " << gPu242->Eval(15) << endl
         << " 235U  " << gU235->Eval(15) << endl
         << " 238U  " << gU238->Eval(15) << endl;
    TCanvas *pC1 = new TCanvas("Pu242", "Pu-242 cross sections", 200, 10, 700, 500);
    gPad->SetTicks(1, 0);
    gPu242->SetTitle("^{242}Pu (n,f); #font[12]{E} [MeV]; #sigma [b]");
    gPu242->GetXaxis()->SetLabelSize(0.05);
    gPu242->GetXaxis()->SetTitleSize(0.07);
    gPu242->GetXaxis()->SetTitleOffset(0.7);
    gPu242->GetXaxis()->SetLimits(0, 20);
    gPu242->GetYaxis()->SetLabelSize(0.05);
    gPu242->GetYaxis()->SetTitleSize(0.07);
    gPu242->GetYaxis()->SetTitleOffset(0.6);
    gPu242->GetXaxis()->SetRangeUser(0, 16);
    Double_t leftmax = 2.4;
    gPu242->GetYaxis()->SetLimits(0, leftmax);
    gPu242->SetLineWidth(lw);
    gPu242->SetLineColor(kBlue);
    gPu242->Draw();

    TGaxis *axis = new TGaxis(16, 0, 16, gPu242->Eval(15), 0, 1, 509,"+L");
//    axis->SetLabelOffset(0.5);
    axis->SetTitle("#frac{#sigma#font[42]{(}#font[12]{E}#font[42]{)}}{#sigma#font[42]{(15MeV)}}");
    axis->SetLabelFont(42);
    axis->SetLabelSize(0.05);
    axis->SetTitleSize(0.07);
    axis->SetTitleOffset(0.8);
    axis->Draw();

    TLine *lh = new TLine(0, gPu242->Eval(15), 16, gPu242->Eval(15));
    lh->Draw();
    TLine *lv = new TLine(15, 0, 15, leftmax);
    lv->Draw();

    TCanvas *pC2 = new TCanvas("U235", "U-235 cross sections", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    gU235->SetTitle("U-235 (n,f); #font[12]{E} / MeV; cross section / barn");
    gU235->Draw();

    TCanvas *pC3 = new TCanvas("U238", "U-238 cross sections", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    gU238->SetTitle("U-238 (n,f); #font[12]{E} / MeV; cross section / barn");
    gU238->Draw();
}


void Plot::Source_E(TGraph *g, Double_t E, Double_t DE, Int_t Variante)
{
    if (Variante == 0)
        Source_E_v0(g, E, DE);
    else
        Source_E_v1(g, E, DE);
}


void Plot::Source_E_v0(TGraph *g, Double_t E, Double_t DE)
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


void Plot::Source_E_v1(TGraph *g, Double_t E, Double_t DE)
{
    char name[64] = "";
    TCanvas *c = new TCanvas("Source_E", "Source_E", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);

    g->SetNameTitle("Source_E", "Neutron Source Energy Spectrum");
    g->GetXaxis()->SetTitle("#font[12]{E} / MeV");
    g->GetXaxis()->SetTitleSize(0.08);
    g->GetXaxis()->SetLabelSize(0.06);
    g->GetXaxis()->SetTitleOffset(0.7);
    g->GetYaxis()->SetTitle("#font[12]{n}(#font[12]{E}) [a.u.]");
    g->GetYaxis()->SetTitleSize(0.08);
    g->GetYaxis()->SetLabelSize(0.06);
    g->GetYaxis()->SetTitleOffset(0.6);
    g->SetLineWidth(lw);

    sprintf(name, "E = %.2f +- %.2f MeV", E, DE);
    TLegend *l = new TLegend(0.4, 0.2, 0.6, 0.3, name);

    g->Draw();
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

    sprintf(name, "%s_%s_QDCeff_%i", FC.c_str(), Setup.c_str(), ch+1);
    TCanvas *pC = new TCanvas(name, name, 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    gPad->SetLogy(1);

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
    Int_t ff = pH->Integral(cut, 4096);
    sprintf(name, "%i Fission fragments", ff);
    TText* FF = new TText();
    FF->SetNDC();
    FF->DrawText(0.5, 0.6, name);
//    TText* alpha = new TText();
//    alpha->SetNDC();
//    alpha->DrawText(0.15, 0.6, "Alphas");

    pC->Modified();
    pC->Update();
    cout << pC->GetName() << "  " << pH1->GetName() << "  " << pH2->GetName() << "  " << fUg->GetName() << endl;
}


void Plot::ExpT(Double_t* expT, Double_t* DexpT, Double_t* simT, Double_t* DsimT)
{
    cout << endl << "Plotting experimental transmission" << endl;
//    for (int i = 0; i < NumCh; i++)
//        cout << uT[i] << " " << DuT[i] << " " << T[i] << " " << DT[i] << endl;
    Double_t X1[] = {0.9, 1.9, 2.9, 3.9, 4.9, 5.9, 6.9, 7.9};
    Double_t X2[] = {1,2,3,4,5,6,7,8};
    Double_t X3[] = {1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 8.1};
    Double_t Xerr[] = {0,0,0,0,0,0,0,0};
    char name[64] = "";


//    for (Int_t i = 0; i < NumCh; i++)
//    {
//        cout << uT[i] << "+-" << DuT[i] << " " << T[i] << "+-" << DT[i] << endl;
//        ut[i] = uT[i]; Dut[i] = DuT[i]; t[i] = T[i]; Dt[i] = DT[i];
//        aT[i] = aT[i] / aT[7] * uT[7];
//        DaT[i] = DaT[i] / aT[7] * uT[7];
//    }

    sprintf(name, "c_%s_ExpT", FC.c_str());
    TCanvas * c1 = new TCanvas(name, "(n,f) rate over SF rate", 200, 10, 700, 500);
    gPad->SetTicks(1, 0);
    TMultiGraph* mg = new TMultiGraph();
    sprintf(name, "g_%s_ExpT", FC.c_str());
    if (strcmp(FC.c_str(), "UFC"))
        mg->SetNameTitle(name, "; Kanal; #font[12]{T} [%]");
    else
        mg->SetNameTitle(name, "n-induced fission count over efficient mass density; Channel; fission yield / a.u.");
    TLegend* l = new TLegend(0.15, 0.4, 0.4, 0.6);

    Double_t scale = (simT[0]+simT[1]+simT[2]+simT[3]+simT[4]+simT[5]+simT[6]+simT[7])
                   / (expT[0]+expT[1]+expT[2]+expT[3]+expT[4]+expT[5]+expT[6]+expT[7]);

    Double_t expTscaled[NumCh];
    Double_t DexpTscaled[NumCh];
    Double_t simTscaled[NumCh];
    Double_t DsimTscaled[NumCh];
    Double_t anaTscaled[NumCh];
    Double_t DanaTscaled[NumCh];
    for (int i = 0; i < NumCh; i++)
    {
        expTscaled[i] = 100 * scale * expT[i];
        DexpTscaled[i] = 100 * scale * DexpT[i];
        simTscaled[i] = 100 * simT[i];
        DsimTscaled[i] = 100 * DsimT[i];
    }
    TGraphErrors *g1 = new TGraphErrors(NumCh, X1, expTscaled, Xerr, DexpTscaled);
    g1->SetName("uExpT");
    g1->SetLineColor(kRed);
    g1->SetMarkerColor(kRed);
    g1->SetLineWidth(lw);
    g1->SetMarkerStyle(21);
    g1->SetMarkerSize(2);
    mg->Add(g1);
    l->AddEntry(g1, "Experiment, stat.", "lp");

    TGraphErrors *g2 = new TGraphErrors(NumCh, X2, simTscaled, Xerr, DsimTscaled);
    g2->SetName("ExpT");
//    g2->SetLineColor(kGreen);
//    g2->SetLineWidth(lw);
    g2->SetMarkerStyle(20);
    g2->SetMarkerSize(2);
    mg->Add(g2);
    l->AddEntry(g2, "Simulation, stat.", "lp");

    TGraphErrors *g3;
    if (strcmp(FC.c_str(), "UFC")) {
        Double_t anaT[] = {094.59, 095.24, 095.89, 096.55, 097.21, 097.87, 098.54, 099.22};
        Double_t DanaT[] = {0.12, 0.10, 0.09, 0.07, 0.05, 0.04, 0.02, 0.01};
        g3 = new TGraphErrors(NumCh, X3, anaT, Xerr, DanaT);
    } else {
        Double_t anaT[] = {0.9465, 0.9529, 0.9594, 0.9658, 0.9724, 0.9789, 0.9856, 0.9922};
        Double_t DanaT[] = {0.0012, 0.0010, 0.0009, 0.0007, 0.0005, 0.0004, 0.0002, 0.0001};
        g3 = new TGraphErrors(NumCh, X3, anaTscaled, Xerr, DanaTscaled);
    }
    g3->SetName("ExpT");
//    g3->SetLineColor(kBlue);
//    g3->SetLineWidth(lw);
    g3->SetMarkerStyle(22);
    g3->SetMarkerSize(2);
    mg->Add(g3);
    l->AddEntry(g3, "Analytisch, syst.", "lp");

    mg->Draw("AP");
//    mg->GetXaxis()->SetRangeUser(0, 9);
    mg->GetXaxis()->SetNdivisions(8);
    mg->GetXaxis()->SetLabelSize(0.05);
    mg->GetXaxis()->SetTitleSize(0.08);
    mg->GetXaxis()->SetTitleOffset(0.6);
//    mg->GetYaxis()->SetRangeUser(0, 1.2);
    mg->GetYaxis()->SetLabelSize(0.05);
    mg->GetYaxis()->SetTitleSize(0.08);
    mg->GetYaxis()->SetTitleOffset(0.6);
    TGaxis *ax = new TGaxis(mg->GetXaxis()->GetXmax(), mg->GetYaxis()->GetXmin(), mg->GetXaxis()->GetXmax(), mg->GetYaxis()->GetXmax(), mg->GetYaxis()->GetXmin() / scale * 10, mg->GetYaxis()->GetXmax() / scale * 10, 509, "+L");
    ax->SetLineColor(kRed);
    ax->SetLabelColor(kRed);
    ax->SetLabelSize(0.05);
    ax->SetTextColor(kRed);
    ax->SetTitle("#frac{#font[12]{C}_{(n,f)}}{#font[12]{C}_{SF}} [10^{-3}]");
    ax->SetTitleSize(0.08);
    ax->SetTitleOffset(0.7);
    ax->Draw();
    l->Draw();
    c1->Modified();
    c1->Update();
}


void Plot::SimF(Double_t* T, Double_t* DT, Double_t* S, Double_t* DS, Double_t* F, Double_t* DF)
{
    cout << endl << "Plotting simulated correction factors" << endl;
    for (int i = 0; i < NumCh; i++)
        cout << F[i] << " " << DF[i] << endl;
    Double_t X1[] = {0.9,1.9,2.9,3.9,4.9,5.9,6.9,7.9};
    Double_t X2[] = {1,2,3,4,5,6,7,8};
    Double_t X3[] = {1.1,2.1,3.1,4.1,5.1,6.1,7.1,8.1};
    Double_t Xerr[] = {0,0,0,0,0,0,0,0};
    char name[64] = "";

    sprintf(name, "%s_SimFac", FC.c_str());
    TCanvas * c1 = new TCanvas(name, "Simulation correction factors", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg = new TMultiGraph();
    sprintf(name, "%s Correction factor simulation; Channel; Factor", FC.c_str());
    mg->SetTitle(name);
    TLegend* l = new TLegend(0.15, 0.4, 0.4, 0.6, "Legend");

    TGraphErrors *g1 = new TGraphErrors(NumCh, X2, T, Xerr, DT);
    g1->SetName("SimT");
    g1->SetMarkerStyle(23);
    g1->SetMarkerColor(kRed);
    g1->SetLineColor(kRed);
    g1->SetLineWidth(lw);
    mg->Add(g1);
    l->AddEntry(g1, "Transmission", "p");

    TGraphErrors *g2 = new TGraphErrors(NumCh, X2, S, Xerr, DS);
    g2->SetName("SimS");
    g2->SetMarkerStyle(22);
    g2->SetMarkerColor(kBlue);
    g2->SetLineColor(kBlue);
    g2->SetLineWidth(lw);
    mg->Add(g2);
    l->AddEntry(g2, "Scattering", "p");

    TGraphErrors *g3 = new TGraphErrors(NumCh, X2, F, Xerr, DF);
    g3->SetName("SimF");
    g3->SetMarkerStyle(33);
    g3->SetMarkerSize(2);
    g3->SetLineWidth(lw);
    mg->Add(g3);
    l->AddEntry(g3, "Correction factor", "p");

    mg->Draw("AP");
    l->Draw();
    c1->Modified();
    c1->Update();
    cout << "Done: Draw correction factor simulation" << endl;
}


void Plot::SimF(Double_t *T, Double_t *DT, Double_t *S, Double_t *DS, Double_t *F, Double_t *DF, Double_t *T2, Double_t *DT2, Double_t *S2, Double_t *DS2, Double_t *F2, Double_t *DF2)
{
    cout << endl << "Plotting simulated correction factors for PuFC and UFC" << endl;

    Double_t X1[] = {0.95,1.95,2.95,3.95,4.95,5.95,6.95,7.95};
    Double_t X2[] = {1,2,3,4,5,6,7,8};
    Double_t X3[] = {1.05,2.05,3.05,4.05,5.05,6.05,7.05,8.05};
    Double_t Xerr[] = {0,0,0,0,0,0,0,0};
    char name[64] = "";

    TCanvas * c1 = new TCanvas("SimFac", "Simulation correction factors", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg = new TMultiGraph();
    mg->SetTitle("Correction factor simulation; Kanal; Faktor");
    TLegend* l = new TLegend(0.15, 0.4, 0.4, 0.6);

    TGraphErrors *g1 = new TGraphErrors(NumCh, X1, T, Xerr, DT);
    g1->SetName("PuFC_SimT");
    g1->SetMarkerStyle(23);
    g1->SetMarkerSize(2);
    g1->SetMarkerColor(kBlue);
    g1->SetLineColor(kBlue);
    g1->SetLineWidth(lw);
    mg->Add(g1);
    l->AddEntry(g1, "PuFC Transmission", "p");

    TGraphErrors *g2 = new TGraphErrors(NumCh, X1, S, Xerr, DS);
    g2->SetName("PuFC_SimS");
    g2->SetMarkerStyle(22);
    g2->SetMarkerSize(2);
    g2->SetMarkerColor(kBlue);
    g2->SetLineColor(kBlue);
    g2->SetLineWidth(lw);
    mg->Add(g2);
    l->AddEntry(g2, "PuFC Streuung", "p");

    TGraphErrors *g3 = new TGraphErrors(NumCh, X1, F, Xerr, DF);
    g3->SetName("PuFC_SimF");
    g3->SetMarkerStyle(21);
    g3->SetMarkerSize(2);
    g3->SetMarkerColor(kBlue);
    g3->SetLineColor(kBlue);
    g3->SetLineWidth(lw);
    mg->Add(g3);
    l->AddEntry(g3, "PuFC Korrektur", "p");

    TGraphErrors *g4 = new TGraphErrors(NumCh, X3, T2, Xerr, DT2);
    g4->SetName("PuFC_SimT");
    g4->SetMarkerStyle(23);
    g4->SetMarkerSize(2);
    g4->SetMarkerColor(kRed);
    g4->SetLineColor(kRed);
    g4->SetLineWidth(lw);
    mg->Add(g4);
    l->AddEntry(g4, "UFC Transmission", "p");

    TGraphErrors *g5 = new TGraphErrors(NumCh, X3, S2, Xerr, DS2);
    g5->SetName("PuFC_SimS");
    g5->SetMarkerStyle(22);
    g5->SetMarkerSize(2);
    g5->SetMarkerColor(kRed);
    g5->SetLineColor(kRed);
    g5->SetLineWidth(lw);
    mg->Add(g5);
    l->AddEntry(g5, "UFC Streuung", "p");

    TGraphErrors *g6 = new TGraphErrors(NumCh, X3, F2, Xerr, DF2);
    g6->SetName("PuFC_SimF");
    g6->SetMarkerStyle(21);
    g6->SetMarkerSize(2);
    g6->SetMarkerColor(kRed);
    g6->SetLineColor(kRed);
    g6->SetLineWidth(lw);
    mg->Add(g6);
    l->AddEntry(g6, "UFC Korrektur", "p");

    mg->Draw("AP");
    mg->GetXaxis()->SetNdivisions(8);
    mg->GetXaxis()->SetLabelSize(0.06);
    mg->GetXaxis()->SetTitleSize(0.08);
    mg->GetXaxis()->SetTitleOffset(0.6);
    mg->GetYaxis()->SetLabelSize(0.06);
    mg->GetYaxis()->SetTitleSize(0.08);
    mg->GetYaxis()->SetTitleOffset(0.7);
    l->Draw();
    c1->Modified();
    c1->Update();
    cout << "Done: Draw correction factor simulation for PuFC and UFC" << endl;
}


void Plot::CompT(Double_t *pExpT, Double_t *DpExpT, Double_t *pSimT, Double_t *DpSimT)
{
    cout << endl << "Plotting transmission comparison..." << endl;
    Double_t X[] = {1,2,3,4,5,6,7,8};
    Double_t Xerr[] = {0,0,0,0,0,0,0,0};
    char name[64] = "";

    // normalization
    Double_t ExpT[NumCh];
    Double_t DExpT[NumCh];
    Double_t norm = pSimT[7] / pExpT[7];
    cout << " Normalization factor " << norm << endl;
    for (Int_t i = 0; i < NumCh; i++)
    {
        ExpT[i] = pExpT[i] * norm;
        DExpT[i] = DpExpT[i] * norm;
    }

    sprintf(name, "%s_CompT", FC.c_str());
    TCanvas * c1 = new TCanvas(name, "Compare transmission", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg = new TMultiGraph();
    sprintf(name, "%s Transmission; Channel; T", FC.c_str());
    mg->SetTitle(name);
    TLegend* l = new TLegend(0.15, 0.4, 0.4, 0.6);

    TGraphErrors *g1 = new TGraphErrors(NumCh, X, ExpT, Xerr, DExpT);
    g1->SetName("ExpT");
    g1->SetLineColor(kRed);
    g1->SetLineWidth(lw);
    mg->Add(g1);
    l->AddEntry(g1, "Experiment", "lp");

    TGraphErrors *g2 = new TGraphErrors(NumCh, X, pSimT, Xerr, DpSimT);
    g2->SetName("SimT");
    g2->SetLineColor(kGreen);
    g2->SetLineWidth(lw);
    mg->Add(g2);
    l->AddEntry(g2, "Simulation", "lp");

    mg->Draw("AP");
    l->Draw();
    c1->Modified();
    c1->Update();
    cout << "Done: Transmission comparison" << endl;
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
    g1->SetMarkerColor(kRed);
    g1->SetMarkerStyle(21);
    g1->SetMarkerSize(0.5);
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


void Plot::Dt(Int_t ch, TH1I *pH, Double_t nf, Double_t Dnf, Int_t l0, Int_t l1, Int_t l2, Int_t l3, string Setup)
{
//    cout << "Plotting Dt..." << endl << "(n,f) " << nf << endl;
    char name[64] = "";
    sprintf(name, "Dt_%i", ch+1);
    TCanvas *C = new TCanvas(name, "Time differecne spectrum", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    Int_t ff = pH->Integral();
    Double_t level = (pH->Integral(l0, l3) - pH->Integral(l1, l2)) / (l1 - l0 + l3 - l2);
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
    sprintf(name, "(n,f): %.0f +- %.0f", nf, Dnf);
    TText *tnf = new TText();
    tnf->SetNDC();
    tnf->DrawText(0.5, 0.8, name);
    sprintf(name, "SF background: %.2f events / ns", level / 4.0);
    TText *tbg = new TText();
    tbg->SetNDC();
//    tbg->DrawText(0.2, 0.2, name);
    sprintf(name, "%i Spaltfragmente", ff);
    TText *tff = new TText();
    tff->SetNDC();
    tff->DrawText(0.4, 0.9, name);
    C->Modified();
    C->Update();
}


void Plot::Stability(Int_t i, Int_t nHist, Double_t *rNIF, Double_t *DrNIF)
{
    Double_t X[nHist];
    Double_t Xerr[nHist];
    for (Int_t j = 0; j < nHist; j++)
    {
        X[j] = j;
        Xerr[j] = 0;
    }
    char name[64] = "";
    sprintf(name, "cStab_%i", i+1);
    TCanvas *C = new TCanvas(name, "Stability", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TGraphErrors *g = new TGraphErrors(nHist, X, rNIF, Xerr, DrNIF);
    sprintf(name, "gStab_%i", i+1);
    g->SetName(name);
    g->SetTitle("(n,f)-Rate Stabilitaet; File nr; Rate [1/s]");
    g->SetLineColor(kRed);
    g->SetLineWidth(lw);
    g->Draw("ap");
    C->Modified();
    C->Update();
}


void Plot::TvsE(Int_t ch, TH2F *pH2TvsE, Double_t binToFmin, Double_t binToFmax, Int_t directN)
{
//    cout << "Going to plot 2D E-t-hist for channel " << ch+1 << endl;
    Double_t ToFmin = pH2TvsE->GetXaxis()->GetBinLowEdge(binToFmin);
    Double_t ToFmax = pH2TvsE->GetXaxis()->GetBinLowEdge(binToFmax + 1);
    Int_t NbinsE = pH2TvsE->GetNbinsY();
    Double_t ScatteredN = pH2TvsE->Integral(binToFmin, binToFmax, 0, NbinsE+1) - directN;

    char name[64] = "";
    sprintf(name, "%s_%s_ToFvsEkin_Ch.%i", FC.c_str(), Setup.c_str(), ch+1);
    TCanvas *C = new TCanvas(name, name, 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    gPad->SetLogz(1);
//    gStyle->SetPalette(kRainBow);
    pH2TvsE->Rebin2D(1, 8);
    sprintf(name, "Kanal %i; #font[12]{t} [ns]; #font[12]{E} [MeV]", ch+1);
    pH2TvsE->SetTitle(name);
    pH2TvsE->SetTitleSize(0.08, "t");
    pH2TvsE->SetStats(0);
    pH2TvsE->GetXaxis()->SetLabelSize(0.06);
    pH2TvsE->GetXaxis()->SetTitleSize(0.07);
    pH2TvsE->GetXaxis()->SetTitleOffset(0.7);
    pH2TvsE->GetYaxis()->SetLabelSize(0.06);
    pH2TvsE->GetYaxis()->SetTitleSize(0.07);
    pH2TvsE->GetYaxis()->SetTitleOffset(0.6);
    pH2TvsE->GetZaxis()->SetLabelSize(0.06);
    pH2TvsE->GetZaxis()->SetLabelOffset(-0.01);

    Double_t xmax = pH2TvsE->GetXaxis()->GetBinLowEdge(pH2TvsE->GetNbinsX()+1);
    Double_t ymax = pH2TvsE->GetYaxis()->GetBinLowEdge(pH2TvsE->GetNbinsY()+1);
    TLine *lLeft = new TLine(ToFmin, 0, ToFmin, ymax);
    lLeft->SetLineColor(kRed);
    TLine *lRight = new TLine(ToFmax, 0, ToFmax, ymax);
    lRight->SetLineColor(kRed);

    pH2TvsE->Draw("colz");
//    lLeft->Draw();
//    lRight->Draw();

    TText *tD = new TText();
    tD->SetNDC();
    sprintf(name, "%i ankommende Neutronen", round_digits(directN, 3));
    tD->DrawText(0.3, 0.8, name);
    TText *tS = new TText();
    tS->SetNDC();
//    char c = '%';
    sprintf(name, "%.0f %% ungestreut", 100.0 * directN / (ScatteredN + (Double_t)directN));
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
    sprintf(name, "%i direct neutrons", round_digits(DirectN, 3));
    tD->DrawText(0.3, 0.8, name);
    TText *tS = new TText();
    tS->SetNDC();
    sprintf(name, "%i scattered neutrons", round_digits(ScatteredN, 3));
    tS->DrawText(0.3, 0.5, name);

    C->Modified();
    C->Update();
}


void Plot::Eproj(Int_t ch, TH1F *pSc, TH1F *pUnsc, Int_t Projectiles)
{
    char name[64] = "";
    sprintf(name, "Neutron energy spectrum channel %i", ch+1);
    TCanvas *C = new TCanvas(name, name, 200, 10, 700, 500);
    gPad->SetTicks(1,1);
    gPad->SetLogy(1);
    pUnsc->Add(pSc, -1);

    // maximum manipulation
    Double_t max = pSc->GetMaximum();
    if (max > pUnsc->GetMaximum())
        pUnsc->SetMaximum(1.1 * max);

    pUnsc->SetStats(0);
    pUnsc->GetXaxis()->SetTitle("#font[12]{E} [MeV]");
    pUnsc->GetYaxis()->SetTitle("#font[12]{N}");//_{eff}
    pUnsc->GetXaxis()->SetTitleSize(0.07);
    pUnsc->GetYaxis()->SetTitleSize(0.07);
    pUnsc->GetXaxis()->SetLabelSize(0.05);
    pUnsc->GetYaxis()->SetLabelSize(0.05);
    pUnsc->SetLineColor(kBlue);
    pUnsc->SetLineWidth(lw);
    pSc->SetLineColor(kRed);
    pSc->SetLineWidth(lw);

    sprintf(name, "%i Neutronen gestartet", round_digits(Projectiles, 3));
    TLegend *l = new TLegend(0.4, 0.6, 0.8, 0.9, name);
    Int_t n = pUnsc->Integral();
    sprintf(name, "%i direkt", round_digits(n, 3));
    l->AddEntry(pUnsc, name);
    n = pSc->Integral();
    sprintf(name, "%i gestreut", round_digits(n, 3));
    l->AddEntry(pSc, name);

    pUnsc->Draw("hist");
    pSc->Draw("same hist");
    l->Draw();
    C->Update();
}



void Plot::Eproj(Int_t ch, TH1F *pSc, TH1F *pUnsc, Int_t Projectiles, TGraph *gSigma)
{
    // Projekt: Merge Eproj and Sigma
    char name[64] = "";
    sprintf(name, "Neutron energy spectrum channel %i", ch+1);
    TCanvas *C = new TCanvas(name, name, 200, 10, 700, 500);
    gPad->SetTicks(1,0);
    gPad->SetLogy(1);
    pUnsc->Add(pSc, -1);

    // maximum manipulation
    Double_t max = pSc->GetMaximum();
    if (max > pUnsc->GetMaximum())
        pUnsc->SetMaximum(1.1 * max);
    max = pUnsc->GetMaximum();

    pUnsc->SetStats(0);
    pUnsc->GetXaxis()->SetTitle("#font[12]{E} / MeV");
    pUnsc->GetYaxis()->SetTitle("#font[12]{N}");
    pUnsc->GetXaxis()->SetTitleSize(0.07);
    pUnsc->GetYaxis()->SetTitleSize(0.07);
    pUnsc->GetXaxis()->SetLabelSize(0.05);
    pUnsc->GetYaxis()->SetLabelSize(0.05);
    pUnsc->SetLineColor(kBlue);
    pUnsc->SetLineWidth(lw);
    pSc->SetLineColor(kRed);
    pSc->SetLineWidth(lw);

    Double_t scale = 1;
    TGaxis *axis = new TGaxis(pUnsc->GetXaxis()->GetXmax(), pUnsc->GetYaxis()->GetXmin() / scale, pUnsc->GetXaxis()->GetXmax(), pUnsc->GetYaxis()->GetXmax() / scale, 0, 1, 509,"+L");
//    axis->Set

    sprintf(name, "%i Neutronen gestartet", round_digits(Projectiles, 3));
    TLegend *l = new TLegend(0.4, 0.6, 0.8, 0.8, name);
    Int_t n = pUnsc->Integral();
    sprintf(name, "%i direkt", round_digits(n, 3));
    l->AddEntry(pUnsc, name);
    n = pSc->Integral();
    sprintf(name, "%i gestreut", round_digits(n, 3));
    l->AddEntry(pSc, name);

    pUnsc->Draw("hist");
    pSc->Draw("same hist");
    gSigma->Draw("same");
    axis->Draw();
    l->Draw();
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
    pUnsc->Add(pSc, -1);

    // maximum manipulation
    Double_t max = pSc->GetMaximum();
    if (max > pUnsc->GetMaximum())
        pUnsc->SetMaximum(1.1 * max);

    pUnsc->SetStats(0);
    pUnsc->GetXaxis()->SetTitle("#font[12]{ToF} / ns");
    pUnsc->GetYaxis()->SetTitle("#font[12]{N}");
    pUnsc->GetXaxis()->SetTitleSize(0.07);
    pUnsc->GetYaxis()->SetTitleSize(0.07);
    pUnsc->GetXaxis()->SetLabelSize(0.05);
    pUnsc->GetYaxis()->SetLabelSize(0.05);
    pUnsc->SetLineColor(kBlue);
    pUnsc->SetLineWidth(lw);
    pSc->SetLineColor(kRed);
    pSc->SetLineWidth(lw);

    sprintf(name, "%i neutrons", round_digits(Projectiles, 3));
    TLegend *l = new TLegend(0.4, 0.6, 0.8, 0.8, name);
    sprintf(name, "%i eff. unsc.", round_digits(pUnsc->Integral(), 3));
    l->AddEntry(pUnsc, name);
    sprintf(name, "%i eff. sc.", round_digits(pSc->Integral(), 3));
    l->AddEntry(pSc, name);
    pUnsc->Draw();
    pSc->SetLineColor(kRed);
    pSc->Draw("same");
    l->Draw();
    C->Update();
}


void Plot::Result(Double_t *uCS, Double_t *DuCS, Double_t *CS, Double_t *DCS)
{
    cout << "Check!" << endl;
    Double_t X1[] = {1,2,3,4,5,6,7,8};
    Double_t X2[] = {1.05,2.05,3.05,4.05,5.05,6.05,7.05,8.05};
    Double_t Xerr[] = {0,0,0,0,0,0,0,0};
    char name[64] = "";

    Double_t scale = 1.E22;
    Double_t uCSscaled[NumCh];
    Double_t DuCSscaled[NumCh];
    Double_t CSscaled[NumCh];
    Double_t DCSscaled[NumCh];
    for (int i = 0; i < NumCh; i++)
    {
        uCSscaled[i] = scale * uCS[i];
        DuCSscaled[i] = scale * DuCS[i];
        CSscaled[i] = scale * CS[i];
        DCSscaled[i] = scale * DCS[i];
    }

    sprintf(name, "c_%s_CS", FC.c_str());
    TCanvas * c1 = new TCanvas(name, "Cross section result", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg = new TMultiGraph();
    sprintf(name, "g_%s_Res", FC.c_str());
    mg->SetNameTitle(name, "N-induzierter Spaltquerschnitt; Kanal; #sigma [b]");
    TLegend* l = new TLegend(0.15, 0.4, 0.4, 0.6);

    TGraphErrors *g1 = new TGraphErrors(NumCh, X1, uCSscaled, Xerr, DuCSscaled);
    g1->SetName("uCS");
    g1->SetMarkerStyle(20);
    g1->SetMarkerColor(kRed);
    g1->SetLineColor(kRed);
    g1->SetLineWidth(lw);
    mg->Add(g1);
    l->AddEntry(g1, "unkorrigiert", "lp");

    TGraphErrors *g2 = new TGraphErrors(NumCh, X2, CSscaled, Xerr, DCSscaled);
    g2->SetName("CS");
    g2->SetMarkerStyle(21);
    g2->SetMarkerColor(kGreen);
    g2->SetLineColor(kGreen);
    g2->SetLineWidth(lw);
    mg->Add(g2);
    l->AddEntry(g2, "T, S korrigiert", "lp");

    mg->Draw("AP");
    mg->GetXaxis()->SetNdivisions(8);
    mg->GetXaxis()->SetLabelSize(0.06);
    mg->GetXaxis()->SetTitleSize(0.08);
    mg->GetXaxis()->SetTitleOffset(0.7);
    mg->GetYaxis()->SetLabelSize(0.06);
    mg->GetYaxis()->SetTitleSize(0.08);
    mg->GetYaxis()->SetTitleOffset(0.8);
    l->Draw();
    c1->Modified();
    c1->Update();
}


void Plot::CalibrateUFC(Double_t* emA, Double_t* DemA, Double_t* Area, Double_t* DArea)
{
    Double_t X[] = {1, 2, 3, 4, 5, 6, 7, 8};
    Double_t Xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};
    Double_t ema[] = {365.2, 396.2, 400.4, 393, 403.9, 403.7, 406.6, 397.1}; // calibration on H19
    Double_t Dema[] = {1.7, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8};

    TCanvas *c1 = new TCanvas("c_U_Cal", "Charakterisierung der Uran-Deposits", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg = new TMultiGraph();
    mg->SetNameTitle("mg_U_Cal", "Charakterisierung der Uran-Deposits; Kanal; #varepsilon#font[12]{m_{A}} [#mug/cm^{2}]");
    TLegend* l = new TLegend(0.6, 0.15, 0.8, 0.4);

    TGraphErrors *g1 = new TGraphErrors(NumCh, X, emA, Xerr, DemA);
    g1->SetName("PTB");
    g1->SetMarkerStyle(20);
    g1->SetMarkerSize(1.5);
    g1->SetMarkerColor(kRed);
    g1->SetLineColor(kRed);
    g1->SetLineWidth(lw);
    mg->Add(g1);
    l->AddEntry(g1, "PTB");

    TGraphErrors *g2 = new TGraphErrors(NumCh, X, ema, Xerr, Dema);
    g2->SetName("H19");
    g2->SetMarkerStyle(20);
    g2->SetMarkerSize(1.5);
    g2->SetMarkerColor(kGreen);
    g2->SetLineColor(kGreen);
    g2->SetLineWidth(lw);
    mg->Add(g2);
    l->AddEntry(g2, "H19");

    mg->Draw("AP");
    mg->GetXaxis()->SetNdivisions(8);
    mg->GetXaxis()->SetLabelSize(0.06);
    mg->GetXaxis()->SetTitleSize(0.08);
    mg->GetXaxis()->SetTitleOffset(0.7);
    mg->GetYaxis()->SetLabelSize(0.06);
    mg->GetYaxis()->SetTitleSize(0.08);
    mg->GetYaxis()->SetTitleOffset(0.8);
    l->Draw();
    c1->Modified();
    c1->Update();

    TCanvas *c2 = new TCanvas("c_U_Area", "Efektive Flaeche", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);

    TGraphErrors *g3 = new TGraphErrors(NumCh, X, Area, Xerr, DArea);
    g3->SetNameTitle("Area", "Effektive Fl#ddot{a}che; Kanal; #font[12]{A} [mm^{2}]");
    g3->SetMarkerStyle(20);
    g3->SetMarkerSize(1.5);
    g3->SetMarkerColor(kBlue);
    g3->SetLineColor(kBlue);
    g3->SetLineWidth(lw);

    g3->Draw("AP");
    g3->GetXaxis()->SetNdivisions(8);
    g3->GetXaxis()->SetLabelSize(0.06);
    g3->GetXaxis()->SetTitleSize(0.08);
    g3->GetXaxis()->SetTitleOffset(0.7);
    g3->GetYaxis()->SetLabelSize(0.06);
    g3->GetYaxis()->SetTitleSize(0.08);
    g3->GetYaxis()->SetTitleOffset(0.8);
    c2->Modified();
    c2->Update();
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


void Plot::SaveToFile(string RootPath, TObject *pObj)
{   //saves a TObject into member TFile *f
    f->ReOpen("UPDATE");
    TDirectory *EvalDir;
    TObject *pGraph;
    pGraph = (TObject*) pObj;
    string GraphName = pGraph->GetName();
    //check if folder path already exists, otherwise create it
    if (f->Get(RootPath.c_str())!=0)
        EvalDir = (TDirectory*) f->Get(RootPath.c_str());
    else
    {   cout << " Creating root directory " << RootPath << endl;
        f->mkdir(RootPath.c_str(), "Folder");
        EvalDir = f->GetDirectory(RootPath.c_str());
    }
    EvalDir->cd();
    if (EvalDir->Get(GraphName.c_str())!=0)
        EvalDir->Delete((GraphName+";*").c_str());
    pGraph->Clone()->Write();
    f->Save(); //file->Close();
    cout << " Saved " << GraphName << endl;
}


Int_t Plot::round_digits(Double_t x, Int_t digits)
{
    Int_t max_digit = floor(log10(x));
    Int_t e = max_digit - digits + 1;
    Double_t f = pow(10, e);
//    cout << x << "  " << digits << "  " << f * floor(x/f+0.5) << endl;
    return f * floor(x / f + 0.5);
}
