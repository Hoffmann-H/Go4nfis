//#include "/gpfs/home/hoffma93/StyleSheets/StyleSheet.C"
//#include "TLatex.h"
#include "TTree.h"
#include "nELBEsim.C"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
//#include "DrawPics.C"
using namespace std;

void DrawVacuum(TFile *f, TFile *g, string FC, Int_t ch)
{
    char name[64] = "";
    sprintf(name, "Hit/CorrLoss/H_k_%s_%i", FC.c_str(), ch+1);
    TH1F *hK = (TH1F*)f->Get(name); if (!hK) cout << "Could not get " << name << endl;
    sprintf(name, "Hit/Transmission/H_T_%s_%i", FC.c_str(), ch+1);
    TH1F *hT = (TH1F*)f->Get(name); if (!hT) cout << "Could not get " << name << endl;
    sprintf(name, "Hit/Correction/H_C_%s_%i", FC.c_str(), ch+1);
    TH1F *hC = (TH1F*)g->Get(name); if (!hC) cout << "Could not get " << name << endl;
    sprintf(name, "H_Div_%s_%i", FC.c_str(), ch+1);
    TH1F *hDiv = (TH1F*)hK->Clone();
    hDiv->Divide(hK, hT, 1.0, 1.0);

    hDiv->SetStats(0);
    hDiv->GetYaxis()->SetTitle("correction factor #font[12]{C}");
    hDiv->GetXaxis()->SetRangeUser(0.1, 20.0);
    hDiv->GetYaxis()->SetRangeUser(0., 1.8);
    hDiv->SetLineColor(kBlue);
    hC->SetLineColor(kRed);
    sprintf(name, "%s channel %i", FC.c_str(), ch+1);
    TLegend *l = new TLegend(0.6, 0.7, 0.9, 0.9, name);
    l->AddEntry(hDiv, "#font[12]{k} / #font[12]{T}");
    l->AddEntry(hC, "#font[12]{C}");

    TCanvas *cV = new TCanvas("cV");
    gPad->SetLogx();
    gPad->SetTicks(1,1);
    hDiv->Draw("hist");
    hC->Draw("same hist");
    l->Draw();
    cV->Update();
}

void DrawFissionEE(TFile *f, string FC, Int_t ch, Int_t sc, Bool_t wcs)
{ /// Draw a single E vs E(t) 2D histogam.
  /// Choose fission chamber, channel, total/scattered, cross section weighting, total(0)/scattered(1)/direct(2).
  /// Beautiful color map, histogram label
    TH2F *h2;
    string TotPath = FC+"/EToFvsEkin/"+FC+"_EToFvsEkin_Ch."+std::to_string(ch+1);
    string ScPath = FC+"/EToFvsEkin/Scattered/"+FC+"_EToFvsEkin_Sc_Ch."+std::to_string(ch+1);
    if (sc == 0)
        h2 = (TH2F*)f->Get(TotPath.c_str()); if (!h2) cout << "Could not get " << TotPath << endl;
    if (sc == 1)
        h2 = (TH2F*)f->Get(ScPath.c_str()); if (!h2) cout << "Could not get " << ScPath << endl;
    if (sc == 2) {
        h2 = (TH2F*)f->Get(TotPath.c_str()); if (!h2) cout << "Could not get " << TotPath << endl;
        TH2F *h2Sc = (TH2F*)f->Get(ScPath.c_str()); if (!h2Sc) cout << "Could not get " << ScPath << endl;
        h2->Add(h2Sc, -1.0);
    }

    // weight with cross section
    TGraphErrors *pGrXS = GetXS("~/TrackLength/FissionXS.root", FC, "JEFF_3.3");
    if (wcs) h2 = WeightWithXS(h2, pGrXS);

    CreateColorGradient();

    TCanvas *cEE = new TCanvas(("cEE"+std::to_string(ch+1)).c_str());
    gPad->SetLogz();
    gPad->SetTicks(1,1);
    gPad->SetRightMargin(0.13);
    h2->SetStats(0);
    h2->GetXaxis()->SetRangeUser(0., 20.);
    h2->GetYaxis()->SetRangeUser(0., 20.);
    h2->SetMinimum(0.5);
    h2->SetMaximum(3e5);
    h2->Draw("COLZ");

    TLatex *Tl = new TLatex(); Tl->SetTextSize(0.07);
    string text = FC+", ch."+std::to_string(ch+1);
    if (sc == 0) text = text + ", sc+unsc";
    if (sc == 1) text = text + ", sc";
    if (sc == 2) text = text + ", unsc";
    if (wcs) text = text + ", #sigma(#font[12]{E})";
    Tl->DrawLatexNDC(0.2, 0.82, ("#font[132]{"+text+"}").c_str());
    gPad->Modified();
    gPad->Update();
}

void DrawFissionEt(TFile *f, string FC, Int_t ch, Int_t sc, Bool_t wcs)
{ /// Draw a single E vs t 2D histogam.
  /// Choose fission chamber, channel, total/scattered, cross section weighting, total(0)/scattered(1)/direct(2).
  /// Beautiful color map, histogram label
    TH2F *h2;
    string TotPath = FC+"/ToFvsEkin/"+FC+"_ToFvsEkin_Ch."+std::to_string(ch+1);
    string ScPath = FC+"/ToFvsEkin/Scattered/"+FC+"_ToFvsEkin_Sc_Ch."+std::to_string(ch+1);
    if (sc == 0)
        h2 = (TH2F*)f->Get(TotPath.c_str()); if (!h2) cout << "Could not get " << TotPath << endl;
    if (sc == 1)
        h2 = (TH2F*)f->Get(ScPath.c_str()); if (!h2) cout << "Could not get " << ScPath << endl;
    if (sc == 2) {
        h2 = (TH2F*)f->Get(TotPath.c_str()); if (!h2) cout << "Could not get " << TotPath << endl;
        TH2F *h2Sc = (TH2F*)f->Get(ScPath.c_str()); if (!h2Sc) cout << "Could not get " << ScPath << endl;
        h2->Add(h2Sc, -1.0);
    }

    // weight with cross section
    TGraphErrors *pGrXS = GetXS("~/TrackLength/FissionXS.root", FC, "JEFF_3.3");
    if (wcs) h2 = WeightWithXS(h2, pGrXS);

    CreateColorGradient();

    TCanvas *cEt = new TCanvas(("cEt"+std::to_string(ch+1)).c_str());
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gPad->SetTicks(1,1);
    gPad->SetRightMargin(0.13);
    h2->SetStats(0);
    h2->GetXaxis()->SetRangeUser(0., 20.);
    h2->GetYaxis()->SetRangeUser(0., 20.);
    h2->SetMinimum(0.5);
    h2->SetMaximum(3e5);
    h2->Draw("COLZ");

    TLatex *Tl = new TLatex(); Tl->SetTextSize(0.07);
    string text = FC+", ch."+std::to_string(ch+1);
    if (sc == 0) text = text + ", sc+unsc";
    if (sc == 1) text = text + ", sc";
    if (sc == 2) text = text + ", unsc";
    if (wcs) text = text + ", #sigma(#font[12]{E})";
    Tl->DrawLatexNDC(0.4, 0.82, ("#font[132]{"+text+"}").c_str());
    gPad->Modified();
    gPad->Update();
}

void Draw()
{
    LoadStyles();
    gROOT->SetStyle("SinglePadStyle");
    gROOT->ForceStyle(kTRUE);
    gStyle->SetLegendFont(132);
    TFile *f = TFile::Open("/home/hoffma93/TrackLength/G4UFCvsH19_1E7_hist_result.root", "READ");
    TFile *g = TFile::Open("/home/hoffma93/TrackLength/G4UFCvsH19_1E7_hist_result.root", "READ");
//    DrawTrackLength(f, "UFC", 0);
    DrawVacuum(f, f, "H19", 0);
}

void DrawWeight(TFile *f, string xvar = "ToF")
{ // xvar = "ToF", "Ekin"
    LoadStyles();
    gROOT->SetStyle("SinglePadStyle");
    gROOT->ForceStyle(kTRUE);
    gStyle->SetLegendFont(132);

    string HistName[4] = {"Weight/H19_Wvs"+xvar+"_Ch.1",
                          "Weight/Scattered/H19_Wvs"+xvar+"_Sc_Ch.1",
                          "Weight/UFC_Wvs"+xvar+"_Ch.1",
                          "Weight/Scattered/UFC_Wvs"+xvar+"_Sc_Ch.1"};
    string Label[4] = {"H19, ch.1, tot.",
                       "H19, ch.1, sc.",
                       "UFC, ch.1, tot.",
                       "UFC, ch.1, sc."};
    Int_t color[4] = {2, 3, 4, 6};

    TH1F *h[4];
    TLegend *l = new TLegend(0.3, 0.2, 0.8, 0.4);

    for (Int_t i = 0; i < 4; i++)
    {
        h[i] = (TH1F*)f->Get(HistName[i].c_str());
        if (!h[i]) cout << "Could not get " << HistName[i] << endl;
//        h[0]->GetXaxis()->SetRangeUser(0.1, 20.);
        h[i]->GetYaxis()->SetTitleSize(0.07);
        h[i]->GetYaxis()->SetLabelSize(0.07);
        h[i]->SetLineColor(color[i]);
        h[i]->SetLineWidth(2);
        l->AddEntry(h[i], Label[i].c_str());
    }
    h[0]->GetXaxis()->SetRangeUser(0.1, 20.);
    h[0]->GetYaxis()->SetRangeUser(0., 2.5);

    new TCanvas();
    gPad->SetTicks(1,1);
    gPad->SetLogx();
    h[0]->Draw("");
    h[1]->Draw("same ");
    h[2]->Draw("same ");
    h[3]->Draw("same ");
    l->Draw();
}

void DrawCF(string StrFileName)
{
    TFile *f = TFile::Open(StrFileName.c_str(), "READ");
    if (!f) cout << "Could not open " << StrFileName << endl;

    string HistName[4] = {"Hit/Correction/H_C_H19_1",
                          "TrackLength/Correction/W_C_H19_1",
                          "Hit/Correction/H_C_UFC_1",
                          "TrackLength/Correction/W_C_UFC_1"};
    string Label[4] = {"H19, ch.1, old",
                       "H19, ch.1, new",
                       "UFC, ch.1, old",
                       "UFC, ch.1, new"};
    Int_t color[4] = {6, 7, 8, 9};

    TH1F *h[4];
    TLegend *l = new TLegend(0.3, 0.2, 0.8, 0.4);

    for (Int_t i = 0; i < 4; i++)
    {
        h[i] = (TH1F*)f->Get(HistName[i].c_str());
        h[0]->GetXaxis()->SetRangeUser(0.1, 20.);
        h[i]->GetYaxis()->SetTitleSize(0.07);
        h[i]->GetYaxis()->SetLabelSize(0.07);
        h[i]->SetLineColor(color[i]);
        h[i]->SetLineWidth(2);
        l->AddEntry(h[i], Label[i].c_str());
    }
    h[0]->GetYaxis()->SetTitle("Correction factor #font[12]{C}");
    h[0]->GetXaxis()->SetRangeUser(0.1, 20.);
    h[0]->GetYaxis()->SetRangeUser(0.8, 1.3);

    new TCanvas();
    gPad->SetTicks(1,1);
    gPad->SetLogx();
    h[0]->Draw("hist");
    h[1]->Draw("same hist");
    h[2]->Draw("same hist");
    h[3]->Draw("same hist");
    l->Draw();
}

void DrawRefH19(string StrFileName, Int_t ch)
{
    char name[64] = "";
    TFile *f = TFile::Open(StrFileName.c_str(), "READ"); if (!f) cout << "Could not open " << StrFileName << endl;
    TFile *g = TFile::Open("~/FC-Analysis/2016.06/results/nfis_2016.06_orig.root", "READ"); if (!g) cout << "Could not open nfis_2016.06.root" << endl;

    sprintf(name, "Hit/Correction/RefH19/H_C_UFC_RefH19_%i", ch+1);
    TH1F *hH = (TH1F*)f->Get(name);
    sprintf(name, "TrackLength/Correction/RefH19/W_C_UFC_RefH19_%i", ch+1);
    TH1F *hW = (TH1F*)f->Get(name);
    sprintf(name, "Analysis/Simulation/Geant4/Correction/G4Correction_UFC_Ch.%i", ch+1);
    TH1D *hO = (TH1D*)g->Get(name); if (!hO) cout << "Could not get " << name << endl;
    TH1D *hH19 = (TH1D*)g->Get("Analysis/Simulation/Geant4/Correction/G4Correction_H19"); if (!hH19) cout << "Could not get H19" << endl;
    hO->Divide(hH19);

    hH->GetYaxis()->SetTitleSize(0.07);
    hH->GetYaxis()->SetLabelSize(0.07);
    hH->SetLineColor(6);
    hW->SetLineColor(7);
    hO->SetLineColor(8);
//    hH->SetLineWidth(2);
    hH->GetYaxis()->SetTitle("Correction factor #font[12]{C}");
    hH->GetXaxis()->SetRangeUser(0.1, 20.);
    hH->GetYaxis()->SetRangeUser(0.8, 1.4);

    sprintf(name, "UFC, ch.%i, ref. H19", ch+1);
    TLegend *l = new TLegend(0.3, 0.2, 0.8, 0.5, name);
    l->AddEntry(hO, "2016");
    l->AddEntry(hH, "2020 Hit Count");
    l->AddEntry(hW, "2020 Track Length");

    new TCanvas();
    gPad->SetTicks(1,1);
    gPad->SetLogx();
    hH->Draw("hist");
    hW->Draw("same hist");
    hO->Draw("same hist");
    l->Draw();
}

void DrawWspectrum(string StrTreeFile = "~/TrackLength/G4UFCvsH19_1E7_tree.root", string FC = "H19", Int_t ch = 0)
{
    string Name = "";
    TFile *f = TFile::Open(StrTreeFile.c_str(), "READ"); if (!f) cout << "Could not open " << StrTreeFile << endl;
    Name = FC+"/Channel_"+std::to_string(ch+1);
    TTree *t = (TTree*)f->Get(Name.c_str()); if (!t) cout << "Could not get " << "H19/CHannel_1" << endl;
    TCanvas *cW = new TCanvas("cW");

    // Total W spectrum
    t->Draw("Weight>>htemp(1000,0,100)");
    TH1F *hTot = (TH1F*)gPad->GetPrimitive("htemp"); if (!hTot) cout << "Could not get " << "htemp" << endl;
    Name = "Wspec_tot_"+FC+std::to_string(ch+1);
    hTot->SetName(Name.c_str());
//    hTot->SetDirectory(0);
    hTot->SetStats(0);
    hTot->SetTitle("; #font[12]{l/V}; #font[12]{N}");
    hTot->SetLineWidth(2);

    // All scattered neutrons' W spectrum
    t->Draw("Weight>>htemp(1000,0,100)", "GunEnergy>En", "same");
    TH1F *hScAll = (TH1F*)gPad->GetPrimitive("htemp"); if (!hScAll) cout << "Could not get " << "htemp" << endl;
    Name = "Wspec_sc_"+FC+std::to_string(ch+1);
    hScAll->SetName("hScAll");
    hScAll->SetLineColor(kOrange);

    // E-dependent scattered neutrons' W spectra
//    TH1F* hSc[4];
//    t->Draw("Weight>>htemp(1000,0,100)", "GunEnergy>En && En<1.0", "same");
//    hSc[0] = (TH1F*)gPad->GetPrimitive("htemp"); if (!hSc[0]) cout << "Could not get " << "htemp" << endl;
//    Name = "Wspec_sc0"+FC+std::to_string(ch+1);
//    hSc[0]->SetName(Name.c_str());
//    hSc[0]->SetLineColor(kRed);
//    t->Draw("Weight>>htemp(1000,0,100)", "GunEnergy>En && En>1.0 && En<2.0", "same");
//    hSc[1] = (TH1F*)gPad->GetPrimitive("htemp"); if (!hSc[1]) cout << "Could not get " << "htemp" << endl;
//    Name = "Wspec_sc1"+FC+std::to_string(ch+1);
//    hSc[1]->SetName(Name.c_str());
//    hSc[1]->SetLineColor(kGreen);
//    t->Draw("Weight>>htemp(1000,0,100)", "GunEnergy>En && En>2.0", "same");
//    hSc[2] = (TH1F*)gPad->GetPrimitive("htemp"); if (!hSc[2]) cout << "Could not get " << "htemp" << endl;
//    Name = "Wspec_sc2"+FC+std::to_string(ch+1);
//    hSc[2]->SetName(Name.c_str());
//    hSc[2]->SetLineColor(kBlue);

    // legend
    Name = FC+", ch."+std::to_string(ch+1);
    TLegend *l = new TLegend(0.7, 0.6, 0.9, 0.9, Name.c_str());
    l->AddEntry(hTot, "sc+unsc");
//    l->AddEntry(hSc[0], "sc, E<1MeV");
//    l->AddEntry(hSc[1], "sc, 1MeV<E<2MeV");
//    l->AddEntry(hSc[2], "sc, E>2MeV");
    l->AddEntry(hScAll, "all sc");
    l->Draw();

    gPad->SetLogy();
    gPad->SetTicks(1,1);
    cW->Update();
}

void DrawWspectra(string StrTreeFile = "~/TrackLength/G4UFCvsH19_1E7_tree.root")
{
    char name[64] = "";
    TFile *f = TFile::Open(StrTreeFile.c_str(), "READ"); if (!f) cout << "Could not open " << StrTreeFile << endl;

    Int_t n = 2;
    string FC[] = {"H19", "UFC"};
    Int_t ch[] = {0, 0};
    TTree *t[n];
    TH1F *h[n];

    TCanvas *cW = new TCanvas("cW");
    for (Int_t i = 0; i < n; i++)
    {
//        Name = FC[i]+"/Channel_"+std::to_string(ch[i]+1);
        sprintf(name, "%s/Channel_%i", FC[i].c_str(), ch[i]+1);
        t[i] = (TTree*)f->Get(name); if (!t[i]) cout << "Could not get " << name << endl;
        t[i]->Draw("Weight>>htemp(1000,0,100)");
        h[i] = (TH1F*)gPad->GetPrimitive("htemp"); if (!h[i]) cout << "Could not get " << "htemp" << endl;
//        Name = "Wspec_"+FC+std::to_string(ch+1);
        sprintf(name, "Wspec_%s_%i", FC[i].c_str(), ch[i]);
        h[i]->SetName(name);
        h[i]->SetLineColor(i+1);
    }
    h[0]->SetTitle("; #font[12]{#Phi} / cm^{2}; #font[12]{n}");
    h[0]->Draw();
    h[1]->Draw("same");
    TLegend *l = new TLegend(0.7, 0.6, 0.9, 0.9);
    for (Int_t i = 0; i < n; i++) {
        sprintf(name,"%s, ch.%i", FC[i].c_str(), ch[i]+1);
        l->AddEntry(h[i], name); }
    l->Draw();
    gPad->SetLogy();
    gPad->SetTicks(1,1);
    cW->Update();
}

void DrawWspectra(TFile *f0, TFile *f1, string FC, Int_t ch, Bool_t sc = 0)
{
    char name[64] = "";

    Int_t n = 2;
    TTree *t[n];
    TH1F *h[3];

    TCanvas *cW = new TCanvas("cW");
    gPad->SetLogx();
    sprintf(name, "%s/Channel_%i", FC.c_str(), ch+1);
    t[0] = (TTree*)f0->Get(name); if (!t[0]) cout << "Could not get " << name << endl;
    t[1] = (TTree*)f1->Get(name); if (!t[1]) cout << "Could not get " << name << endl;

    t[0]->Draw("Weight>>htemp(1000,0,2.325127)", sc?"GunEnergy>En":"");
    h[0] = (TH1F*)gPad->GetPrimitive("htemp"); if (!h[0]) cout << "Could not get " << "htemp" << endl;
    h[0]->SetName("hPhi");
    h[0]->SetLineColor(1);

    t[1]->Draw("Weight>>htemp(1000,0,100)", sc?"GunEnergy>En":"");
    h[1] = (TH1F*)gPad->GetPrimitive("htemp"); if (!h[1]) cout << "Could not get " << "htemp" << endl;
    h[1]->SetName("hTheta");
    h[0]->SetTitle("; #font[12]{#Phi} / cm^{-2}; #font[12]{n}");
    h[2] = (TH1F*) h[0]->Clone();
    h[2]->SetLineColor(2);
    for (Int_t bin = 0; bin < h[1]->GetNbinsX()+2; bin++)
    {
        h[2]->SetBinContent(bin, h[1]->GetBinContent(bin));
        h[2]->SetBinError(bin, h[1]->GetBinError(bin));
    }

    h[0]->GetXaxis()->SetRangeUser(0.01, 1.);
    h[0]->Draw("hist");
    h[0]->Draw("same e");
    h[2]->Draw("same hist");
    h[2]->Draw("same e");
    sprintf(name, "%s, ch.%i%s", FC.c_str(), ch+1, sc?", sc":"");
    TLegend *l = new TLegend(0.7, 0.6, 0.9, 0.9, name);
    l->AddEntry(h[0], "#Phi #propto #font[12]{l_{i}}");
    l->AddEntry(h[2], "#Phi #propto 1/cos(#font[12]{#Theta_{i}})");
    l->Draw();
    gPad->SetLogy();
    gPad->SetTicks(1,1);
    cW->Update();
}

void CompSimRes(TFile *f0, string descr0, TFile *f1, string descr1, string FC, Int_t ch, string plot = "dual", string var = "Ekin", Bool_t sc = 0)
{ /// Plot a 1D variable for 2 simulations, given their hist files and descriptions.
  /// Plot as "dual", quotient "div" or difference "diff".
  /// Variables "Ekin", "ToF".
    Double_t RebinFactor = 1;
    char name[64] = "";
    if (!strcmp(var.c_str(), "Ekin")) {
        RebinFactor = 1;
        if (sc) sprintf(name, "%s/nEnergy/Scattered/%s_Ekin_Sc_Ch.%i", FC.c_str(), FC.c_str(), ch+1);
        else    sprintf(name, "%s/nEnergy/%s_Ekin_Ch.%i", FC.c_str(), FC.c_str(), ch+1);
    } else if (!strcmp(var.c_str(), "ToF")) {
        RebinFactor = 100;
        if (sc) sprintf(name, "%s/ToF/Scattered/%s_ToF_Sc_Ch.%i", FC.c_str(), FC.c_str(), ch+1);
        else    sprintf(name, "%s/ToF/%s_ToF_Ch.%i", FC.c_str(), FC.c_str(), ch+1);
    } else cout << "Unknown variable " << var << endl;
    TH1F *h0 = (TH1F*)f0->Get(name); if (!h0) cout << "Could not get " << name << " from " << f0->GetName() << endl;
    TH1F *h1 = (TH1F*)f1->Get(name); if (!h1) cout << "Could not get " << name << " from " << f1->GetName() << endl;
    h0->SetStats(0);
    h0->SetLineColor(kRed);
    sprintf(name, "%s, ch.%i%s, %s", FC.c_str(), ch+1, sc?", sc":"", var.c_str());
    TLegend *l = new TLegend(0.7, 0.6, 0.9, 0.9, name);
    TCanvas *cS = new TCanvas("cS");
    if (!strcmp(plot.c_str(), "dual"))
    {
        h0->Rebin(RebinFactor);
        h0->Scale(1./RebinFactor);
        h1->Rebin(RebinFactor);
        h1->Scale(1./RebinFactor);
        h0->Draw("hist");
        h1->Draw("same hist");
        l->AddEntry(h0, descr0.c_str());
        l->AddEntry(h1, descr1.c_str());
    } else if (!strcmp(plot.c_str(), "diff"))
    {
        sprintf(name, "#Delta %s", h0->GetYaxis()->GetTitle());
        h0->GetYaxis()->SetTitle(name);
        h0->Add(h1, -1.0);
        h0->Rebin(RebinFactor);
        h0->Scale(1./RebinFactor);
        h0->Draw();
        sprintf(name, "(%s) - (%s)", descr0.c_str(), descr1.c_str());
        l->AddEntry(h0, name);
    } else if (!strcmp(plot.c_str(), "div"))
    {
        sprintf(name, "#font[12]{N}_{%s} / #font[12]{N}_{%s}", descr0.c_str(), descr1.c_str());
        h0->GetYaxis()->SetTitle(name);
        h0->Divide(h1);
        h0->Rebin(RebinFactor);
        h0->Scale(1./RebinFactor);
        h0->Draw();
        sprintf(name, "(%s) / (%s)", descr0.c_str(), descr1.c_str());
        l->AddEntry(h0, name);
    } else cout << "Unknown plotting option " << plot << endl;
    l->Draw();
    cS->Update();
}

void CompProjections(string FC, Int_t ch, Bool_t sc = 0)
{ /// Compare spectra projections for simulations
    string Files[] = {"~/TrackLength/G4UFCvsH19_23082016_hist.root",
                      "~/TrackLength/G4UFCvsH19_v10.2_hist.root",
                      "~/TrackLength/G4UFCvsH19_v10.2_hist_w.root",
                      "~/TrackLength/G4UFCvsH19_1E9_hist.root",
                      "~/TrackLength/G4UFCvsH19_1E9_hist_w.root",
                      "~/TrackLength/MCNP_0313_filled.root"};
    string Descr[] = {"G4 2016", "G4 10.2", "G4 10.2 TL", "G4 10.5", "G4 10.5 TL", "MCNP6"};
    TH2F *h2[6];
    TH1F *hEToF[6];
    TGraphErrors *gEToF[6];
//    TH1D *hEkin[6];
    uint N = 6;
    Double_t RebinFactor = 1;
    char name[64] = "";
    for (Int_t i = 0; i < N; i++)
    {
        TFile *f = TFile::Open(Files[i].c_str());
        TH1I *hSrc = (TH1I*)f->Get("Source/nEnergy/Source_Ekin");
        if (sc) sprintf(name, "%s/EToFvsEkin_wXS/Scattered/%s_EToFvsEkin_wXS_Sc_Ch.%i", FC.c_str(), FC.c_str(), ch+1);
        else    sprintf(name, "%s/EToFvsEkin_wXS/%s_EToFvsEkin_wXS_Ch.%i", FC.c_str(), FC.c_str(), ch+1);
        h2[i] = (TH2F*) f->Get(name); if (!h2[i]) cout << "Could not get " << name << endl;
        h2[i]->SetDirectory(0);
        sprintf(name, "px_%s_%i", FC.c_str(), ch+1);
        hEToF[i] = (TH1F*)h2[i]->ProjectionX(name, 0, -1, "e");
        hEToF[i]->Scale(1.0 / hSrc->Integral());
        hEToF[i]->SetDirectory(0);
        gEToF[i] = TH1toTGraphError(hEToF[i]);
    }
    Double_t A = 43; // [A] = cm^2. UFC.
    hEToF[0]->Scale(1.0 / A);
    hEToF[1]->Scale(1.0 / A);

    hEToF[0]->SetStats(0);
    sprintf(name, "%s, ch.%i%s", FC.c_str(), ch+1, sc?", sc":"");
    TLegend *l = new TLegend(0.7, 0.6, 0.9, 0.9, name);
    TCanvas *cS = new TCanvas("cS");
    for (Int_t i = 0; i < N; i++)
    {
        hEToF[i]->SetLineColor(i+2);
        hEToF[i]->Draw("same hist");
        l->AddEntry(hEToF[i], Descr[i].c_str());
        gEToF[i]->SetLineColor(i+2);
        gEToF[i]->SetFillColorAlpha(i+2, 0.5);
//        gEToF[i]->Draw("same 3");
    }

    l->Draw();
    cS->Update();
}

void CompGraph(string HistoPath, string HistoName, Int_t N = 2)
{
    string Files[] = {"~/TrackLength/G4UFCvsH19_1E7_hist.root",
                      "~/TrackLength/G4UFCvsH19_1E7_hist_w.root", // vorlaeufig Test mit 1E7
                      "~/TrackLength/G4UFCvsH19_23082016_hist.root",
                      "~/TrackLength/G4UFCvsH19_v10.2_hist.root",
                      "~/TrackLength/G4UFCvsH19_v10.2_hist_w.root",
                      "~/TrackLength/G4UFCvsH19_1E9_hist.root",
                      "~/TrackLength/G4UFCvsH19_1E9_hist_w.root",
                      "~/TrackLength/MCNP_0313_filled.root"};
    string Descr[] = {"G4 2016", "G4 10.2", "G4 10.2 TL", "G4 10.5", "G4 10.5 TL", "MCNP6"};
    TCanvas *cGr = new TCanvas("cGr");
    gPad->SetTicks(1, 1);
    for (uint i = 0; i < N; i++)
    {
        TFile *f = TFile::Open(Files[i].c_str(), "READ"); if (!f) cerr << "Could not open " << Files[i] << endl;
        string Name = HistoPath + "Graph/Gr_" + HistoName;
        TGraphErrors *g = (TGraphErrors*) f->Get(Name.c_str()); if (!g) cerr << "Could not get " << Name << endl;

    }
}

void CompCF(string FC, Int_t ch, TFile *f0, string descr0, TFile *f1, string descr1, TFile *f2 = 0, string descr2 = "")
{
/// Plot correction factor C(E(t)) from histogram restult files *_hist.root or *_hist_w.root
/// FC = "H19", "UFC", "UFC_RefH19"
    string HistPath, HistName, Name; // Changed from hists to TGraphs!
    if (strcmp(FC.c_str(), "UFC") && strcmp(FC.c_str(), "H19")) {
        HistPath = "UFC_RefH19/";
        HistName = "C_UFC_RefH19_" + std::to_string(ch+1);
    } else {
        HistPath = FC + "/Correction/";
        HistName = "G4Correction_" + FC + "_Ch." + std::to_string(ch+1);
    }
    TCanvas *cTL = new TCanvas("cTL");
    gPad->SetLogx();
    gPad->SetTicks(1,1);
    Name = FC + ", ch." + std::to_string(ch+1);
    TLegend *l = new TLegend(0.6, 0.6, 0.9, 0.9, Name.c_str());

    TMultiGraph *mg = new TMultiGraph();
    mg->GetXaxis()->SetRangeUser(0.1, 20.0);
//    mg->GetYaxis()->SetRangeUser(0.8, 1.5);
    mg->GetYaxis()->SetTitle("correction factor #font[12]{C}");

    Name = HistPath + HistName;
    TGraphErrors *g0 = (TGraphErrors*)f0->Get(Name.c_str()); if (!g0) cout << "Could not get " << Name << " from " << f0->GetName() << endl;
    g0->SetLineColor(kBlue);
    g0->SetFillColorAlpha(kBlue, 0.5);
    mg->Add(g0);
    l->AddEntry(g0, descr0.c_str());

    Name = HistPath + HistName;
    TGraphErrors *g1 = (TGraphErrors*)f1->Get(Name.c_str()); if (!g1) cout << "Could not get " << Name << " from " << f1->GetName() << endl;
    g1->SetLineColor(kRed);
    g1->SetFillColorAlpha(kRed, 0.5);
    mg->Add(g1);
    l->AddEntry(g1, descr1.c_str());

    if (f2) {
        Name = HistPath + HistName;
        TGraphErrors *g2 = (TGraphErrors*)f2->Get(Name.c_str()); if (!g2) cout << "Could not get " << Name << " from " << f2->GetName() << endl;
        g2->SetLineColor(kRed);
        g2->SetFillColorAlpha(kRed, 0.5);
        mg->Add(g2);
        l->AddEntry(g2, descr2.c_str());
    }

    mg->Draw("ae0");
    l->Draw();
    cTL->Update();
}

void CompCF_UFC_H19(Int_t var, TFile *f0, string descr0, TFile *f1, string descr1, string descrFile = "TL")
{
/// Plot correction factor C(E(t)) for H19 ch 1 and UFC ch 1 together
    char name[64] = "";
    string ArrStrVar[] = {"TransMis", "Scattering", "Correction"};
    string ArrStrVarName[] = {"Transmission factor  #it{T}",
                            "Correlation factor  #it{k}",
                            "Correction factor  #it{C}"};
    string ArrStrSymbol[] = {"T", "k", "C"};

    sprintf(name, "H19/%s/G4%s_H19_Ch.1", ArrStrVar[var].c_str(), ArrStrVar[var].c_str());
    TGraphErrors *g0H19 = (TGraphErrors*)f0->Get(name); if (!g0H19) cerr << "Could not get " << name << endl;
    TGraphErrors *g1H19 = (TGraphErrors*)f1->Get(name); if (!g1H19) cerr << "Could not get " << name << endl;

    sprintf(name, "UFC/%s/G4%s_UFC_Ch.1", ArrStrVar[var].c_str(), ArrStrVar[var].c_str());
    TGraphErrors *g0UFC = (TGraphErrors*)f0->Get(name); if (!g0UFC) cerr << "Could not get " << name << endl;
    TGraphErrors *g1UFC = (TGraphErrors*)f1->Get(name); if (!g1UFC) cerr << "Could not get " << name << endl;

    g0H19->SetLineWidth(1);
    g0H19->SetLineColor(kBlack);
    g0H19->SetMarkerColor(kBlack);
    g0H19->SetFillColorAlpha(kBlack, 0.2);
    g1H19->SetLineWidth(1);
    g1H19->SetLineColor(kRed);
    g1H19->SetMarkerColor(kRed);
    g1H19->SetFillColorAlpha(kRed, 0.2);
    g0UFC->SetLineWidth(1);
    g0UFC->SetLineColor(kBlue);
    g0UFC->SetMarkerColor(kBlue);
    g0UFC->SetFillColorAlpha(kBlue, 0.2);
    g1UFC->SetLineWidth(1);
    g1UFC->SetLineColor(kMagenta);
    g1UFC->SetMarkerColor(kMagenta);
    g1UFC->SetFillColorAlpha(kMagenta, 0.2);


    string Title = "#it{n}ELBE FC Correction";
    TMultiGraph *pMG = new TMultiGraph(ArrStrVarName[var].c_str(), Title.c_str());
    pMG->Add(g0H19);
    pMG->Add(g1H19);
    pMG->Add(g0UFC);
    pMG->Add(g1UFC);


    TCanvas *cUFCH19 = new TCanvas("cUFCH19", "Canvas displaying the correction factors", 1, 1, 1915, 1120);
    gPad->SetLogx();
    pMG->Draw("AL3");
    pMG->GetXaxis()->SetTitle("#it{E(t)} / MeV");
    pMG->GetXaxis()->SetRangeUser(0.1, 10.);

    Double_t yRange[3][2] = {{0.5, 1.0}, {0.75, 1.05}, {0.7, 1.4}};
    pMG->GetYaxis()->SetTitle(ArrStrVarName[var].c_str());
//    pMG->GetYaxis()->SetTitle(ArrStrVarName[var].c_str());
    pMG->GetYaxis()->SetRangeUser(yRange[var][0], yRange[var][1]);
    pMG->GetYaxis()->SetNdivisions(505);

    TLegend *l = new TLegend(0.25, 0.18, 0.85, 0.4);
    l->SetNColumns(3);
    l->SetLineColor(kWhite);
    l->AddEntry((TObject*)0, descr0.c_str(), "");
    l->AddEntry(g0H19, "H19 channel 1");
    l->AddEntry(g0UFC, "UFC channel 1");
    l->AddEntry((TObject*)0, descr1.c_str(), "");
    l->AddEntry(g1H19, "H19 channel 1");
    l->AddEntry(g1UFC, "UFC channel 1");
    l->Draw();

    cUFCH19->Update();

    string Filename = "pic/" + ArrStrSymbol[var] + "_" + descrFile + "_UFCH19.pdf";
    cUFCH19->SaveAs(Filename.c_str());
}

void CompAllCF(Int_t var, string StrFC, Int_t ch = 0)
{ /// Draw TGraphErrors for UFCvsH19 for different simulation versions, with and without track length.
  ///    StrVar: "TransMis", "Scattering", "Correction"
  ///    StrFC: "UFC", "H19", "RefH19"
    string StrPath = "~/TrackLength/";
    string ArrStrFile[] = {"G4UFCvsH19_23082016_hist.root",
                           "G4UFCvsH19_v10.2_hist.root",
                           "G4UFCvsH19_v10.2_hist_w.root",
                           "G4UFCvsH19_1E9_hist.root",
                           "G4UFCvsH19_1E9_hist_w.root",
                           "MCNP_0313_filled.root"};
    string ArrStrName[] = {"G4 10.2, 2016",
                         "G4 10.2, no TL",
                         "G4 10.2, with TL",
                         "G4 10.5, no TL",
                         "G4 10.5, with TL",
                         "MCNP6 1.1b"};
    string ArrStrVar[] = {"TransMis", "Scattering", "Correction"};
    string ArrStrVarName[] = {"Transmission factor  #it{T}",
                            "Correlation factor  #it{k}",
                            "Correction factor  #it{C}"};
    string ArrStrSymbol[] = {"T", "k", "C"};
//    uint Select[] = {0, 1, 3}; // version
//    uint N = 3;
    uint Select[] = {3, 4}; // TL
    uint N = 2;
//    uint Select[] = {1, 2, 4}; // version & TL
//    uint N = 3;

    Int_t Colors[] = {600, 867, 418, 632, 887, 820, 616, 807, 1, 921};
    TColor *Col_Pantone301C = new TColor (5000, 0,     0.345, 0.612, "Pantone301C");
    TColor *Col_Pantone717C = new TColor (5001, 0.843, 0.459, 0.251, "Pantone717C");

    //define TLegend
    TLegend *l1 = new TLegend(0.2, 0.15 , 0.9, 0.4);
    l1->SetTextSize(0.04);
    l1->SetNColumns(N);
    l1->SetFillStyle(0);
    l1->SetLineWidth(0);
    l1->SetLineColor(kWhite);

    string Title = "#it{n}ELBE FC " + ArrStrVarName[var];
    TMultiGraph *pMG = new TMultiGraph(ArrStrVar[var].c_str(), Title.c_str());

    string Filename, ObjName;
    TGraphErrors *pGraph;
    for (uint j = 0; j < N; j++)
    {
        uint i = Select[j];
        Filename = StrPath + ArrStrFile[i];
        TFile *f = TFile::Open(Filename.c_str(), "READ"); if (!f) cerr << "Could not open " << Filename << endl;
        ObjName = StrFC + "/" + ArrStrVar[var] + "/G4" + ArrStrVar[var] + "_" + StrFC + "_Ch." + std::to_string(ch+1);
        pGraph = (TGraphErrors*) f->Get(ObjName.c_str()); if (!pGraph) cerr << "Could not get " << ObjName << " from " << f->GetName() << endl;
        pGraph->SetLineWidth(1);
        pGraph->SetLineColor(Colors[i]);
        pGraph->SetMarkerColor(Colors[i]);
        pGraph->SetFillColorAlpha(Colors[i],0.2);
        pMG->Add(pGraph);//, "L3");
        l1->AddEntry(pGraph, ArrStrName[i].c_str(), "lf");

    }
    TCanvas *c = new TCanvas("TransMis", "Canvas displaying the correction factors", 1, 1, 1915, 1120);
    c->cd();
    c->cd(1)->SetTicks(1,1);
    pMG->Draw("AL3");
//    pMG->Draw("AL");
    pMG->GetXaxis()->SetTitle("#it{E(t)} / MeV");
    pMG->GetXaxis()->SetRangeUser(0.1, 10.);

    Double_t yRange[3][2] = {{0.5, 1.0}, {0.75, 1.05}, {0.7, 1.4}};
    pMG->GetYaxis()->SetTitle(ArrStrVarName[var].c_str());
    pMG->GetYaxis()->SetRangeUser(yRange[var][0], yRange[var][1]);
    pMG->GetYaxis()->SetNdivisions(505);

    c->SetLogx();
    c->Modified();
    c->Update();

    l1->Draw();
    TLatex latex;
        latex.SetTextSize(0.07);
        latex.SetTextAlign(11); //align at top
        latex.DrawLatexNDC(.825,.9,"#font[32]{#color[5001]{n}}#font[62]{#color[5000]{ELBE}}");
    TPave *pPave = new TPave(0.945, 0.9, 0.95, 0.91, 4, "NDC");
        pPave->SetFillColor(5001);
        pPave->Draw("NDC");

    c->Modified();
    c->Update();
//    Filename = "pic/" + ArrStrSymbol[var] + "_v_" + StrFC + "_" + std::to_string(ch+1) + ".pdf";
    Filename = "pic/" + ArrStrSymbol[var] + "_TL_" + StrFC + "_" + std::to_string(ch+1) + ".pdf";
//    Filename = "pic/" + ArrStrSymbol[var] + "_" + StrFC + "_" + std::to_string(ch+1) + ".pdf";
    c->SaveAs(Filename.c_str());
}

void DrawProjections(TFile *f, string StrFC, Int_t ch)
{ /// Draw E and E(t) projection of direct neutrons.
    string HistName = StrFC + "/EToFvsEkin/" + StrFC + "_EToFvsEkin_Ch." + std::to_string(ch+1);
    TH2F* h2Dir = (TH2F*) f->Get(HistName.c_str()); if (!h2Dir) cerr << "Could not get " << HistName << endl;
    HistName = StrFC + "/EToFvsEkin/Scattered/" + StrFC + "_EToFvsEkin_Sc_Ch." + std::to_string(ch+1);
    TH2F* h2Sc = (TH2F*) f->Get(HistName.c_str()); if (!h2Sc) cerr << "Could not get " << HistName << endl;
    h2Dir->Add(h2Sc, -1.0);
    TH1F* pEt = (TH1F*)h2Dir->ProjectionX();
    TH1F* pE = (TH1F*)h2Dir->ProjectionY();
    new TCanvas();
    pE->Draw("hist");
    pEt->SetLineColor(kRed);
    pEt->Draw("same hist");
    TLegend *l = new TLegend(0.3, 0.6, 0.5, 0.8);
    l->AddEntry(pE, "#font[12]{#Phi(E})");
    l->AddEntry(pEt, "#font[12]{#Phi(E(t))}");
    l->Draw();
}

void DrawCfProjections(TFile *f, string StrFC, Int_t ch)
{ /// Compare 3 ways of handling direct neutrons when calculating C(E(t)).
    // Get Direct 2D hist
    string HistName = StrFC + "/EToFvsEkin/" + StrFC + "_EToFvsEkin_Ch." + std::to_string(ch+1);
    TH2F* h2Dir = (TH2F*) f->Get(HistName.c_str()); if (!h2Dir) cerr << "Could not get " << HistName << endl;
    HistName = StrFC + "/EToFvsEkin/Scattered/" + StrFC + "_EToFvsEkin_Sc_Ch." + std::to_string(ch+1);
    TH2F* h2Sc = (TH2F*) f->Get(HistName.c_str()); if (!h2Sc) cerr << "Could not get " << HistName << endl;
    h2Dir->Add(h2Sc, -1.0);
    // Get cross section
    TGraphErrors *pGrXS = GetXS("~/TrackLength/FissionXS.root", StrFC, "JEFF_3.3");

    // way 1: cs(E) - pEt
    TH2F* h2FCDir = WeightWithXS(h2Dir, pGrXS);
    TH1F* pXS = (TH1F*)h2FCDir->ProjectionX();
    // way 2: pEt - cs(E(t))
    TH1F* pEt = (TH1F*)h2Dir->ProjectionX();
    TH1F* hEt = WeightTH1withXS(pEt, pGrXS);
    // way 3: pE - cs(E)
    TH1F* pE = (TH1F*)h2Dir->ProjectionY();
    TH1F* hE = WeightTH1withXS(pE, pGrXS);

    TCanvas *cP = new TCanvas();
//    pXS->SetLineColor(kBlue);
//    pXS->Draw("hist");
    pEt->SetLineColor(kRed);
    pEt->Draw("same hist");
    pE->SetLineColor(kOrange);
    pE->Draw("same hist");
    TLegend *l = new TLegend(0.3, 0.6, 0.5, 0.8);
//    l->AddEntry(pXS, "#sigma(E) - px");
    l->AddEntry(pEt, "px - #sigma(E(t))");
    l->AddEntry(pE, "py - #sigma(E)");
    l->Draw();
    cout << hEt->Integral() << " " << hE->Integral() << " " << pXS->Integral() << endl;
    cP->Update();
}

void drawTL()
{
    LoadStyles();
    gROOT->SetStyle("SinglePadStyle");
    gROOT->ForceStyle(kTRUE);
    gStyle->SetLegendFont(132);

    /// simulation files
//    TFile *fs1 = TFile::Open("/home/hoffma93/bigdata/nELBE/G4TransMis/20200128/results/G4UFCvsH19_28012020.root", "READ");
//    TFile *fs2 = TFile::Open("/home/hoffma93/bigdata/nELBE/G4TransMis/20200131/G4UFCvsH19_31012020.root", "READ");
//    TFile *fs3 = TFile::Open("/home/hoffma93/bigdata/nELBE/G4TransMis/20200211/results/G4TransMis_11022020.root", "READ");
//    TFile *fs4 = TFile::Open("/home/hoffma93/bigdata/nELBE/G4TransMis/20200217/results/G4UFCvsH19_17022020.root", "READ");
//    TFile *fv = TFile::Open("~/TrackLength/G4UFCvsH19_vac_1E7_hist.root");
    /// histogram files
    TFile *fh0 = TFile::Open("~/TrackLength/G4UFCvsH19_23082016_hist.root", "READ");
//    TFile *fh1 = TFile::Open("~/TrackLength/G4UFCvsH19_v10.2_hist.root", "READ");
//    TFile *fh2 = TFile::Open("~/TrackLength/G4UFCvsH19_v10.2_hist_w.root", "READ");
    TFile *fh3 = TFile::Open("~/TrackLength/G4UFCvsH19_1E9_hist.root", "READ");
//    TFile *fh4 = TFile::Open("~/TrackLength/G4UFCvsH19_1E9_hist_w.root", "READ");
//    TFile *fh5 = TFile::Open("~/TrackLength/MCNP_0313.root", "READ");

//    DrawWeight(fh2, "kin");
//    DrawWspectrum(StrFileInPath, "UFC", 1);
//    DrawWspectra(StrFileInPath);
//    DrawRefH19(File3, 7);

//    DrawWspectra(fs4, fs3, "UFC", 0, 1);

//    DrawFissionEE(f, "UFC", 0, 2, 0);
//    DrawFissionEt(f, "UFC", 0, 2, 0);
//    DrawProjections(fv, "UFC", 7);
//    DrawCfProjections(f, "UFC", 0);
//    CompSimRes(fh3, "2020", fh7, "2016", "UFC", 2, "dual", "Ekin", 1);
//    CompProjections("UFC", 0, 0);
//    CompCF("H19", 7, fh3, "G4 10.5 no TL", fh4, "G4 10.5 with TL");
    CompCF_UFC_H19(0, fh0, "G4 10.2, 2016:", fh3, "G4 10.5, no TL:", "v");
//    CompAllCF(2, "H19", 0);
}
