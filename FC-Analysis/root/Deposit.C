#ifndef DEPOSIT_H
#define DEPOSIT_H
#include "fstream"
#include "TH2.h"
#include "TAxis.h"
#include "TEllipse.h"
#include "TArc.h"
#define GridX 25
#define GridY 25

TH2D* ReadTxtMatrix(string isotope, Int_t ch)
{
    string file = isotope[0] == 'U' ?
                "/home/hoffma93/Experiment/Autoradiographien/U235/Result_U#"+to_string(ch+1)+".txt" :
                "/home/hoffma93/Experiment/Autoradiographien/Pu242/Result_"+to_string(ch+1)+".txt";
    std::ifstream ifs (file.c_str(), std::ifstream::in);
    if (!ifs.is_open())
        return 0;
    TH2D *pH2 = new TH2D("name", "title", GridX+2, -80/25., 80/25.*26, GridY+2, -80/25., 80/25.*26); // GridX, 0, 80, GridY, 0, 80);
    pH2->SetName((isotope+"_"+to_string(ch+1)).c_str());
    string line;
    std::getline(ifs, line);
    Int_t x = 0, y;
     while (std::getline(ifs, line))
     {
         std::stringstream s(line);
         string dummy; s >> dummy;
         Double_t z;
         y = 0;
         while (s >> z) {
             pH2->SetBinContent(x+2, y+2, z);
             y++; }
         x++;
     }
     ifs.close();
     pH2->SetBinContent(25, 25, 0.5*(pH2->GetBinContent(24, 25)+pH2->GetBinContent(25, 24)));
     Double_t val = -1;
     for (Int_t b = 1; b <= GridX+2; b++)
     {
         pH2->SetBinContent(b, 1, val);
         pH2->SetBinContent(GridX+2, b, val);
         pH2->SetBinContent(b, GridY+2, val);
         pH2->SetBinContent(1, b, val);
     }
     pH2->Scale(1.0 / pH2->GetMaximum());
//     pH2->Scale(1.0/pH2->GetMaximum());
     return pH2;
}

Double_t AnalyzeDeposit(TH2D *h, Double_t beam_radius = 28, Double_t outer_radius = 38)
{
    Double_t xCenter = 40, yCenter = 40;
    Double_t sDep = 0, sBeam = 0, sBg = 0;
    Int_t nDep = 0, nBeam = 0, nBg = 0;
    Double_t c, x, y, r;
    for (Int_t xBin = 1; xBin < h->GetXaxis()->GetNbins()+1; xBin++)
    {   for (Int_t yBin = 1; yBin < h->GetYaxis()->GetNbins()+1; yBin++)
        {
            c = h->GetBinContent(xBin, yBin);
            if (c == 0) // Exclude entries of zero
                continue;
            x = h->GetXaxis()->GetBinCenter(xBin) - xCenter;
            y = h->GetYaxis()->GetBinCenter(yBin) - yCenter;
            r = sqrt(x*x + y*y);
            if (r < outer_radius)
            {
                nDep += 1;
                sDep += c;
            } else {
                nBg += 1;
                sBg += c;
            }
            if (r < beam_radius)
            {
                nBeam += 1;
                sBeam += c;
            }
        }
    }

    Double_t Fg = sDep - nDep * sBg / nBg;
    Double_t Beam = sBeam - nBeam * sBg / nBg;

//    cout << "sDep " << sDep << endl <<
//            "nDep " << nDep << endl <<
//            "sBg " << sBg << endl <<
//            "nBg " << nBg << endl <<
//            "sBeam " << sBeam << endl <<
//            "nBeam " << nBeam << endl <<
//            "Fg " << Fg << endl;

    return TMath::Pi() * pow(beam_radius, 2) * Fg / Beam;
}

void Deposit(Int_t ch, string isotope = "U235")
{
    TH2D *srfc = ReadTxtMatrix(isotope, ch);
    cout << ch+1 << " " << AnalyzeDeposit(srfc, 27.16, 45) << endl;
}

void DrawDeposit(Int_t ch, string isotope = "U235")
{
    TH2D *srfc = ReadTxtMatrix(isotope, ch);
    gStyle->SetPalette(kRainBow);
//    TColor::InvertPalette();
//    srfc->SetMinimum(-0.0);
    srfc->GetZaxis()->SetRangeUser(0, 1);
    srfc->GetXaxis()->SetRangeUser(-0.1, 80.1);
//    srfc->GetYaxis()->SetRangeUser(0, 80);
    srfc->SetStats(0);
    srfc->SetTitle("; #font[12]{x} / mm; #font[12]{y} / mm");
    Int_t w = 1000, h = 1000;
    TCanvas *c1 = new TCanvas("c1", "c1", w, h);
    gPad->SetTicks(1,1);
    c1->SetWindowSize(w + (w - c1->GetWw()), h + (h - c1->GetWh()));
//    srfc->Draw("SURF1Z FB BB");
    srfc->Draw("CONT4Z");
    srfc->Draw("CONT0Z");
    TEllipse *e = new TEllipse(40, 40, 28, 28);
    e->SetLineColor(kGray);
    e->SetFillColorAlpha(0, 0);
    e->Draw();
    TArc *a1 = new TArc(40, 40, 45, 27, 63);
    a1->SetLineColor(kGray);
    a1->SetNoEdges();
    a1->SetFillColorAlpha(0, 0);
    a1->Draw();
    TArc *a2 = new TArc(40, 40, 45, 117, 153);
    a2->SetLineColor(kGray);
    a2->SetNoEdges();
    a2->SetFillColorAlpha(0, 0);
    a2->Draw();
    TArc *a3 = new TArc(40, 40, 45, 207, 243);
    a3->SetLineColor(kGray);
    a3->SetNoEdges();
    a3->SetFillColorAlpha(0, 0);
    a3->Draw();
    TArc *a4 = new TArc(40, 40, 45, 297, 333);
    a4->SetLineColor(kGray);
    a4->SetNoEdges();
    a4->SetFillColorAlpha(0, 0);
    a4->Draw("same");
    TBox *b = new TBox(-1, -1, 81, 81);
    b->SetFillStyle(0);
    b->SetLineColor(kWhite);
    b->SetLineWidth(12);
    b->Draw();
    char name[64] = "";
    sprintf(name, "~/Pictures/nU/%s_Dep_%i.pdf", isotope.c_str(), ch+1);
//    c1->SaveAs(name);
}

void DrawDeposits(string isotope = "U235")
{
    gStyle->SetPalette(kRainBow);
//    TColor::InvertPalette();
    Int_t w = 1100, h = 1100;
    TCanvas *c1 = new TCanvas("c1", "c1", w, h);
    c1->SetWindowSize(w + (w - c1->GetWw()), h + (h - c1->GetWh()));
    c1->Divide(3, 3);
    Int_t nr[] = {0, 2, 3, 4, 5, 6, 7, 8};
    Double_t RightMargin[] = {0.01, 0.1, 0.2, 0.01, 0.1, 0.2, 0.01, 0.1, 0.2};
    Double_t TopMargin[] = {0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.01, 0.01, 0.01};
    Double_t LeftMargin[] = {0.2, 0.1, 0.01, 0.2, 0.1, 0.01, 0.2, 0.1, 0};
    Double_t BottomMargin[] = {0.01, 0.01, 0.01, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2};
    TH2D *srfc[9];
    for (Int_t i = 0; i < 8; i++)
    {
        srfc[i] = ReadTxtMatrix(isotope, nr[i]);
        //    srfc[i]->SetMinimum(-0.0);
        srfc[i]->GetZaxis()->SetRangeUser(0, 1);
        srfc[i]->GetXaxis()->SetRangeUser(-0.1, 80.1);
        //    srfc[i]->GetYaxis()->SetRangeUser(0, 80);
        srfc[i]->SetStats(0);
        srfc[i]->SetTitle(";;");
        srfc[i]->GetXaxis()->SetLabelSize(0.0);
        srfc[i]->GetXaxis()->SetNdivisions(405);
        srfc[i]->GetYaxis()->SetLabelSize(0.0);
        srfc[i]->GetYaxis()->SetNdivisions(405);
        if (i == 3) {
            srfc[i]->SetTitle(";; #font[12]{y} / mm");
            srfc[i]->GetYaxis()->SetLabelSize(0.07);
            srfc[i]->GetYaxis()->CenterTitle();
            srfc[i]->GetYaxis()->SetTitleSize(0.07); }
        if (i == 7) {
            srfc[i]->SetTitle("; #font[12]{x} / mm;");
            srfc[i]->GetXaxis()->SetLabelSize(0.07);
            srfc[i]->GetXaxis()->CenterTitle();
            srfc[i]->GetXaxis()->SetTitleSize(0.07); }
        if (i == 0 || i == 6)
            srfc[i]->GetYaxis()->SetLabelSize(0.07);
        if (i == 6 || i == 8)
            srfc[i]->GetXaxis()->SetLabelSize(0.07);


        c1->cd(i+1);
        gPad->SetTicks(1,1);
        gPad->SetRightMargin(RightMargin[i]);
        gPad->SetTopMargin(TopMargin[i]);
        gPad->SetLeftMargin(LeftMargin[i]);
        gPad->SetBottomMargin(BottomMargin[i]);

//        srfc[i]->Draw("SURF1Z FB BB");
        srfc[i]->Draw("CONT4");
        srfc[i]->Draw("CONT0");
        TBox *b = new TBox(-1, -1, 81, 81);
        b->SetFillStyle(0);
        b->SetLineColor(kWhite);
        b->SetLineWidth(4);
        b->Draw();
    }
    srfc[8] = (TH2D*)srfc[0]->Clone((isotope+"_tot").c_str());
    srfc[8]->Add(srfc[1]);
    srfc[8]->Add(srfc[2]);
    srfc[8]->Add(srfc[3]);
    srfc[8]->Add(srfc[4]);
    srfc[8]->Add(srfc[5]);
    srfc[8]->Add(srfc[6]);
    srfc[8]->Add(srfc[7]);
    srfc[8]->Scale(0.125);//srfc[8]->GetMaximum());
    srfc[8]->GetZaxis()->SetRangeUser(0, 1);
    srfc[8]->GetXaxis()->SetRangeUser(-0.1, 80.1);
    srfc[8]->SetStats(0);
    srfc[8]->SetTitle(";;");
    srfc[8]->GetXaxis()->SetLabelSize(0.07);
    srfc[8]->GetXaxis()->SetNdivisions(405);
    srfc[8]->GetYaxis()->SetLabelSize(0.0);
    srfc[8]->GetYaxis()->SetNdivisions(405);
    srfc[8]->GetZaxis()->SetLabelSize(0.07);
    c1->cd(9);
    gPad->SetTicks(1,1);
    gPad->SetRightMargin(RightMargin[8]);
    gPad->SetTopMargin(TopMargin[8]);
    gPad->SetLeftMargin(LeftMargin[8]);
    gPad->SetBottomMargin(BottomMargin[8]);
    srfc[8]->Draw("CONT4Z");
    srfc[8]->Draw("CONT0Z");
    TBox *b = new TBox(-1, -1, 81, 81);
    b->SetFillStyle(0);
    b->SetLineColor(kWhite);
    b->SetLineWidth(4);
    b->Draw();

    TEllipse *e = new TEllipse(40, 40, 28, 28);
    e->SetLineColor(kGray);
    e->SetFillColorAlpha(0, 0);
    e->Draw();
    TArc *a1 = new TArc(40, 40, 45, 27, 63);
    a1->SetLineColor(kGray);
    a1->SetNoEdges();
    a1->SetFillColorAlpha(0, 0);
    a1->Draw();
    TArc *a2 = new TArc(40, 40, 45, 117, 153);
    a2->SetLineColor(kGray);
    a2->SetNoEdges();
    a2->SetFillColorAlpha(0, 0);
    a2->Draw();
    TArc *a3 = new TArc(40, 40, 45, 207, 243);
    a3->SetLineColor(kGray);
    a3->SetNoEdges();
    a3->SetFillColorAlpha(0, 0);
    a3->Draw();
    TArc *a4 = new TArc(40, 40, 45, 297, 333);
    a4->SetLineColor(kGray);
    a4->SetNoEdges();
    a4->SetFillColorAlpha(0, 0);
    a4->Draw("same");

    // Save
//    char name[64] = "";
//    sprintf(name, "~/Pictures/nU/%s_Dep_%i.pdf", isotope.c_str(), ch+1);
//    c1->SaveAs(name);
}

void Deposit()
{
//    for (Int_t i = 0; i < 8; i++)
//        Deposit(i, "U235");
//    DrawDeposit(0);
    DrawDeposits();
}
#endif
