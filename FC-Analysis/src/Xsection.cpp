#include "Xsection.h"

using namespace std;

Xsection::Xsection()
{
    CommentFlag = kFALSE;

    // physics parameters
    u = 1.660539E-24; // atomic mass unit in g
    PuLit = 2120E-25; // 242Pu neutron-induced fission cross section in mm^2, original in mb
    DPuLit = 35E-25;
    ULit = 2.1E-22; // 235U nif cs in mm^2, original in b
    DULit = 0.0315E-22;
    PuSFT2 = 6.77E10 * 365.24*24*60*60; // Pu-242 spontaneaus fission half-life period in s
    DPuSFT2 = 7E8 * 365.24*24*60*60;
    Yield = 2.15882E4;
    DYield = 543.3;
    MonitorNIF = 27492079+33478370+30916955+54792679;
    DMonitorNIF = MonitorNIF * 0.0014;
    MonitorSB = 3623069+28646614+4757385;
    DMonitorSB = MonitorSB * 0.0015;
    MonitorUNIF = 9509138+12371470+12876188+27941726+9962706;
    DMonitorUNIF = MonitorUNIF * 0.0014;
    MonitorUSB = 13815794;
    DMonitorUSB = MonitorUSB * 0.0014;
    AreaPuFC = 4300; // in mm^2
    DAreaPuFC = 120;
    mPu = 0.03724; // in g
    DmPu = 0.00001;
    mU = 0.15801;
    DmU = 0.00001;

    pHNIF = new Hist("/home/hoffma93/Go4nfis/offline/results/NIF.root", "NIF");
    pHNIF->SetNeutronField(Yield, DYield, MonitorNIF, DMonitorNIF, 1500, 1);
    pHSB  = new Hist("/home/hoffma93/Go4nfis/offline/results/SB.root", "SB");
    pHSB->SetNeutronField(Yield, DYield, MonitorSB, DMonitorSB, 1500, 1);
    pHSF  = new Hist("/home/hoffma93/Go4nfis/offline/results/SF.root", "SF");
    pHNIF->DoAnalyzeDt();
    pHNIF->DoAnalyzeQDC();
    pHSB ->DoAnalyzeDt();
    pHSB ->DoAnalyzeQDC();
    pHSF ->DoAnalyzeDt();
    pHSF ->DoAnalyzeQDC();

//    pHUNIF = new Hist("/home/hoffma93/Go4nfis/offline/results/UFC_NIF.root", "NIF", "UFC");
//    pHUNIF->SetNeutronField(Yield, DYield, MonitorUNIF, DMonitorUNIF, 1500, 1);
//    pHUSB = new Hist("/home/hoffma93/Go4nfis/offline/results/UFC_SB.root", "SB", "UFC");
//    pHUSB->SetNeutronField(Yield, DYield, MonitorUSB, DMonitorUSB, 1500, 1);
//    pHUNIF->DoAnalyzeDt();
//    pHUNIF->DoAnalyzeQDC();
//    pHUSB->DoAnalyzeDt();
//    pHUSB->DoAnalyzeQDC();

    CalculateThresholds();
    PrintInScat();
    ScatCorrNIF();
    CalculateNPu();
//    Evaluation1();
//    Evaluation2();
//    Evaluation3();

//    CompareFiles("/home/hoffma93/Go4nfis/offline/results/test", 191, 204);
}
void Xsection::CalculateThresholds()
{ // Merge QDC spectra analyses. Calculate average extrema positions
//    cout << "Common cut positions" << endl;
    Double_t PuFC_Cut[NumHist];
    Double_t PuFC_DCut[NumHist];
    Double_t UFC_Cut[NumHist];
    Double_t UFC_DCut[NumHist];
    Double_t NIF_Ped, NIF_Cut, NIF_Max,
             SB_Ped, SB_Cut, SB_Max,
             SF_Ped, SF_Cut, SF_Max;
    /// method 1: Allow different relative minimum positions in deposits
    for (int i = 0; i < NumHist; i++)
    {
        // PuFC
        NIF_Cut = pHNIF->CutQDC[i];
        SB_Cut = pHSB->CutQDC[i];
        SF_Cut = pHSF->CutQDC[i];
        Double_t counts[] = {pHNIF->GetNumberEvents(i), pHSB->GetNumberEvents(i), pHSF->GetNumberEvents(i)};
//        {pHNIF->t_live*(pHNIF->SFRate[i]+pHNIF->NIFRate[i]), // number of FF events for relative weighting
//                             pHSB->t_live*(pHSB->NIFRate[i]+pHSB->NIFRate[i]),
//                             pHSF->t_live*(pHSF->SFRate[i]+pHSF->NIFRate[i])};//
        Double_t sum = counts[0] + counts[1] + counts[2];
        Double_t w[] = {counts[0]/sum, counts[1]/sum, counts[2]/sum}; // weights
//        cout << "ch " << i+1 << "  " << w[0] << "  " << w[1] << "  " << w[2] << endl;
        PuFC_Cut[i] = w[0] * NIF_Cut + w[1] * SB_Cut + w[2] * SF_Cut;
        PuFC_DCut[i] = sqrt( w[0] * pow(NIF_Cut - PuFC_Cut[i], 2) +
                             w[1] * pow(SB_Cut - PuFC_Cut[i], 2) +
                             w[2] * pow(SF_Cut - PuFC_Cut[i], 2) );
//        cout << "      " << PuFC_Cut[i] << "+-" << PuFC_DCut[i] << endl;
        // UFC...
    }
    // Draw
    double x[] = {1, 2, 3, 4, 5, 6, 7, 8};
    double xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};

    TGraphErrors* g1 = new TGraphErrors(NumHist, x, PuFC_Cut, xerr, PuFC_DCut);
    g1->SetNameTitle("gCut", "PuFC combined QDC minimum position vs channel");
    SaveToFile("Analysis/QDC", g1);
}

void Xsection::CompareFiles(string path, Int_t start, Int_t stop)
{
    TCanvas* c1 = new TCanvas("SFRates", "SFRates", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg1 = new TMultiGraph();
    mg1->SetTitle("Spontaneaus fission detection rates; Deposit; Rate / Hz");
    TLegend* l1 = new TLegend(0.6, 0.2, 0.85, 0.40, "File nr");

    char fname[64] = "";
    Double_t x[] = {1,2,3,4,5,6,7,8};
    Double_t xerr[] = {0,0,0,0,0,0,0,0};
    Color_t color[] = {kRed, kGreen, kBlue, kYellow, kBlue, kYellow, kOrange, kRed, kGreen, kBlue, kGray, kYellow, kMagenta, kOrange};

    for ( int i = start; i < stop; i++ )
    {
        sprintf(fname, "%s/SB%i.root", path.c_str(), i);
        Hist* pH = new Hist(fname, "SB");
        pH->DoAnalyzeDt();
        cout << "nr " << i << ",  t_live/t_real=" << pH->t_live / pH->t_real << endl;
        for (int j = 0; j < 8; j++)
        {
            cout << "  ch " << j+1 << ",  SF rate=" << pH->SFRate[j] << endl;
        }
        TGraphErrors* ge = new TGraphErrors(NumHist, x, pH->SFRate, xerr, pH->DSFRate);
//        sprintf(fname, "SF_Rate_%i", i);
//        ge->SetNameTitle(fname, fname);
        ge->SetLineWidth(2);
        ge->SetLineColor(color[i-start]);
        mg1->Add(ge);
        sprintf(fname, "%i", i);
        l1->AddEntry(ge, fname, "lp");
    }
    mg1->Draw("AP");
    l1->Draw();
//    mg1->GetYaxis()->SetRangeUser(0, 2);
    c1->Modified();
    c1->Update();
}

Xsection::~Xsection()
{

}

void Xsection::PrintInScat()
{
    cout << "=== In-scattering corrections ===" << endl;
    cout << "Neutron-induced fissions per monitor count" << endl;
    cout << "Channel nr, Free, Shadow bar, Ratio " << endl;
    for(int i = 0; i < 8; i++)
    {
        cout << /*"  Channel " <<*/ i+1 /*<< endl*/;
        cout << ",  " << pHNIF->NIFRate[i] / MonitorNIF * pHNIF->t_real/*
             << " +- " << pHNIF->DNIFRate[i] / MonitorNIF * pHNIF->t_real << endl*/;
        cout << ",  " << pHSB->NIFRate[i] / MonitorSB * pHSB->t_real/*
             << " +- " << pHSB->DNIFRate[i] / MonitorSB * pHSB->t_real << endl*/;
//        cout << ",  " << pHNIF->NIFRate[i] / MonitorNIF * pHNIF->t_real - pHSB->NIFRate[i] / MonitorSB * pHSB->t_real
//             << " +- " << sqrt( pow(pHNIF->DNIFRate[i] / MonitorNIF * pHNIF->t_real, 2) +
//                                pow(pHSB->DNIFRate[i] / MonitorSB * pHSB->t_real, 2) ) /*<< endl*/;
        cout << ",  " << pHSB->NIFRate[i] / MonitorNIF * pHNIF->t_real / pHNIF->NIFRate[i] * MonitorSB / pHSB->t_real
             << " +- " << sqrt( pow(pHSB->DNIFRate[i]/pHNIF->NIFRate[i], 2) +
                                pow(pHNIF->DNIFRate[i]*pHSB->NIFRate[i]/pow(pHNIF->NIFRate[i], 2), 2) ) << endl;
    }
}

void Xsection::ScatCorrNIF()
{ // Calculate In-scattering corrected NIF rate normed to neutron flux
    for (int i = 0; i < NumHist; i++)
    {
        NIFRate[i] = pHNIF->NIFRate[i] - pHSB->NIFRate[i] * pHNIF->NeutronFlux[i] / pHSB->NeutronFlux[i];
        DNIFRate[i] = sqrt( pow(pHNIF->DNIFRate[i], 2) +
                            pow(pHSB->DNIFRate[i] * pHNIF->NeutronFlux[i] / pHSB->NeutronFlux[i], 2) +
                            pow(pHSB->NIFRate[i] * pHNIF->DNeutronFlux[i] / pHSB->NeutronFlux[i], 2) +
                            pow(pHSB->NIFRate[i] * pHNIF->NeutronFlux[i] * pHSB->DNeutronFlux[i] / pow(pHSB->NeutronFlux[i], 2), 2) );
    }
}

void Xsection::CalculateNPu()
{
    Double_t nPuEff[NumHist];
    Double_t DnPuEff[NumHist];
    for(int i_ch = 0; i_ch < NumHist; i_ch++)
    {   // Calculate SF rate in 1/s and N(242Pu) for each plate.
        cout << "Spontaneaus fission rates in s^-1 for channel " << i_ch + 1 << endl;
        // NIF measurement:
        cout << "    NIF: " << pHNIF->SFRate[i_ch] << " +- " << pHNIF->DSFRate[i_ch] << endl;
        cout << "    SB:  " << pHSB->SFRate[i_ch] << " +- " << pHSB->DSFRate[i_ch] << endl;
        cout << "    SF:  " << pHSF->SFRate[i_ch] << " +- " << pHSF->DSFRate[i_ch] << endl;
        // Combine 3 measurements
        Double_t sum = 1/pHNIF->DSFRate[i_ch] + 1/pHSB->DSFRate[i_ch] + 1/pHSF->DSFRate[i_ch];
        Double_t w[] = {1/pHNIF->DSFRate[i_ch]/sum, // weighting rates according to inverse uncertainty
                        1/pHSB->DSFRate[i_ch]/sum,
                        1/pHSF->DSFRate[i_ch]/sum};
        SFRate[i_ch] = w[0] * pHNIF->SFRate[i_ch] +
                 w[1] * pHSB->SFRate[i_ch] +
                 w[2] * pHSF->SFRate[i_ch];
        DSFRate[i_ch] = sqrt( pow(w[0] * pHNIF->DSFRate[i_ch], 2) +
                              pow(w[1] * pHSB->DSFRate[i_ch], 2) +
                              pow(w[2] * pHSF->DSFRate[i_ch], 2) );
        cout << "    All: " << SFRate[i_ch] << " +- " << DSFRate[i_ch] << endl;

        // Calculate effective N(242Pu)
        nPuEff[i_ch] = SFRate[i_ch] * PuSFT2 / log(2.0); // SFRate, PuSFT2: both in seconds!
        DnPuEff[i_ch] = sqrt( pow(DSFRate[i_ch] * PuSFT2 / log(2.0) / pHSF->eInt[i_ch], 2) +
                           pow(SFRate[i_ch] * DPuSFT2 / log(2.0) / pHSF->eInt[i_ch], 2) );
        // Calculate N(242Pu)
        nPu[i_ch] = nPuEff[i_ch] / pHSF->eInt[i_ch];
        DnPu[i_ch] = sqrt( pow(DnPuEff[i_ch] / pHSF->eInt[i_ch], 2) +
                           pow(nPuEff[i_ch] * pHSF->DeInt[i_ch] / pow(pHSF->eInt[i_ch], 2), 2) );
        cout << "    #Pu: " << nPu[i_ch] << " +- " << DnPu[i_ch] << ", eff=" << pHSF->eInt[i_ch] << endl;
    }

    // Draw...
    double x[] = {1, 2, 3, 4, 5, 6, 7, 8};
    double xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};

    TGraphErrors* g1 = new TGraphErrors(NumHist, x, SFRate, xerr, DSFRate);
    g1->SetNameTitle("SF_Rate", "PuFC combined SF detection rate");
    SaveToFile("Analysis/SF", g1);

    TGraphErrors* g2 = new TGraphErrors(NumHist, x, nPuEff, xerr, DnPuEff);
    g2->SetNameTitle("NPuEff", "Effective number of 242Pu atoms from SF");
    SaveToFile("Analysis/SF", g2);

    TGraphErrors* g3 = new TGraphErrors(NumHist, x, nPu, xerr, DnPu);
    g3->SetNameTitle("NPu", "Number of 242Pu atoms from spontaneaus fission");
    SaveToFile("Analysis/SF", g3);
//*/
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////  Evaluation 1: #FF, neutron flux, #SF, T2_SF  -->  #Pu, sigma_nfis                         ////
////////////////////////////////////////////////////////////////////////////////////////////////////

void Xsection::Evaluation1()
{
    cout << endl << "=== Evaluation 1 ===" << endl;
//    CalculateNPu();
    CalculateCrossSection();
}


void Xsection::CalculateCrossSection()
{
    cout << "Cross section" << endl;
    for(int i_ch = 0; i_ch < NumHist; i_ch ++)
    {
        Double_t FluxNIF = pHNIF->NeutronFlux[i_ch];
        Double_t DFluxNIF = pHNIF->DNeutronFlux[i_ch];

        // cross section = reaction rate / (number of target atoms * flux)
        Double_t Xs = NIFRate[i_ch] / (nPu[i_ch] * FluxNIF);
        Double_t DXs = sqrt( pow(DNIFRate[i_ch] / (nPu[i_ch] * FluxNIF), 2) +
                               pow(NIFRate[i_ch] * DnPu[i_ch] / (pow(nPu[i_ch], 2) * FluxNIF), 2) +
                               pow(NIFRate[i_ch] * DFluxNIF / (nPu[i_ch] * pow(FluxNIF, 2)), 2) );

        cout << "XS(ch " << i_ch+1 << ") = " << Xs << " +- " << DXs << " mm^2" << endl;
        cout << "         = " << Xs * 1.E22 << " +- " << DXs * 1.E22 << " barn" << endl;
        XSec[i_ch] = Xs;
        DXSec[i_ch] = DXs;
    }
    // Calculate cross section average across channels
    Double_t sum = 0;
    Double_t Dsum = 0;
    for (int i = 0; i < NumHist; i++)
    {
        sum += XSec[i];
        Dsum += pow(DXSec[i], 2);
    }
    CrossSection = sum / NumHist;
    DCrossSection = sqrt(Dsum/NumHist);
    cout << "Combined:  " << CrossSection * 1.E22 << " barn +- " << 100 * DCrossSection / CrossSection << " %" << endl;

    Double_t x[] = {1, 2, 3, 4, 5, 6, 7, 8};
    Double_t xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};
    TGraphErrors* g3 = new TGraphErrors(NumHist, x, XSec, xerr, DXSec);
    g3->SetNameTitle("XS_SF", "Cross section (PuFC alone)");
    SaveToFile("Analysis/Evaluation", g3);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
////  Evaluation 2: #FF, neutron flux, sigma_nfis(Lit.), T2_SF  -->  #Pu                        ////
////////////////////////////////////////////////////////////////////////////////////////////////////
void Xsection::Evaluation2()
{
    cout << endl << "=== Evaluation 2 ===" << endl;
    Double_t nPu[NumHist];
    Double_t DnPu[NumHist];
    Double_t sum = 0;
    Double_t Dsum = 0;
    cout << "Effective number of Pu-242 atoms calculated with literature value" << endl <<
            "sigma_nfis = " << PuLit*1.E22 << " +- " << DPuLit*1.E22 << " barn" << endl;
    for (int i = 0; i < NumHist; i++)
    {
        nPu[i] = NIFRate[i] / (pHNIF->NeutronFlux[i] * PuLit);
        DnPu[i] = sqrt( pow(DNIFRate[i] / (pHNIF->NeutronFlux[i] * PuLit), 2) +
                        pow(NIFRate[i] * pHNIF->DNeutronFlux[i] / ( pow(pHNIF->NeutronFlux[i], 2) * PuLit ), 2) +
                        pow(NIFRate[i] * DPuLit / ( pHNIF->NeutronFlux[i] * pow(PuLit, 2) ), 2) );
        cout << "N(ch " << i+1 << ") = " << nPu[i] << /*" +- " << DnPu[i] << */endl;
        sum += nPu[i];
        Dsum += pow(DnPu[i], 2);
    }
    cout << "Sum = " << sum << " +- " << sqrt(Dsum) << endl;

    double x[] = {1, 2, 3, 4, 5, 6, 7, 8};
    double xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};

    TGraphErrors* g1 = new TGraphErrors(NumHist, x, nPu, xerr, DnPu);
    g1->SetNameTitle("NPuLit", "Effective number of 242Pu atoms from given cross section");
    SaveToFile("Analysis/Evaluation", g1);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
////  Evaluation 3: #FF, neutron flux, sigma_nfis(Lit.), #SF, T2_SF  -->  efficiency            ////
////////////////////////////////////////////////////////////////////////////////////////////////////
void Xsection::Evaluation3()
{
    cout << endl << "=== Evaluation 3 ===" << endl;
    cout << "Detection efficiencies" << endl;

    Double_t m = 242.059 * u; // mass of one 242Pu atom in g
    Double_t N = mPu / m; // total number of 242Pu atoms assuming given mass
    Double_t DN = DmPu / m;

    Double_t sumSF = 0;
    Double_t D2sumSF = 0;
    for (int i = 0; i < NumHist; i++)
    {
        sumSF += SFRate[i];
        D2sumSF += pow(DSFRate[i], 2);
    }

    eSF = sumSF * PuSFT2 / (N * log(2.0));
    DeSF = sqrt( D2sumSF*pow(PuSFT2 / (N * log(2.0)), 2) +
                         pow(sumSF * DPuSFT2 / (N * log(2.0)), 2) +
                         pow(sumSF * PuSFT2 * DN / (N * N * log(2.0)), 2) );
    cout << "average SF eff.: " << eSF << " +- " << DeSF << endl;
    cout << "ch  NIF  NIF/SF" << endl;
    for (int i = 0; i < NumHist; i++)
    {
        Double_t n = N * SFRate[i] / sumSF;
        Double_t Dn = sqrt( pow(DN * SFRate[i] / sumSF, 2) +
                            pow(N * DSFRate[i] / sumSF, 2) +
                  D2sumSF * pow(N * SFRate[i] / (sumSF*sumSF), 2) );
        eNIF[i] = NIFRate[i] / (n * PuLit * pHNIF->NeutronFlux[i]);
        DeNIF[i] = eNIF[i] * sqrt( pow(DNIFRate[i] / NIFRate[i], 2) +
                                   pow(Dn / n, 2) +
                                   pow(DPuLit / PuLit, 2) +
                                   pow(pHNIF->DNeutronFlux[i] / pHNIF->NeutronFlux[i], 2) );

        eRel[i] = NIFRate[i] / (nPu[i] * PuLit * pHNIF->NeutronFlux[i]);
        DeRel[i] = sqrt( pow(DNIFRate[i] / NIFRate[i], 2) +
                         pow(DnPu[i] / nPu[i], 2) +
                         pow(DPuLit / PuLit, 2) +
                         pow(pHNIF->DNeutronFlux[i] / pHNIF->NeutronFlux[i], 2)
                       ) * eRel[i];
//        eNIF[i] = eRel[i] * eSF;
//        DeNIF[i] = sqrt( pow(DeRel[i] * eSF, 2) +
//                         pow(eRel[i] * DeSF, 2) );
        cout << " " << i+1 << "  " << eNIF[i] << "+-" << DeNIF[i] << "  " << eRel[i] << "+-" << DeRel[i] << endl;
    }
    Double_t x[] = {1, 2, 3, 4, 5, 6, 7, 8};
    Double_t xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};
    Double_t y[] = {eSF, eSF, eSF, eSF, eSF, eSF, eSF, eSF};
    Double_t yerr[] = {DeSF, DeSF, DeSF, DeSF, DeSF, DeSF, DeSF, DeSF};
//    cout << "checkpoint" << endl;
    TGraphErrors* g1 = new TGraphErrors(NumHist, x, eNIF, xerr, DeNIF);
    g1->SetNameTitle("eNIF", "Neutron-induced fission detection efficiency");
    SaveToFile("Analysis/Evaluation", g1);
    TGraphErrors* g2 = new TGraphErrors(NumHist, x, y, xerr, yerr);
    g2->SetNameTitle("eSF", "Spontaneaus fission detection efficiency");
    SaveToFile("Analysis/Evaluation", g2);
    TGraphErrors* g3 = new TGraphErrors(NumHist, x, eRel, xerr, DeRel);
    g3->SetNameTitle("eRel", "NIF to SF detection efficiency ratio");
    SaveToFile("Analysis/Evaluation", g3);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
////  Evaluation 4: #FF(Pu), #FF(U), neutron flux ratio, #SF, T2_SF, #U  -->  sigma_nfis        ////
////////////////////////////////////////////////////////////////////////////////////////////////////
void Xsection::Evaluation4()
{

}

void Xsection::SaveToFile(string path, TObject *pObj)
{   //saves a TObject into the selected file fname
    TFile* file = TFile::Open("/home/hoffma93/Go4nfis/offline/results/Evaluation.root", "UPDATE");
    TDirectory *EvalDir;
    TObject *pGraph;
    pGraph = (TObject*) pObj;
    string GraphName = pGraph->GetName();
    //check if folder path already exists, otherwise create it
    if (file->Get(path.c_str())!=0)
        EvalDir = (TDirectory*) file->Get(path.c_str());
    else
    {   file->mkdir(path.c_str(), "Folder containing offline Analysis objects");
        EvalDir = file->GetDirectory(path.c_str());
    }
    EvalDir->cd();
    if (EvalDir->Get(GraphName.c_str())!=0)
        EvalDir->Delete((GraphName+";*").c_str());
    pGraph->Clone()->Write();
    file->Save(); file->Close();
}
