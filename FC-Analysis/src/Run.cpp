#include "Run.h"

using namespace std;

Run::Run(string FC, string setup, Int_t nr, Bool_t draw)
{
    this->FC = FC;
    UFC = strcmp(FC.c_str(), "PuFC");
    Setup = setup;
    std::stringstream s;
    s << FC << "_" << setup << "_" << nr+1;
    Name = s.str();
    CommentFlag = kFALSE;
    DrawSingle = draw;
    DrawMulti = kFALSE;
    nHist = 0;
    t_live = 0;
    t_real = 0;
    for (Int_t i = 0; i < NumCh; i++)
    {
        nNIF[i] = 0;
        DnNIF[i] = 0;
        nSF[i] = 0;
        DnSF[i] = 0;
    }
    cout << endl << "Created Run " << Name << endl;
}


Run::~Run()
{

}

Bool_t Run::IsForeground()
{
    return !strcmp(Setup.c_str(), "Open") || !strcmp(Setup.c_str(), "FG");
}


Double_t Run::SetNeutronField(Double_t monitor, Double_t unc_rel, Double_t t_real, Double_t l, Double_t Dl)
{
    cout << endl << "Getting neutron monitor data for " << Name << endl;
    Monitor = monitor;
    t_mon = t_real;
    Double_t Yield = 2.15882E4;
    Double_t DYield = 543.3;
    Double_t r;
    if (CommentFlag)
         cout << " t_real  = " << t_mon << endl
              << " Ch   n-Fluence[mm^-2]   n-Flux[mm^-2 s^-1]" << endl;
    for (int i = 0; i < NumCh; i++)
    {
        if (UFC)
            r = l + 119 - 10.8 * i; // UFC: 5mm distance, plates 0.4mm each
        else
            r = l + 191 - 20.8 * i; // PuFC: 10mm distance, plates 0.4mm each
        NeutronFlux[i] = Yield * monitor / (r*r * t_mon);
        DNeutronFlux[i] = sqrt( pow(DYield * monitor / (r*r * t_mon), 2) +
                                pow(Yield * unc_rel * monitor / (r*r * t_mon), 2) +
                                pow(Yield * monitor * 2*Dl / (r*r*r * t_mon), 2) );
        DstatNeutronFlux[i] = Yield * sqrt(monitor) / (r*r * t_mon); // exact stat. uncertainty would be sqrt(uncorrected monitor)
        if(CommentFlag)
            cout << " " << i+1 << "   " << NeutronFlux[i] * t_mon << " +- " << DNeutronFlux[i] * t_mon
                 << "   " << NeutronFlux[i] << " +- " << DNeutronFlux[i] << " +- " << DstatNeutronFlux[i] << endl;
    }
    cout << "Done: Get neutron monitor data" << endl;
    return monitor;
}


void Run::GetHist(Hist *hist)
{
//    if (CommentFlag)
        cout << endl << Name << " getting " << hist->Name << endl << " Ch   (n,f)   SF   t_live" << endl;
    for (Int_t i = 0; i < NumCh; i++)
    {
        nNIF[i] += hist->nNIF[i];
        DnNIF[i] += hist->DnNIF[i];
        nSF[i] += hist->nSF[i];
        DnSF[i] += hist->DnSF[i];
        t_live = hist->t_live;
//        if (CommentFlag)
            cout << " " << i+1 << "   " << nNIF[i] << "+-" << DnNIF[i] << "   " << nSF[i] << "+-" << DnSF[i] << "   " << t_live << endl;
    }
}


void Run::SetNatoms(Double_t *nAt, Double_t *DnAt)
{
    if (CommentFlag)
        cout << "Run " << Name << " setting atom numbers" << endl << " Ch   atoms" << endl;
    for (Int_t i = 0; i < NumHist; i++)
    {
        nAtoms[i] = nAt[i];
        DnAtoms[i] = DnAt[i];
        if (CommentFlag)
            cout << " " << i+1 << "   " << nAtoms[i] << "+-" << DnAtoms[i] << endl;
    }
}


//void Run::AnalyzeDt(Int_t i, Int_t lim0, Int_t lim1, Int_t lim2, Int_t lim3)
//{
//    nNIF[i] = 0;
//    Double_t D2nNIF = 0;
//    nSF[i] = 0;
//    Double_t D2nSF = 0;
//    for (Int_t j = 0; j < nHist; j++)
//    {
//        pH[j]->AnalyzeDtPeak(i, lim0, lim1, lim2, lim3);
//        nNIF[i] += pH[j]->nNIF[i];
//        D2nNIF += pow(DnNIF[i], 2);
//        nSF[i] += pH[j]->nSF[i];
//        D2nSF += pow(DnSF[i], 2);
//    }
//}

void Run::CrossSection(Int_t i)
{
//    cout << nNIF[i] << "+-" << DnNIF[i] << endl;
    uncCS[i] = 1.E22 * nNIF[i] / (t_live * nAtoms[i] * NeutronFlux[i]);
    DuncCS[i] = uncCS[i] * sqrt( pow(DnNIF[i] / nNIF[i], 2) +
                                       pow(DnAtoms[i] / nAtoms[i], 2) +
                                       pow(DstatNeutronFlux[i] / NeutronFlux[i], 2) );
    if (CommentFlag)
        cout << "  Cross section. Channel " << i+1 << ", run " << Name << ", nNIF " << nNIF[i] <<
                ", t_live " << t_live << ", atoms " << nAtoms[i] << ", Flux " << NeutronFlux[i] <<
                ", sigma " << uncCS[i] << "+-" << DuncCS[i] << endl;
}
