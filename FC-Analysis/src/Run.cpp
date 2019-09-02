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
    CommentFlag = kTRUE;
    DrawSingle = draw;
    DrawMulti = kFALSE;
//    nHist = 0;
    t_live = 0;
    t_real = 0;
//    for (Int_t i = 0; i < NumCh; i++)
//    {
//        nNIF[i] = 0;
//        DnNIF[i] = 0;
//        nSF[i] = 0;
//        DnSF[i] = 0;
//    }
    cout << endl << "Created Run " << Name << endl;
}


Run::~Run()
{

}

Bool_t Run::IsForeground()
{
    return !strcmp(Setup.c_str(), "Open") || !strcmp(Setup.c_str(), "FG");
}


Double_t Run::SetNeutronField(Double_t monitor, Double_t unc_rel, Double_t t_real, Double_t d_start, Double_t t_start/*, Double_t l, Double_t Dl*/)
{
    cout << endl << "Getting neutron monitor data for " << Name << endl;
    Monitor = monitor;
    DMonitor = sqrt(monitor);//unc_rel * monitor;
    t_mon = TMath::Abs(t_real);
    TTimeStamp* ts = new TTimeStamp((UInt_t)d_start, (UInt_t)t_start, 0u, kTRUE, -7200);
    Int_t t = ts->GetSec();
    tStart = t_real < 0 ? t + t_real : t;
    tStop = t_real < 0 ? t : t + t_real;

    Double_t Yield = 2.158821152E4;
    Double_t DYield = 543.3;
    Double_t l = 1500;
    Double_t Dl = 1;
    Double_t r;
    if (CommentFlag)
         cout << " t_real  = " << t_mon << endl
              << " Ch\tn-Fluence[mm^-2]\t n-Flux[mm^-2 s^-1]" << endl;
    for (int i = 0; i < NumCh; i++)
    {
        if (UFC)
            r = l + 134.2 - 10.5 * i; // UFC: 5mm distance, plates 0.4mm each
        else
            r = l + 209.2 - 20.5 * i; // PuFC: 10mm distance, plates 0.4mm each
        NeutronFlux[i] = Yield * monitor / (r*r * t_mon);
        DNeutronFlux[i] = sqrt( pow(DYield * monitor / (r*r * t_mon), 2) +
                                pow(Yield * unc_rel * monitor / (r*r * t_mon), 2) +
                                pow(Yield * monitor * 2*Dl / (r*r*r * t_mon), 2) );
        DstatNeutronFlux[i] = Yield * sqrt(monitor) / (r*r * t_mon); // exact stat. uncertainty would be sqrt(uncorrected monitor)
        if(CommentFlag)
            cout << " " << i+1 << "\t" << NeutronFlux[i] * t_mon << " +- " << DNeutronFlux[i] * t_mon
                 << "\t " << NeutronFlux[i] << " +- " << DNeutronFlux[i] << " +- " << DstatNeutronFlux[i] << endl;
    }
    cout << "Done: Get neutron monitor data" << endl;
    return monitor;
}


void Run::SetToF(string name)
{
    std::stringstream s;
    s << "/home/hoffma93/Programme/Go4nfis/offline/results/" << name << ".root";
    cout << endl << Name << " going to create ToF instance " << s.str() << endl;
    pToF = new ToF(s.str(), FC, Setup, name);
}


void Run::SetLimits(Int_t *pl0, Int_t *pl1, Int_t *pl2, Int_t *pl3)
{
    pToF->SetLimits(pl0, pl1, pl2, pl3);
}


Double_t Run::GetnfoverPhi(Int_t i)
{
    return pToF->GetnfRate(i) / NeutronFlux[i];
}


Double_t Run::GetDnfoverPhi(Int_t i)
{
    return sqrt( pow(pToF->GetDnfRate(i) / NeutronFlux[i], 2) +
                 pow(pToF->GetnfRate(i) * DstatNeutronFlux[i], 2) / pow(NeutronFlux[i], 4) );
}


//void Run::GetHist(Hist *hist)
//{
////    if (CommentFlag)
//        cout << endl << Name << " getting " << hist->Name << endl << " Ch   (n,f)   SF   t_live" << endl;
//    for (Int_t i = 0; i < NumCh; i++)
//    {
//        nNIF[i] += hist->nNIF[i];
//        DnNIF[i] += hist->DnNIF[i];
//        nSF[i] += hist->nSF[i];
//        DnSF[i] += hist->DnSF[i];
//        t_live += hist->t_live;
////        if (CommentFlag)
//            cout << " " << i+1 << "   " << nNIF[i] << "+-" << DnNIF[i] << "   " << nSF[i] << "+-" << DnSF[i] << "   " << t_live << endl;
//    }
//}


//void Run::SetNatoms(Double_t *nAt, Double_t *DnAt)
//{
//    if (CommentFlag)
//        cout << "Run " << Name << " setting atom numbers" << endl << " Ch   atoms" << endl;
//    for (Int_t i = 0; i < NumHist; i++)
//    {
//        nAtoms[i] = nAt[i];
//        DnAtoms[i] = DnAt[i];
//        if (CommentFlag)
//            cout << " " << i+1 << "   " << nAtoms[i] << "+-" << DnAtoms[i] << endl;
//    }
//}


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

//void Run::CrossSection(Int_t i)
//{
////    cout << nNIF[i] << "+-" << DnNIF[i] << endl;
//    uncCS[i] = 1.E22 * nNIF[i] / (t_live * nAtoms[i] * NeutronFlux[i]);
//    DuncCS[i] = uncCS[i] * sqrt( pow(DnNIF[i] / nNIF[i], 2) +
//                                       pow(DnAtoms[i] / nAtoms[i], 2) +
//                                       pow(DstatNeutronFlux[i] / NeutronFlux[i], 2) );
//    if (CommentFlag)
//        cout << "  Cross section. Channel " << i+1 << ", run " << Name << ", nNIF " << nNIF[i] <<
//                ", t_live " << t_live << ", atoms " << nAtoms[i] << ", Flux " << NeutronFlux[i] <<
//                ", sigma " << uncCS[i] << "+-" << DuncCS[i] << endl;
//}
