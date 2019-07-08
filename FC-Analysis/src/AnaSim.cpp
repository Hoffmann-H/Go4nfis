#include "AnaSim.h"
using namespace std;

AnaSim::AnaSim(string fc, Bool_t use_track_id, Plot *p)
{
    cout << endl << "Creating simulation analysis instance for " << fc << endl;

    FC = fc;
    PuFC = strcmp(fc.c_str(), "UFC");
    if (PuFC)
    {
        FgPath = "/home/hoffma93/Programme/Geant4-Work/builds/G4PuFCvsH19/results/4_ene/PuFC_Open_5E7.root";
        BgPath = "/home/hoffma93/Programme/Geant4-Work/builds/G4PuFCvsH19/results/4_ene/PuFC_SB_5E7.root";
    } else {
        FgPath = "/home/hoffma93/Programme/Geant4-Work/builds/G4PuFCvsH19/results/4_ene/UFC_Open_5E7.root";
        BgPath = "/home/hoffma93/Programme/Geant4-Work/builds/G4PuFCvsH19/results/4_ene/UFC_SB_5E7.root";
    }

    tID = use_track_id;
    DoneCorrections = kFALSE;

    DrawSingle = kTRUE;
    DrawMulti = kTRUE;
    if (p != 0)
        plot = p;
//    else {
        DrawSingle = kFALSE;
        DrawMulti = kFALSE;
//    }
    CommentFlag = kFALSE;

    cout << "Done: Create simulation analysis instance for " << FC << ", " << Setup << endl;
}


AnaSim::~AnaSim()
{

}


void AnaSim::Corrections()
{// Calculate correction factors for data analysis from simulation.
    // T := Direct incident n / emitted neutrons towards deposit
    // S := fissions induced by direct / induced fissions (both in time gate)
    //    = direct fissions / ( direct fissions + scattered fissions )

    Fg = new Sim(FgPath, FC, "Open", tID, DrawSingle);
    Fg->Calculate();
    cout << endl << "Correction factors..." << endl;
//    char name[64] = "";

//    cout << "Input for correction factors" << endl;
    for (int i = 0; i < NumCh; i++)
    {
//        cout << " " << Fg->nDirect[i] << " " << Fg->DnDirect[i] << " " << Fg->nProj[i] << " " << Fg->effDirect[i] << " " << Fg->effScat[i] << endl;
        //// Calculate Transmission and Scattering correction factor
        T[i] = Fg->nDirect[i] / Fg->nProj[i];
        DT[i] = Fg->DnDirect[i] / Fg->nProj[i];
        S[i] = Fg->effDirect[i] / (Fg->effScat[i] + Fg->effDirect[i]);

        //// Calculate T&S's correction factor
        F[i] = Fg->nProj[i] / (Fg->effScat[i] + Fg->effDirect[i]) * Fg->effDirect[i] / Fg->nDirect[i];   // == S[i] / T[i]
    }
        //// Calculate correlated uncertainties
        Uncertainties();

    for (int i = 0; i < NumCh; i++)
    {

        cout << " " << i+1 << "  " << T[i] << "+-" << DT[i] << //endl;
                              "  " << S[i] << "+-" << DS[i] << //"+-" << sqrt(D2S_sys) <<
                              "  " << F[i] << "+-" << DF[i] << endl;//"+-" << sqrt(D2S_sys) / T[i] << endl;
    }

    DoneCorrections = kTRUE;
    cout << "Done: Corrections" << endl;
    if (!DrawSingle)
        return;
    plot->SimF(T, DT, S, DS, F, DF);
}


void AnaSim::Uncertainties()
{
    Double_t D2F[NumCh]; Double_t D2S[NumCh];
    for (Int_t i = 0; i < NumCh; i++)
    { D2F[i] = 0; D2S[i] = 0; }

    Double_t E, w, Dw, n;
    for (Int_t binE = Fg->binEmin; binE <= Fg->binEmax; binE++) // Iterate over all energy bins containing direct neutrons
    {
        E = Fg->pH2TvsE[0][0]->GetYaxis()->GetBinCenter(binE);
//        cout << E << " ";
        Fg->relSigma(E, w, Dw);
        for (Int_t i = 0; i < NumCh; i++) // Iterate over channels
        {
            n = tID ? Fg->pH1Eproj[i][0]->GetBinContent(binE) - Fg->pH1Eproj[i][1]->GetBinContent(binE)
                    : Fg->pH1Eproj[i][0]->GetBinContent(binE);
            D2F[i] += n * pow(dF_dnE(w, i), 2); // (sqrt(n_E) * dF/dn_E)^2
            D2S[i] += n * pow(dS_dnE(w, i), 2);

//            cout << n << " ";
        }
//        cout << w << endl;
    }
    for (Int_t i = 0; i < NumCh; i++)
    {
        Double_t dF_dSc = - Fg->effDirect[i] * Fg->nProj[i]
                          / pow(Fg->effScat[i] + Fg->effDirect[i], 2)
                          / Fg->nDirect[i]; // Derivation of correction factor F wrt number of scattered fissions Fg->effScat[i].
        DF[i] = sqrt(D2F[i] + pow(Fg->DeffScat[i] * dF_dSc, 2));
        Double_t dS_dSc = - Fg->effDirect[i]
                          / pow(Fg->effScat[i] + Fg->effDirect[i], 2);
        DS[i] = sqrt(D2S[i] + pow(Fg->DeffScat[i] * dS_dSc, 2));

//        cout << Fg->DeffScat[i] << "   " << sqrt(D2F[i]) << endl;
    }
    return;
}


Double_t AnaSim::dF_dnE(Double_t w, Int_t ch)
{ // Returns the derivation of correction factor F wrt n_E
    return Fg->nProj[ch] * ( (w * Fg->nDirect[ch] - Fg->effDirect[ch]) * Fg->effScat[ch] - pow(Fg->effDirect[ch], 2) )
                         / ( pow((Fg->effScat[ch] + Fg->effDirect[ch]) * Fg->nDirect[ch], 2) );
}


Double_t AnaSim::dF_dwE(Double_t n, Int_t ch)
{ // Returns the deviation of correction factor F wrt w_E
    return 0;
}


Double_t AnaSim::dS_dnE(Double_t w, Int_t ch)
{ // Returns the derivation of scattering factor S wrt n_E
    return w * Fg->effScat[ch] / pow(Fg->effDirect[ch] + Fg->effScat[ch], 2);
}


Double_t AnaSim::dS_dwE(Double_t n, Int_t ch)
{ // Returns the deviation of scattering factor S wrt w_E
    return 0;
}


void AnaSim::ShadowCone()
{
    Bg = new Sim(BgPath, FC, "SB", tID, DrawSingle);
    Bg->Calculate();
    cout << endl << "Analyzing shadow cone simulation..." << endl;
    cout << " Ch   S(shadow cone)" << endl;
    for (Int_t i = 0; i < NumCh; i++)
    {
        SC[i] = 1.0 - (Bg->effDirect[i] + Bg->effScat[i]) / (Fg->effDirect[i] + Fg->effScat[i]) * Fg->nProj[i] / Bg->nProj[i];
        DSC[i] = (1.0 - SC[i]) * sqrt( (pow(Bg->DeffDirect[i], 2) + pow(Bg->DeffScat[i], 2)) / pow(Bg->effDirect[i] + Bg->effScat[i], 2) +
                               (pow(Fg->DeffDirect[i], 2) + pow(Fg->DeffScat[i], 2)) / pow(Fg->effDirect[i] + Fg->effScat[i], 4) );
        Double_t DSC2 = sqrt( pow(Bg->DeffDirect[i] / (Fg->effDirect[i] + Fg->effScat[i]) * Fg->nProj[i] / Bg->nProj[i], 2) +
                              pow(Bg->DeffScat[i] / (Fg->effDirect[i] + Fg->effScat[i]) * Fg->nProj[i] / Bg->nProj[i], 2) +
                              pow(Fg->DeffDirect[i] * (Bg->effDirect[i] + Bg->effScat[i]) / pow(Fg->effDirect[i] + Fg->effScat[i], 2) * Fg->nProj[i] / Bg->nProj[i], 2) +
                              pow(Fg->DeffScat[i] * (Bg->effDirect[i] + Bg->effScat[i]) / pow(Fg->effDirect[i] + Fg->effScat[i], 2) * Fg->nProj[i] / Bg->nProj[i], 2) );
        cout << " " << i+1 << "   " << SC[i] << "+-" << DSC[i] << "(" << DSC2 << ")" << endl;
    }
    cout << "Done: shadow cone simulation" << endl;
}


//void AnaSim::SaveToFile(string RootPath, TObject *pObj)
//{   //saves a TObject into member TFile *f
//    f->ReOpen("UPDATE");
//    TDirectory *EvalDir;
//    TObject *pGraph;
//    pGraph = (TObject*) pObj;
//    string GraphName = pGraph->GetName();
//    //check if folder path already exists, otherwise create it
//    if (f->Get(RootPath.c_str())!=0)
//        EvalDir = (TDirectory*) f->Get(RootPath.c_str());
//    else
//    {   cout << " Creating root directory " << RootPath << endl;
//        f->mkdir(RootPath.c_str(), "Folder");
//        EvalDir = f->GetDirectory(RootPath.c_str());
//    }
//    EvalDir->cd();
//    if (EvalDir->Get(GraphName.c_str())!=0)
//        EvalDir->Delete((GraphName+";*").c_str());
//    pGraph->Clone()->Write();
//    f->Save(); //file->Close();
//    cout << " Saved " << GraphName << endl;
//}


TH1D* AnaSim::CopyRange(TH1D* pH, char* name, Double_t x0, Double_t x1, Double_t yoffset)
{
    Double_t offset = pH->GetBinLowEdge(0);
    Double_t ChPerBin = pH->GetBinWidth(0);
    int b0 = (x0 - offset) / ChPerBin;
    int b1 = (x1 - offset) / ChPerBin + 1;
    TH1D* pH2 = new TH1D(name, name, b1 - b0, pH->GetBinLowEdge(b0), pH->GetBinLowEdge(b1));
    for (int bin = 0; bin < pH2->GetNbinsX(); bin++)
        pH2->SetBinContent(bin + 1, pH->GetBinContent(bin + b0) + yoffset);
    return pH2;
}

//*/
