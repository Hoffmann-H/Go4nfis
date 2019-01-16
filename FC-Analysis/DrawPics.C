

int TreatFile(string file_name)
{
    TCanvas* pCRawQDC[8];
//    TCanvas* pCAnaQDC[8];
    TFile* f = TFile::Open(file_name.c_str());
    /// Draw QDC fits
    char hname[64] = "";
    for (int i = 0; i < 8; i++)
    {
        /// raw QDC spectrum
        // open histogram
        sprintf(hname, "/Histograms/Raw/QDC/low/H1RawQDCl_%i", i + 1);
        TH1I *pHRawQDC = (TH1I*)f->Get(hname);
        // open fit
        sprintf(hname, "/Histograms/Raw/QDC/low/fit/fPed_%i", i+1);
        TF1* fPed = (TF1*)f->Get(hname);
        sprintf(hname, "/Histograms/Raw/QDC/low/fit/fCut_%i", i+1);
        TF1* fCut = (TF1*)f->Get(hname);
        sprintf(hname, "/Histograms/Raw/QDC/low/fit/fMax_%i", i+1);
        TF1* fMax = (TF1*)f->Get(hname);
        // draw together
        pCRawQDC[i] = new TCanvas();
        pHRawQDC->Draw();
        fPed->Draw("same");
        fCut->Draw("same");
        fMax->Draw("same");

        /// self-triggered ToF-gated QDC
        // ...
//        TH1I *pHAnaQDC;



    }
    return 1;
}

int DrawPics()
{
    TreatFile("/home/hoffma93/Go4nfis/offline/results/NIF.root");
//    TFile* f = TFile::Open("/home/hoffma93/Go4nfis/offline/test/NIF.root");

    return 1;
}
