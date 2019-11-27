#ifndef DEPOSIT_H
#define DEPOSIT_H
#include "fstream"
#define GridX 25
#define GridY 25
#define ScaleX 4
#define ScaleY 2
#define SizeX 200
#define SizeY 100
#define MaxZ 25

TH2D* ReadTxtMatrix(string file)
{
    std::ifstream ifs (file.c_str(), std::ifstream::in);
    if (!ifs.is_open())
        return 0;
    TH2D *pH2 = new TH2D("name", "title", GridX, 0, 80, GridY, 0, 80);
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
             pH2->SetBinContent(x+1, y+1, z);
             y++; }
         x++;
     }
     ifs.close();
     pH2->Scale(1.0/pH2->GetMaximum());
     return pH2;
}

Double_t Shift(TH2D *srfc, Double_t column, Int_t row)
{
    Double_t x = (column - 0.5 * (SizeX - GridX * ScaleX)) / (Double_t)ScaleX;
    Double_t y = (row - 0.5 * (SizeY - GridY * ScaleY)) / (Double_t)ScaleY;
    Double_t z;
    if (x < 0 || x > GridX || y < 0 || y > GridY)
        z = 0.0;
    else
        z = srfc->GetBinContent((Int_t)x, (Int_t)y);
//    if (z) cout << x << " " << y << " " << z << endl;
    return 20.0 * (100.0 - MaxZ * z) / (200.0 - MaxZ * z);
}

string Line(TH2D *srfc, Int_t row)
{
    std::stringstream s;
    Int_t N[10];
    for (Int_t column = 0; column < 10; column++) {
        N[column] = (Int_t)TMath::Floor(TMath::Abs(1000.0*sin(SizeY*column+row))) % 10;
        s << N[column]; }
    for (Int_t column = 10; column < SizeX; column++)
    {
        Double_t c = (Double_t)column;
        while(c >= 10)
            c -= (Int_t)Shift(srfc, column, row);
        s << N[(Int_t)c];
//        s << (Int_t)Shift(srfc, column, row)-1;
    }
    return s.str();
}

void Deposit(string isotope = "Pu242", Int_t ch = 0)
{
    TH2D *srfc = ReadTxtMatrix("/home/hoffma93/Experiment/Autoradiographien/"+isotope+"/Result_"+to_string(ch+1)+".txt");
    for (Int_t row = 0; row < SizeY; row++)
        cout << Line(srfc, row) << endl;

//    new TCanvas();
//    srfc->Draw("colz");
}
#endif
