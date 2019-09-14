#ifndef FC_H
#define FC_H

Double_t PeakCenter(Int_t ch, string FC = "PuFC")
{
    if (!strcmp(FC.c_str(), "PuFC"))
    { // PuFC
        Double_t m[] = {139.7705, 128.2225, 125.879, 129.663, 129.5655, 128.589, 127.466, 127.49};
        return m[ch];
    }
    else
    {
        Double_t m[] = {273.4865, 268.945, 260.5955, 265.2345, 265.21, 264.3555, 263.501, 263.086};
        return m[ch];
    }
}

Double_t Distance(Int_t ch, string FC = "PuFC")
{
    if (!strcmp(FC.c_str(), "PuFC"))
    { // PuFC
        return 1709.2 - ch * 20.5;
    } else { // UFC
        return 1634.2 - ch * 10.5;
    }
}

Double_t SolidAngle(Double_t Angle)
{
    return 2*TMath::Pi() * (1.0 - cos(Angle));
}

Double_t SolidAngle(Double_t D, Double_t R)
{
    return 2*TMath::Pi() * (1.0 - cos(atan(R / D)));
}

Double_t SolidAngle(Int_t ch, string FC = "PuFC")
{
    Double_t DepositRadius = 37.0;
    return SolidAngle(Distance(ch, FC), DepositRadius);
}

void FC()
{
    string FC = "PuFC";
    for (Int_t i = 0; i < 8; i++)
        cout << " " << i+1 << "   " << Distance(i, FC) <<  "   " << SolidAngle(i, FC) << endl;
}

#endif
