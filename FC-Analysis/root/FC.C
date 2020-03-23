#ifndef FC_H
#define FC_H
#include <fstream>
/// hard coded ToF gates
//#define GATE 10
#define LEFT 15
#define RIGHT 40
#define TAIL 40
using namespace std;

Int_t Left(string FC)
{
    if (!strcmp(FC.c_str(), "PuFC"))
        return 15;//7;//
    else // UFC
        return 15;
}

Int_t Right(string FC)
{
    if (!strcmp(FC.c_str(), "PuFC"))
        return 35;//7;//
    else // UFC
        return 40;
}

Int_t Tail(string FC)
{
    if (!strcmp(FC.c_str(), "PuFC"))
        return 40;
    else // UFC
        return 50;
}

Double_t PeakCenter(Int_t ch, string FC = "PuFC")
{
    if (!strcmp(FC.c_str(), "PuFC"))
    { // PuFC
        Double_t m[] = {139.7705, 128.2225, 125.879, 129.663, 129.5655, 128.589, 127.466, 127.49};
//        Double_t m[] = {32.9355, 32.5355, 32.1294, 31.7217, 31.3162, 30.9138, 30.5155, 30.1205};
        return m[ch];
    }
    else
    {
        Double_t m[] = {273.4865, 268.945, 260.5955, 265.2345, 265.21, 264.3555, 263.501, 263.086};
//        Double_t m[] = {31.7416, 31.5622, 31.286, 31.2567, 30.8611, 30.6516, 30.4355, 30.2321};
        return m[ch];
    }
}

Int_t Gate_0(Int_t ch, string FC = "PuFC")
{
    return 62;
}

Int_t Gate_a(Int_t ch, string FC = "PuFC", Int_t l = 15)
{
    return (Int_t)(PeakCenter(ch, FC)+0.5) - (l ? l : Left(FC));
}

Int_t Gate_1(Int_t ch, string FC = "PuFC", Int_t l = 0)
{
    return (Int_t)(PeakCenter(ch, FC)+0.5) - (l ? l : Left(FC));
}

Int_t Gate_2(Int_t ch, string FC = "PuFC", Int_t r = 0)
{
    return (Int_t)(PeakCenter(ch, FC)+0.5) + (r ? r : Right(FC));
}

Int_t Gate_b(Int_t ch, string FC = "PuFC", Int_t t = 0)
{
    return (Int_t)(PeakCenter(ch, FC)+0.5) + (t ? t : Tail(FC));
}

Int_t Gate_3(Int_t ch, string FC = "PuFC")
{
    return 402;
}

Double_t DepositRadius(Int_t ch = 0, string FC = "PuFC")
{
    return 37.0; // mm
}

Double_t DepositArea(Int_t ch = 0, string FC = "PuFC")
{
    return TMath::Pi() * pow(DepositRadius(ch, FC), 2); // mm^2
}

Double_t Distance(Int_t ch, string FC = "PuFC")
{ // in mm
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

Double_t SolidAngle(Int_t ch, string FC)
{
    return SolidAngle(Distance(ch, FC), DepositRadius(ch, FC));
}

Double_t QDCcut(Int_t ch, string FC)
{
    if (strcmp(FC.c_str(), "PuFC"))
    { // UFC
        Double_t qdc_min[]   = {220.396, 500.964, 219.016, 214.101, 199.19,  171.929, 169.046, 121.56}; // UFC
        return qdc_min[ch];
    } else { // PuFC
        Double_t qdc_min[]   = {899.24,  853.668, 895.652, 849.393, 1046.41, 891.396, 906.123, 837.486}; // PuFC
        return qdc_min[ch];
    }
}

Int_t digit(Double_t err)
{ // return Nachkommastellen to print uncertainty err
    Int_t Digit = (Int_t)-TMath::Floor(TMath::Log10(err));
    if (TMath::Floor(err * pow(10, Digit)) < 3)
        Digit++;
    return Digit;
}

Double_t round(Double_t val, Int_t digit)
{
    return  TMath::Floor(val * pow(10, digit) + 0.5) / pow(10, digit);
}

string pm(Double_t val, Double_t err)
{
    Int_t Digit = digit(err);
    char name[64] = "";
    char format[64] = "";
    sprintf(format, "%s%i%s%i%s", "%.", Digit > 0 ? Digit : 0, "f #pm %.", Digit > 0 ? Digit : 0, "f");
    sprintf(name, format, round(val, Digit), round(err, Digit));
    string s(name);
    return s;
}

string br(Double_t val, Double_t err)
{ // brackets notation. Does not make sense for large numbers: 11590000(28)
    Int_t Digit = digit(err);
    char name[64] = "";
    char format[64] = "";
    sprintf(format, "%s%i%s", "%.", Digit > 0 ? Digit : 0, "f(%i)");
    sprintf(name, format, round(val, Digit), (Int_t)(err * pow(10, Digit) + 0.5));
    string s(name);
    return s;
}

void FC()
{
    string FC = "PuFC";
    for (Int_t i = 0; i < 8; i++)
        cout << " " << i+1 << "   " << Distance(i, FC) <<  "   " << SolidAngle(i, FC) <<  "   " << Gate_1(i, FC) <<  "   " << Gate_2(i, FC) << endl;
//    cout << br(1.159E+19, 2.7813E+15) << endl;
    cout << Gate_b(0) << endl;
}

#endif
