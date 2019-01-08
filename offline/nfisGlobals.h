#ifndef NFISGLOBALS_H
#define NFISGLOBALS_H

#define NumQDC 8                // Number of QDC channel
#define NumQDCReg 4096          // Number of QDC register
#define NumTDC 32               // Number of TDC channel
#define NumVeto 5               // Number of Veto channel
#define NumScaler 64            // Number of Scaler channel
#define NumHZDRFC 8             // Number of HZDR fission chamber channels

#define PI 3.141592654

#define MAXEVENTS 33            //Maximum number of trigger-evnts / readout-evnt
#define MAXHITS 32              //Maximum number of TDC hits per event & channel
#define MAXSTRINGLENGTH 256     //Maximum length of strings
#define MAXLINELENGTH 1024
#define MAXSCALERCHANNEL 32500  // Maximum channel in scaler histogram

#define TDCCUTOFF -600000       // define cut-off value for TDC histograms


#endif // NFISGLOBALS_H
