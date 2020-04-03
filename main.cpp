#include<iostream>
#include "Basics.h"
//#include "ExtendedCk.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TFile.h"
#include "TH1F.h"
//#include "Worksheet.h"

void ComparePionPion(){
    TGraph* grCATS = Basics_PiPiCATS(1,1);
    TGraph* grTheory = Basics_PiPiTheory(1,1);
    TFile* OutputFile = new TFile("ComparePionPion.root","recreate");
    grCATS->Write();
    grTheory->Write();

    delete grCATS;
    delete grTheory;
    delete OutputFile;
}
void DoSomeActuallPiPiFitting(char *argv[]){
    unsigned* settings = new unsigned[3];
    //settings[0] = atoi(argv[1]);
    //settings[1] = atoi(argv[2]);
    //settings[2] = atoi(argv[3]);
    settings[0] = 0;
    settings[1] = 1;

    settings[2] = 0;
    Basics_PiPiCATS2(1,1,settings);

    settings[2] = 1;
    Basics_PiPiCATS2(1,1,settings);

    //settings[2] = 2;
    //Basics_PiPiCATS2(1,1,settings);

    delete [] settings;
}
void FitpLCF(){
     //TH1F* pLCF = Worksheet_ProtonLambda();
     //TFile* OutputFile2 = new TFile("ProtonLambdaUsmaniCF.root","recreate");
     //pLCF->Write();

     //delete pLCF;
     //delete OutputFile2;
}

int main(int argc, char *argv[]){

    printf("\nHello to this lovely CATS tutorial!\n");
    printf(" To find a bug: continue with this tutorial\n");
    printf("-------------------------------------------\n");

    //ComparePionPion();

DoSomeActuallPiPiFitting(argv);
//Compare_AOD_nano();

    //FitpLCF();
    //Basics_ProtonLambda();

    //Ck_pL_Ledni_Usmani();
    //Worksheet_SetUp_pp();
    //Worksheet_SetUp_pXi();
    //Worksheet_SetUp_pSigma();

    return 0;
}
