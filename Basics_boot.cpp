
#include "Basics.h"
#include "CATS.h"
#include "CATStools.h"
#include "CATSconstants.h"
#include "TGraph.h"
#include "TString.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include <vector>
#include "DLM_CkDecomposition.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "TDirectory.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TList.h"
#include "TStyle.h"
#include "TH1.h"
#include <iostream>
#include <cmath>
#include "DLM_Random.h"
#include "TNtuple.h"
#include "TFrame.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TMath.h"


CATS* KITTY_CATS_FIT_PL;

double Basics_Potential_Usmani(double* pars){
    double r = pars[0];
    //pars[1] is the momentum
    //Values for the potential
    double vbar = 6.2;//6.2
    double vsigma = 0.25;//0.25

    double x=r*0.7;
    double vc = pars[3]/(1+exp((r-pars[4])/pars[5]));
    double tpi = (1.0+3.0/x+3.0/(x*x)) * (exp(-x)/x) * pow(1.-exp(-2.*r*r),2.);

    double v = 0.;

    double& Spin = pars[2];
    if (Spin == 0) v = vc - (vbar + 0.75*vsigma)*tpi*tpi;//Usmani singlet
    else if (Spin == 1)  v = vc - (vbar - 0.25*vsigma)*tpi*tpi;//Usmani triplet
    else printf ("wrong polarization\n");

    return v;
}
double Basics_PionCoulombCorrection(const bool& Identical){

}

//computes the expectation based on quantum statistics only
//the formula used is from chapter 4.1 in Phys. Rev. C 96 (2017) 064908 (ATLAS paper on pipi correlations in p-Pb)
TGraph* Basics_PiPiTheory(const bool& Identical, const bool& WithCoulomb){
    //load from a Mathematica output file
    FILE *InFileCBE;
    const TString CBEname = "../Mathematica/tab_txCBE.dat";
    InFileCBE = fopen(CBEname.Data(), "r");
    if(!InFileCBE){
        printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", CBEname.Data());
        return NULL;
    }
    char*  cdummy = new char [512];
    float kstar;
    float ckval;
    TGraph gCBE;
    gCBE.SetName("gCBE");
    unsigned NumPointsCBE=0;
    while(!feof(InFileCBE)){
        if(!fgets(cdummy, 511, InFileCBE)) continue;
        sscanf(cdummy, "%f %f",&kstar, &ckval);
        gCBE.SetPoint(NumPointsCBE,kstar,ckval);
        NumPointsCBE++;
    }
    fclose(InFileCBE);

    FILE *InFileK;
    const TString Kname = "../Mathematica/tab_txCoulombSame.dat";
    InFileK = fopen(Kname.Data(), "r");
    TGraph gK;
    gK.SetName("gK");
    unsigned NumPointsK=0;
    if(!InFileK){
        printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", Kname.Data());
        return NULL;
    }
    while(!feof(InFileK)){
        if(!fgets(cdummy, 511, InFileK)) continue;
        sscanf(cdummy, "%f %f",&kstar, &ckval);
        gK.SetPoint(NumPointsK,kstar,ckval);
        NumPointsK++;
    }
    fclose(InFileK);

    TGraph* gCk = new TGraph();
    gCk->SetName("gPiPiTheory");
    double dkstar,dckval;
    for(unsigned uBin=0; uBin<NumPointsCBE; uBin++){
        gCBE.GetPoint(uBin,dkstar,dckval);
        gCk->SetPoint(uBin,dkstar,dckval*gK.Eval(dkstar));
    }

    delete [] cdummy;
    return gCk;
}

TGraph* Basics_PiPiCATS(const bool& Identical, const bool& WithCoulomb){

    const unsigned NumMomBins = 200;
    const double kMin = 0;
    const double kMax = 100;

    CATS PionKitty;
    //(#bins,min,max) or (#bins,BinRangeArray[] as in ROOT)
    PionKitty.SetMomBins(NumMomBins,kMin,kMax);

    //object for the parameters to be used by the source function
    CATSparameters SOURCE_PARS(CATSparameters::tSource,1,true);
    //set the first and only par (source size)
    SOURCE_PARS.SetParameter(0,2.2);
    //say to CATS which Source function to use, and with which parameter set
    PionKitty.SetAnaSource(GaussSource,SOURCE_PARS);
    PionKitty.SetUseAnalyticSource(true);

    //standard settings for a CATS object which has no strong interaction potential included
    PionKitty.SetNumChannels(1);

    //which channel, how many PWs
    PionKitty.SetNumPW(0,0);
    //which channel, spin value
    PionKitty.SetSpin(0,0);
    //which channel, weight
    PionKitty.SetChannelWeight(0, 1);

    //include the coulomb interaction. Q1Q2 is the multiplied charge numbers of the two particles
    PionKitty.SetQ1Q2(1);
    //the PdgId is needed when using a source from a transport model
    //actually the only important thing here is if the particles are identical
    if(Identical) PionKitty.SetPdgId(211, 211);
    else PionKitty.SetPdgId(-211, 211);

    if(WithCoulomb) PionKitty.SetRedMass( 0.5*Mass_Pic );
    else if(Identical) PionKitty.SetRedMass( 0.5*Mass_Pi0 );
    else PionKitty.SetRedMass( (Mass_Pi0*Mass_Pic)/(Mass_Pi0+Mass_Pic) );

    PionKitty.KillTheCat();

    TGraph* grCk = new TGraph();
    grCk->SetName("grCk");
    grCk->Set(NumMomBins);
    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
        grCk->SetPoint(uMom,PionKitty.GetMomentum(uMom),PionKitty.GetCorrFun(uMom));
    }
    //Test
    //std::vector<double> LambdaHousehold = CalculateLambdaParam();

    return grCk;
}

double CATS_FIT_PL(double* x, double* pars){
    double& Norm = pars[0];
    double& LambdaPar = pars[1];
    double& SourceSize = pars[2];
    //set the radius to the fit value
    //the last parameter says its a small change of the radius, which does not require a new computing grid (saves time)
    //however, the last step requires making sure a good initial value of the radius
    KITTY_CATS_FIT_PL->SetAnaSource(0,SourceSize,true);
    //useful tip: this makes CATS to shut up and not flood your screen. Only errors will be displayed. Use nSilent to suppress even those
    KITTY_CATS_FIT_PL->SetNotifications(CATS::nError);
    KITTY_CATS_FIT_PL->KillTheCat();
    return Norm*(LambdaPar*KITTY_CATS_FIT_PL->EvalCorrFun(*x)+1.-LambdaPar);
}

void Basics_ProtonLambda(){

    const unsigned NumMomBins = 60;
    const double kMin = 0;
    const double kMax = 240;

    //object for the parameters to be used by the source function
    CATSparameters SOURCE_PARS(CATSparameters::tSource,1,true);
    SOURCE_PARS.SetParameter(0,1.3);

    CATS Kitty_pL;
    Kitty_pL.SetMomBins(NumMomBins,kMin,kMax);
    Kitty_pL.SetAnaSource(GaussSource,SOURCE_PARS);
    Kitty_pL.SetAnaSource(0,SOURCE_PARS.GetParameter(0));
    Kitty_pL.SetUseAnalyticSource(true);
    Kitty_pL.SetMomentumDependentSource(false);
    Kitty_pL.SetThetaDependentSource(false);
    //should you include in the result any bins, where the Schroedinger solver failed
    Kitty_pL.SetExcludeFailedBins(false);
    Kitty_pL.SetQ1Q2(0);
    Kitty_pL.SetPdgId(2212, 3122);
    Kitty_pL.SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );
    Kitty_pL.SetNumChannels(2);
    Kitty_pL.SetNumPW(0,1);
    Kitty_pL.SetNumPW(1,1);
    Kitty_pL.SetSpin(0,0);
    Kitty_pL.SetSpin(1,1);
    Kitty_pL.SetChannelWeight(0, 1./4.);
    Kitty_pL.SetChannelWeight(1, 3./4.);

    CATSparameters POT_PARS_1S0(CATSparameters::tPotential,4,true);
    POT_PARS_1S0.SetParameter(0,0);
    POT_PARS_1S0.SetParameter(1,2137);
    POT_PARS_1S0.SetParameter(2,0.5);
    POT_PARS_1S0.SetParameter(3,0.2);
    CATSparameters POT_PARS_3S1(CATSparameters::tPotential,4,true);
    POT_PARS_3S1.SetParameter(0,1);
    POT_PARS_3S1.SetParameter(1,2137);
    POT_PARS_3S1.SetParameter(2,0.5);
    POT_PARS_3S1.SetParameter(3,0.2);
    //WhichChannel,WhichPW,PotentialFunction,Parameters
    Kitty_pL.SetShortRangePotential(0,0,Basics_Potential_Usmani,POT_PARS_1S0);
    Kitty_pL.SetShortRangePotential(1,0,Basics_Potential_Usmani,POT_PARS_3S1);
    Kitty_pL.SetMaxNumThreads(4);
    Kitty_pL.KillTheCat();

    KITTY_CATS_FIT_PL = &Kitty_pL;

    TFile* fInput = new TFile("../Files/DummyProtonLambda.root","read");
    TH1F* hInput = (TH1F*)fInput->Get("hDummyProtonLambda");
    TF1* FITTER = new TF1("FITTER",CATS_FIT_PL,kMin,kMax,6);
    FITTER->SetParameter(0,1);FITTER->SetParLimits(0,0.5,2);
    FITTER->SetParameter(1,0.5);FITTER->SetParLimits(1,0,1);
    FITTER->SetParameter(2,SOURCE_PARS.GetParameter(0));FITTER->SetParLimits(2,0.5,3);
    //FITTER->SetParameter(3,POT_PARS_1S0.GetParameter(1));FITTER->SetParLimits(3,0.99*POT_PARS_1S0.GetParameter(1),1.01*POT_PARS_1S0.GetParameter(1));
    //FITTER->SetParameter(4,POT_PARS_1S0.GetParameter(2));FITTER->SetParLimits(4,0.99*POT_PARS_1S0.GetParameter(2),1.01*POT_PARS_1S0.GetParameter(2));
    //FITTER->SetParameter(5,POT_PARS_1S0.GetParameter(3));FITTER->SetParLimits(5,0.99*POT_PARS_1S0.GetParameter(3),1.01*POT_PARS_1S0.GetParameter(3));

    //FITTER->FixParameter(1,0.6);
    FITTER->FixParameter(3,POT_PARS_1S0.GetParameter(1));
    FITTER->FixParameter(4,POT_PARS_1S0.GetParameter(2));
    FITTER->FixParameter(5,POT_PARS_1S0.GetParameter(3));


    hInput->Fit(FITTER,"S, N, R, M");

    delete FITTER;
    delete fInput;
}















//Outside definition in order for the following to work
    DLM_CkDecomposition* CkDe_PiPi;
    CATS* KITTY_CATS_FIT_PiPi;

void SetTheLambdaParamsToDML_CkDecompObj (DLM_CkDecomposition* CkDe_PiPi, bool verbose=true) {
    //Generate all the Lambdaparams
    //Add here all contributions
    //Exlcude the contribution form the pi_prim---pi_prim case e.g. last the last entry in the vec.
    //The first entry are the misid


    std::vector<double> LambdaHousehold2 = CalculateLambdaParam();
    std::vector<double> LambdaHousehold = CalculateLambdaParam2();

    CkDe_PiPi->AddContribution(1,LambdaHousehold[2],DLM_CkDecomposition::cFake);
    CkDe_PiPi->AddContribution(0,LambdaHousehold[1],DLM_CkDecomposition::cFeedDown);//Here Lambda from the Feeddown (the pT integrated one)
}
    //CkDe_PiPi->AddContribution(0,lambda_feed,DLM_CkDecomposition::cFeedDown);//Here Lambda form the Feeddown (guess the integrated one)
    //CkDe_PiPi->AddContribution(1,lambda_misid,DLM_CkDecomposition::cFake);//Here Lambda form the misid (here to maybe?

//First try of a TF1 object which gets then passed to CATS obj, for the Fitting procedure...
double MickeyFitter_pipi(double* x, double* par){
    double& MOM = *x;
    KITTY_CATS_FIT_PiPi->SetAnaSource(0,par[0],true);
    if(KITTY_CATS_FIT_PiPi->GetNumSourcePars()>1){
        KITTY_CATS_FIT_PiPi->SetAnaSource(1,par[1],true);
    }
    KITTY_CATS_FIT_PiPi->KillTheCat();
    CkDe_PiPi->Update(true);
    double CkVal;
    CkVal = CkDe_PiPi->EvalCk(MOM);
    double BaseLine = par[2]*(1.+par[3]*MOM+par[4]*MOM*MOM);
    return CkVal*BaseLine;
}
//FakeFitting: All params will be fixed, just to draw quickly the baselines
double MickeyFitter_pipi_Baseline(double* x, double* par){
    double& MOM = *x;
    double BaseLine = par[0]*(1.+par[1]*MOM+par[2]*MOM*MOM);
    return BaseLine;
}

//VariationSettings[0] = 0/1 default data only/all variations
//VariationSettings[1] = 0/1 bootstrap only/also systematics
//VariationSettings[2] = WhichIteration to perform, the 0th is the default
void Basics_PiPiCATS2(const bool& Identical, const bool& WithCoulomb, unsigned* VariationSettings){

    //make sure we have different seeds for each iteration. This way we can run the variations in parallel
    TRandom3 rangen(VariationSettings[2]+1);
    const bool DefaultVariation = VariationSettings[2]==0;

    // Set-up
    CATS PionKitty;

    const double sourcsize = 1.7;
    bool fUseResoMod = true;
    string SOURCE = "GaussCore";

    const double kMin = 0;
    const double kMax = 300;
    const double kWidth = 4;
    const unsigned NumMomBins = TMath::Nint( (kMax-kMin)/kWidth );

    const double kMinVars[] = {0, 8, 16};
    const double kMaxVars[] = {224, 248, 276};

    //(#bins,min,max) or (#bins,BinRangeArray[] as in ROOT)
    PionKitty.SetMomBins(NumMomBins,kMin,kMax);
    //object for the parameters to be used by the source function
    CATSparameters SOURCE_PARS(CATSparameters::tSource,1,true);
    //set the first and only par (source size)
    SOURCE_PARS.SetParameter(0,sourcsize); //a bit lower
    //say to CATS which Source function to use, and with which parameter set
    //SetUp for the SourceMod by Reso:
    DLM_CleverMcLevyResoTM CleverMcLevyResoTM;
    if (fUseResoMod) {
        if(SOURCE=="GaussCore") CleverMcLevyResoTM.InitStability(1,2-1e-6,2+1e-6);
        else CleverMcLevyResoTM.InitStability(21,1,2);
        CleverMcLevyResoTM.InitScale(38,0.15,2.0);
        CleverMcLevyResoTM.InitRad(257*2,0,64);
        CleverMcLevyResoTM.InitType(2);
        CleverMcLevyResoTM.SetUpReso(0,0.61);//fraction of resonances feeding into particle 1
        CleverMcLevyResoTM.SetUpReso(1,0.61);//fraction of resonances feeding into particle 2
        const double k_CutOff = 200;
        Float_t k_D;
        Float_t fP1;
        Float_t fP2;
        Float_t fM1;
        Float_t fM2;
        Float_t Tau1;
        Float_t Tau2;
        Float_t AngleRcP1;
        Float_t AngleRcP2;
        Float_t AngleP1P2;
        DLM_Random RanGen(11);//"DLM_Random.h"
        double RanVal1;
        double RanVal2;
        double RanVal3;
        //give the path to the file I will sent you
        //TFile* F_EposDisto_p_pReso = new TFile("../Files/ForMax_pi_piReso_withOmega.root");
        TFile* F_EposDisto_p_pReso = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/Max/ForMax_pi_piReso_noOmega.root");
        TNtuple* T_EposDisto_p_pReso = (TNtuple*)F_EposDisto_p_pReso->Get("InfoTuple_ClosePairs");
        unsigned N_EposDisto_p_pReso = T_EposDisto_p_pReso->GetEntries();
        T_EposDisto_p_pReso->SetBranchAddress("k_D",&k_D);
        T_EposDisto_p_pReso->SetBranchAddress("P1",&fP1);
        T_EposDisto_p_pReso->SetBranchAddress("P2",&fP2);
        T_EposDisto_p_pReso->SetBranchAddress("M1",&fM1);
        T_EposDisto_p_pReso->SetBranchAddress("M2",&fM2);
        T_EposDisto_p_pReso->SetBranchAddress("Tau1",&Tau1);
        T_EposDisto_p_pReso->SetBranchAddress("Tau2",&Tau2);
        T_EposDisto_p_pReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
        T_EposDisto_p_pReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
        T_EposDisto_p_pReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
        for(unsigned uEntry=0; uEntry<N_EposDisto_p_pReso; uEntry++){
            //by default Tau and Mass are read out from the file. We do NOT use that, instead we use the avg. cTau and Mass, which are set below
            T_EposDisto_p_pReso->GetEntry(uEntry);
            Tau1 = 0;
                //treat the omega separately
                if(fM2>782&&fM2<783){
                    fM2 = 782.6;
                    Tau2 = 23.24;
                }
                //the avg. values below should be computed for all resonances besides omega
                else{
                    Tau2 = 1.5;
                    fM2 = 1124;
                }
            if(k_D>k_CutOff) continue;
            RanVal1 = RanGen.Exponential(fM2/(fP2*Tau2));
            CleverMcLevyResoTM.AddBGT_PR(RanVal1,-cos(AngleRcP2));
            CleverMcLevyResoTM.AddBGT_RP(RanVal1,cos(AngleRcP2));
        }
        delete F_EposDisto_p_pReso;
        //TFile* F_EposDisto_pReso_pReso = new TFile("../Files/ForMax_piReso_piReso_withOmega.root");
        TFile* F_EposDisto_pReso_pReso = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/Max/ForMax_piReso_piReso_noOmega.root");
        TNtuple* T_EposDisto_pReso_pReso = (TNtuple*)F_EposDisto_pReso_pReso->Get("InfoTuple_ClosePairs");
        unsigned N_EposDisto_pReso_pReso = T_EposDisto_pReso_pReso->GetEntries();
        T_EposDisto_pReso_pReso->SetBranchAddress("k_D",&k_D);
        T_EposDisto_pReso_pReso->SetBranchAddress("P1",&fP1);
        T_EposDisto_pReso_pReso->SetBranchAddress("P2",&fP2);
        T_EposDisto_pReso_pReso->SetBranchAddress("M1",&fM1);
        T_EposDisto_pReso_pReso->SetBranchAddress("M2",&fM2);
        T_EposDisto_pReso_pReso->SetBranchAddress("Tau1",&Tau1);
        T_EposDisto_pReso_pReso->SetBranchAddress("Tau2",&Tau2);
        T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
        T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
        T_EposDisto_pReso_pReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
        for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_pReso; uEntry++){
            T_EposDisto_pReso_pReso->GetEntry(uEntry);
            //treat the omega separately
                if(fM1>782&&fM1<783){
                    fM1 = 782.6;
                    Tau1 = 23.24;
                }
                //the avg. values below should be computed for all resonances besides omega
                else{
                    Tau1 = 1.5;
                    fM1 = 1124;
                }
                //treat the omega separately
                if(fM2>782&&fM2<783){
                    fM2 = 782.6;
                    Tau2 = 23.24;
                }
                //the avg. values below should be computed for all resonances besides omega
                else{
                    Tau2 = 1.5;
                    fM2 = 1124;

                }
            if(k_D>k_CutOff) continue;
            RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
            RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
            CleverMcLevyResoTM.AddBGT_RR(RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));
        }
        delete F_EposDisto_pReso_pReso;
        //how many MC iterations to do to build up S(r) for a single set of source parameters
        if(SOURCE=="GaussCore") CleverMcLevyResoTM.InitNumMcIter(1000000);
        else CleverMcLevyResoTM.InitNumMcIter(100000);
        PionKitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyResoTM, 2);
        PionKitty.SetAnaSource(0,sourcsize);
        PionKitty.SetAnaSource(1,2.0);
        PionKitty.SetUseAnalyticSource(true);
        PionKitty.SetAutoNormSource(false);
    }
    else {
        //PionKitty.SetAnaSource(GaussSource,SOURCE_PARS); //Gaussian source
        PionKitty.SetAnaSource(CauchySource, SOURCE_PARS);	//Cauchy source
        PionKitty.SetUseAnalyticSource(true);
    }

    //standard settings for a CATS object which has no strong interaction potential included
    PionKitty.SetNumChannels(1);

    //which channel, how many PWs
    PionKitty.SetNumPW(0,0);
    //which channel, spin value
    PionKitty.SetSpin(0,0);
    //which channel, weight
    PionKitty.SetChannelWeight(0, 1);//In pipi only one!

    //include the coulomb interaction. Q1Q2 is the multiplied charge numbers of the two particles
    PionKitty.SetQ1Q2(1);
    //the PdgId is needed when using a source from a transport model
    //actually the only important thing here is if the particles are identical
    if(Identical) PionKitty.SetPdgId(211, 211);
    else PionKitty.SetPdgId(-211, 211);

    if(WithCoulomb) PionKitty.SetRedMass( 0.5*Mass_Pic );
    else if(Identical) PionKitty.SetRedMass( 0.5*Mass_Pi0 );
    else PionKitty.SetRedMass( (Mass_Pi0*Mass_Pic)/(Mass_Pi0+Mass_Pic) );

    PionKitty.KillTheCat();
    PionKitty.SetNotifications(CATS::nWarning);
    PionKitty.SetNotifications(CATS::nError);
    DLM_Ck pipiCk (PionKitty.GetNumSourcePars(),0,PionKitty); //behaves like a histo
    pipiCk.Update(kTRUE);

    //Set-up DLM_Dekomp
    TFile* FileROOT = new TFile("../Files/AnalysisResultsMC.root", "read");
    if (!FileROOT) {std::cout << "no file " << FileROOT<< std::endl;}
    TList *Results;
    auto* dir = FileROOT->GetDirectory(Form("%sResults%s", "MB", "0"));
    if (!dir) {std::cout << "no dir " << dir<< std::endl;}
    auto *dirResults = (TList*)dir->Get(Form("%sResults%s",  "MB", "0"));
    if (!dirResults) {std::cout << "no dirResults " << dirResults<< std::endl;}
    Results = (TList*)dirResults->FindObject("PairQA");
    if (!Results) {std::cout << "no Results " << Results<< std::endl;}
    TString FolderName = Form("QA_Particle%i_Particle%i", 0, 0);
    TList *PartList;
    PartList = (TList*) Results->FindObject(FolderName.Data());
    if (!PartList) {std::cout << "no PartList " << PartList<< std::endl;}
    TH2F* histo = (TH2F*) PartList->FindObject(Form("MomentumResolutionSE_Particle0_Particle0"));
    if (!histo) {std::cout << "no histo " << histo<< std::endl;}
    //Change of Units GeV -> MeV
    double Xmin = histo->GetXaxis()->GetXmin()*1000;
    double Xmax = histo->GetXaxis()->GetXmax()*1000;
    double Ymin = histo->GetYaxis()->GetXmin()*1000;
    double Ymax = histo->GetYaxis()->GetXmax()*1000;
    histo->GetXaxis()->SetLimits(Xmin, Xmax);
    histo->GetYaxis()->SetLimits(Ymin, Ymax);
    TString Name = histo->GetName();
    gROOT->cd();
    TH2F *histoCopy = (TH2F*)histo->Clone("histoCopy");
    delete FileROOT;
    histoCopy->SetName(Name);

    DLM_CkDecomposition CkDe_PiPiRealOne ("PipiDecomp", 2, pipiCk, histoCopy);
    SetTheLambdaParamsToDML_CkDecompObj(&CkDe_PiPiRealOne);
    CkDe_PiPiRealOne.Update(kTRUE);
    KITTY_CATS_FIT_PiPi = &PionKitty;
    CkDe_PiPi= &CkDe_PiPiRealOne;
    KITTY_CATS_FIT_PiPi->KillTheCat();
    CkDe_PiPi->Update(kTRUE);
    //TFile* fInput = new TFile("../Files/CFOutput_pipi_MB_mult_2_kT_2.root","read");//Put here the Input from the GentleFEMTO but remember in GeV and as a TGraph? But here is a TF1*
    //TH1F* hInput = (TH1F*)fInput->Get("hCkTotNormWeightpipiMeV"); //Can use the structure for the plotting macro here!!!
 	std::vector<TH1D*> fGentleCF;
 	std::vector<TH1D*> fGentleCF_Rebinned;
    std::vector<TH1D*> fGentleBuffer;
    std::vector<TF1*> fFITTERBuffer;
    std::vector<TF1*> fFITTERBASEBuffer;
    std::vector<TGraph*> fgFITTERBuffer;
    std::vector<TGraph*> fgFITTERBASEBuffer;
    std::vector<int> MultBins;
 	MultBins.push_back(0);
  	MultBins.push_back(18);
  	MultBins.push_back(30);
  	MultBins.push_back(999);
    std::vector<double> kTBins;
    kTBins.push_back(0.15);
    kTBins.push_back(0.3);
    kTBins.push_back(0.5);
    kTBins.push_back(0.7);
    kTBins.push_back(0.9);
    kTBins.push_back(1.5);
    int iVAR = -1;
    if(VariationSettings[0]) iVAR = rangen.Integer(10);
    if(VariationSettings[0]&&VariationSettings[2]) iVAR = 0;
    for (unsigned int i = 0; i < MultBins.size() - 1; ++i) {
        TString namePart1 = Form("mult:%i-%i_",MultBins[i],MultBins[i+1]);
        for (unsigned int j = 0; j < kTBins.size() - 1; ++j) {
            //fGentleBuffer = ExtractRelevantPlotsGF(Form("../Files/CFOutput_ForDimi_pipiVar%i_MB_mult_%i_kT_%i.root",iVAR,i+1,j+1),j,i,iVAR);
            if(iVAR>=0) fGentleBuffer = ExtractRelevantPlotsGF(Form("../Files/AOD/CFOutput_pipiVar%i_MB_mult_%i_kT_%i.root",iVAR,i+1,j+1),j,i,iVAR);
            else fGentleBuffer = ExtractRelevantPlotsGF(Form("../Files/AOD/CFOutput_pipi_MB_mult_%i_kT_%i.root",i+1,j+1),j,i,iVAR);
            TH1D* dummy = fGentleBuffer[0];
            dummy->SetTitle(Form("#splitline{%i < mult < %i}{%0.2f < k_{T} #leq %.2f [GeV/#it{c}]}",MultBins[i],MultBins[i+1]+1,kTBins[j],kTBins[j+1]));
            dummy->SetName(TString::Format("hFit_Mult%u_kT%u",i,j));
            fGentleCF.emplace_back(dummy);
            fGentleCF_Rebinned.emplace_back(fGentleBuffer[1]);
        }
    }
    int Index = 0;
    for ( auto hfit : fGentleCF) {
        //TH1D* hFit = (TH1D*)it->Clone(TString::Format("hFit_Mult%u_kT%u",i,j));

        //bootstrap here
        float Value;
        float Error;
        float RanVal;
        for(unsigned uBin=0; uBin<hfit->GetNbinsX(); uBin++){
            if(DefaultVariation) break;
            Value = hfit->GetBinContent(uBin+1);
            Error = hfit->GetBinError(uBin+1);
            RanVal = rangen.Gaus(Value,Error);
            hfit->SetBinContent(uBin+1,RanVal);
            hfit->SetBinError(uBin+1,Error);
        }
        Index++;
		std::cout<<"At Index: "<<Index<<"\n";
		string FitterName = Form("FITTER_%i",Index);
		string FitterNameBase = Form("FITTER_%i_Baseline",Index);
		string gFitterName = Form("gFITTER_%i",Index);
		string gFitterNameBase = Form("gFITTER_%i_Baseline",Index);
		double kFitMin = kMinVars[rangen.Uniform(3)];
		double kFitMax = kMaxVars[rangen.Uniform(3)];
		if(DefaultVariation){
            kFitMin = kMinVars[1];
            kFitMax = kMaxVars[1];
		}
		TF1* FITTER = new TF1(FitterName.c_str(),MickeyFitter_pipi,kFitMin,kFitMax,5);
		FITTER->SetParameter(0,sourcsize);FITTER->SetParLimits(0,0.5,3);
		FITTER->FixParameter(1,2);
		FITTER->SetParameter(2,1);
		FITTER->SetParameter(3,0);FITTER->SetParLimits(3,-1e-2,1e-2);
		FITTER->SetParameter(4,0);FITTER->SetParLimits(4,-1e-3,1e-3);
		/*FITTER->FixParameter(0,sourcsize);
		FITTER->FixParameter(1,2);
		FITTER->FixParameter(2,1);
		FITTER->FixParameter(3,0);
		FITTER->FixParameter(4,0);*/
		string ItName =  Form("%s_%i",hfit->GetName(),Index);
		//it->SetName(ItName.c_str());
		hfit->Fit(FITTER,"Q, S, R, M");
        TF1* FITTERBASE = new TF1(FitterNameBase.c_str(),MickeyFitter_pipi_Baseline,kFitMin,kFitMax,3);
		FITTERBASE->FixParameter(0,FITTER->GetParameter(2));
		FITTERBASE->FixParameter(1,FITTER->GetParameter(3));
		FITTERBASE->FixParameter(2,FITTER->GetParameter(4));
		fFITTERBuffer.emplace_back(FITTER);
		fFITTERBASEBuffer.emplace_back(FITTERBASE);
		TGraph* gFITTER = new TGraph();
		gFITTER->SetName(gFitterName.c_str());
		TGraph* gFITTERBASE = new TGraph();
		gFITTERBASE->SetName(gFitterNameBase.c_str());
        unsigned uPoint=0;
		for(double kVal=kWidth*0.5; kVal<=kMax; kVal+=kWidth){
            gFITTER->SetPoint(uPoint,kVal,FITTER->Eval(kVal));
            gFITTERBASE->SetPoint(uPoint,kVal,FITTERBASE->Eval(kVal));
            uPoint++;
		}
		fgFITTERBuffer.emplace_back(gFITTER);
		fgFITTERBASEBuffer.emplace_back(gFITTERBASE);
	}
	TFile* OutputFile = new TFile(TString::Format("../Output/OUTPUT_%u_%u_%u.root",VariationSettings[0],VariationSettings[1],VariationSettings[2]),"recreate");
    int Index2 = 0;
    std::vector<float>  ParYValueSource;
    std::vector<float>  ParYErrorSource;
    std::vector<float>  ParXValueSource;
    std::vector<float>  ParXErrorSource;
    for ( auto it : fGentleCF) {
        fFITTERBASEBuffer[Index2]->Write(fFITTERBASEBuffer[Index2]->GetName());
        fFITTERBuffer[Index2]->Write(fFITTERBuffer[Index2]->GetName());
        fgFITTERBASEBuffer[Index2]->Write(fgFITTERBASEBuffer[Index2]->GetName());
        fgFITTERBuffer[Index2]->Write(fgFITTERBuffer[Index2]->GetName());
        ParYValueSource.push_back(fFITTERBuffer[Index2]->GetParameter(0));
	    ParYErrorSource.push_back(fFITTERBuffer[Index2]->GetParError(0));
	    ParXValueSource.push_back(Index2);
	    ParXErrorSource.push_back(0.);

	    string Cname = Form("%i",Index2);
        auto* C = new TCanvas(Cname.c_str(),Cname.c_str());
        gPad->SetTicks();
        gStyle->SetTitleY(0.9);
        gPad->Update();
	    if (true) {
            it->GetXaxis()->SetTitle("k* [MeV/#it{c}]");
            it->GetYaxis()->SetTitle("C(k*)");
        }
        else {
            it->GetXaxis()->SetTitle("#it{p}_{T} [GeV/#it{c}]");
            it->GetYaxis()->SetTitle("Entries");
	    }
        it->GetYaxis()->SetRangeUser(0.6, 1.8);
        it->GetXaxis()->SetRangeUser(kMin, kMax);
        it->SetTitle(Form("#splitline{%s}{<r> = %f #pm %f}",it->GetTitle(),fFITTERBuffer[Index2]->GetParameter(0),fFITTERBuffer[Index2]->GetParError(0)));
	    it->Draw("pe");
	    //it->SetTitle(Titlebuffer);
	    //fFITTERBuffer[Index2]->Draw("lsame");
	    TF1* f0 = fFITTERBASEBuffer[Index2];
	    f0->SetLineColor(kBlue);
	    f0->Draw("same");
        auto legend2 = new TLegend(0.51,0.52,0.89,0.73);
        legend2->AddEntry(it,"CF","lep");
	    legend2->AddEntry(fFITTERBuffer[Index2],"Full Fit","l");
	    legend2->AddEntry(f0,"Baseline","l");
   	    legend2->SetBorderSize(0);
	    legend2->Draw("Same");
	    //C->Write(it->GetName());
	    string name = (it)->GetName();
        auto savenamepart="Plot_withOmegaInLambda400/"+name;
	    auto savename=savenamepart+".pdf";
        C->Print(savename.c_str());
	    it->Write(it->GetName());
        Index2++;
	}
    //printf("I am here\n");
    auto c1 = new TCanvas("caa1","A Simple Graph with error bars",200,10,700,500);
    c1->SetFillColor(42);
    c1->SetGrid();
    c1->GetFrame()->SetFillColor(21);
    c1->GetFrame()->SetBorderSize(12);
    std::vector<float>  ParYValueSourceSplit;
    std::vector<float>  ParYErrorSourceSplit;
    std::vector<float>  ParXValueSourceSplit;
    std::vector<float>  ParXErrorSourceSplit;
    int iter5=0;
    int track=1;
    //printf("I am here\n");
    for (unsigned int i = 0; i < MultBins.size() - 1; ++i) {
    for (unsigned int j = 0; j < kTBins.size() - 1; ++j) {
        //printf("I am here1\n");
            ParXValueSourceSplit.push_back(j);
       //printf("I am here2\n");
            ParYValueSourceSplit.push_back(ParYValueSource[iter5]);
       //printf("I am here3\n");
            ParXErrorSourceSplit.push_back(0);
       //printf("I am here4\n");
            ParYErrorSourceSplit.push_back(ParYErrorSource[iter5]);
       //printf("I am here5\n");
        //   printf("ParXValueSourceSplit[j] %f",ParXValueSourceSplit[j]);
          // printf("ParXErrorSourceSplit[j] %f",ParXErrorSourceSplit[j]);
           //printf("ParYValueSourceSplit[j] %f",ParYValueSourceSplit[j]);
           //printf("ParYErrorSourceSplit[j] %f",ParYErrorSourceSplit[j]);
            iter5++;
       //printf("I am here\n");
        }
        auto gr = new TGraphErrors(5,&ParXValueSourceSplit[0],&ParYValueSourceSplit[0],&ParXErrorSourceSplit[0],&ParYErrorSourceSplit[0]);
       //printf("I am here\n");
        gr->SetTitle("TGraphErrors Example");
        gr->SetMarkerColor(4);
        gr->SetMarkerStyle(21);
        gr->Draw("ALP");
        string grName = Form("FitRadiusValues_%i",track);
        gr->Write(grName.c_str());
        ParXValueSourceSplit.clear();
        ParYValueSourceSplit.clear();
        ParYErrorSourceSplit.clear();
        ParXErrorSourceSplit.clear();
       //printf("I am here\n");
        track++;
    }
    //printf("I am here\n");
    OutputFile->Write();
    //delete fFITTERBuffer;
    //delete fGentleCF;
    //delete OutputFile;


   // TGraph* grCk = new TGraph();
   // grCk->SetName("grCk");
   // grCk->Set(NumMomBins);
   // for(unsigned uMom=0; uMom<NumMomBins; uMom++){
   //     grCk->SetPoint(uMom,PionKitty.GetMomentum(uMom),PionKitty.GetCorrFun(uMom));
   // }

   // return grCk;

   //delete FITTER;
   //delete FITTERBASE;
   //delete C;
   //delete legend2;
   //delete c1;
   //delete gr;
   //delete FileROOT;
   //delete OutputFile;

}

void PlotTheResults(unsigned* VariationSettings){

    double kMin = 0;
    double kMax = 252;
    double kWidth = 4;
    unsigned kBins = TMath::Nint((kMax-kMin)/(kWidth));

    const unsigned NumMultBins = 3;
    const unsigned NumKtBins = 5;

    const unsigned NumVariations = VariationSettings[2];
    TString BaseFileName = TString::Format("OUTPUT_%u_%u_",VariationSettings[0],VariationSettings[1]);

    TH2F*** hCkFit = new TH2F** [NumMultBins];
    TH2F*** hCkBl = new TH2F** [NumMultBins];
    TH1F*** hRadius = new TH1F** [NumMultBins];
    for(unsigned uMult=0; uMult<NumMultBins; uMult++){
        hCkFit[uMult] = new TH2F* [NumKtBins];
        hCkBl[uMult] = new TH2F* [NumKtBins];
        hRadius[uMult] = new TH1F* [NumKtBins];
        for(unsigned uKt=0; uKt<NumKtBins; uKt++){
            hCkFit[uMult][uKt] = new TH2F(TString::Format("hCkFit_Mult%u_kT%u",uMult,uKt),TString::Format("hCkFit_Mult%u_kT%u",uMult,uKt),kBins,kMin,kMax,4096,0.5,2.5);
            hCkBl[uMult][uKt] = new TH2F(TString::Format("hCkBl_Mult%u_kT%u",uMult,uKt),TString::Format("hCkBl_Mult%u_kT%u",uMult,uKt),kBins,kMin,kMax,4096,0.5,2.5);
            hRadius[uMult][uKt] = new TH1F(TString::Format("hRadius_Mult%u_kT%u",uMult,uKt),TString::Format("hRadius_Mult%u_kT%u",uMult,uKt),4096*2,0.5,3.0);
        }
    }

    for(unsigned uVar=0; uVar<NumVariations; uVar++){
        TString InputFileName = BaseFileName+TString::Format("%u.root",uVar);
        TFile InputFile(InputFileName,"read");
        for(unsigned uMult=0; uMult<NumMultBins; uMult++){
            TGraph* gRadii = (TGraph*)InputFile.Get(TString::Format("FitRadiusValues_%u",uMult+1));
            for(unsigned uKt=0; uKt<NumKtBins; uKt++){
                TString InputFitName = TString::Format("gFITTER_%u",uKt+uMult*NumKtBins+1);
                TString InputBlName = TString::Format("gFITTER_%u_Baseline",uKt+uMult*NumKtBins+1);
                TGraph* gFITTER = (TGraph*)InputFile.Get(InputFitName);
                TGraph* gFITTERBASE = (TGraph*)InputFile.Get(InputBlName);
                for(unsigned uBin=0; uBin<kBins; uBin++){
                    double Momentum = hCkFit[uMult][uKt]->GetBinCenter(uBin+1);
                    hCkFit[uMult][uKt]->Fill(Momentum,gFITTER->Eval(Momentum));
                    hCkBl[uMult][uKt]->Fill(Momentum,gFITTERBASE->Eval(Momentum));
                }
                double RADIUS;
                double DUMMY;
                gRadii->GetPoint(uKt,DUMMY,RADIUS);
                hRadius[uMult][uKt]->Fill(RADIUS);
            }
        }
    }

}


std::vector<double> CalculateLambdaParam(double pur1, double pur2, bool SamePart, bool extinput, bool verbose) {
/*The average including under and overflow of the purity api is: 0.925035
The average including under and overflow of the purity pi is: 0.923983
For pi+ taken from the archvied version of fitDCA_FractionFitter_SigmaK0 24.03.2020
Renorm factor = 0.996427
 frac pri  = 0.937751
 frac secLam  = 0.00652577
 frac secSig  = 0.0019332
 frac secSigM  = 0.00321461
 frac secKS  = 0.0361176
 frac secKch  = 0.014458
 tot sum fractions  = 1
 frac misid  = -1.03972e-07*/
	double checksum = 0;
	int multiplicator = 0;
	//Structure of fractions = (prim, second*[size-2], misid)
	std::vector<double> fractions;
	if(!extinput) {
		fractions.push_back(0.937751*(1-0.122048));	//pri //withomega
		//fractions.push_back(0.937751*(1-0.075493));	//pri
		fractions.push_back(0.00652577);//Lam
		fractions.push_back(0.014458);	//Kch
		fractions.push_back(0.00321461);//SM
		fractions.push_back(0.0019332); //SP
		fractions.push_back(0.0361176); //KS
                // fractions.push_back(0.937751*0.075493);  //prim-stronglonglived
		fractions.push_back(0.937751*(1-0.122048));	//pri //withomega
		fractions.push_back(-1.03972e-07);//misid
	}
	double OneOverNorm=0;
	if(SamePart) pur2=pur1;
	//reverse vector
	std::reverse(fractions.begin(), fractions.end());
	for( int unsigned iter3=0; iter3<fractions.size(); iter3++) {
		if (verbose) printf("fractions[%u] should be reversed: %f \n", iter3, fractions[iter3]);
	}
	//fill misid-misid case actuall also contains all other cases prim/weak-misid and the other way around
	std::vector<double> LambdaParam;
	LambdaParam.push_back((pur1)*fractions[0]*(1.-pur2)*fractions[0]);
	//initialise normalisation
	//Skip misid
        for(unsigned int iter=0; iter<fractions.size(); iter++) {
		if (verbose) printf("iter: %u \n",iter);
		for (unsigned int iter2=1; iter2 < iter+1; iter2++) {
			(iter +1 == fractions.size() && iter2 != iter) ? multiplicator=2 : multiplicator=1; // take care of the case pi_prim---pi_weak AND pi_weak---pi_prim hence *2 in this case...
  			double Lambda = pur1*fractions[iter]*pur2*fractions[iter2]*multiplicator;
			if (verbose) printf("iter: %u iter2: %u Lambda %f multiplicator: %i \n",iter,iter2,Lambda,multiplicator);
			LambdaParam.push_back(Lambda);
			OneOverNorm += Lambda; //doesn't care of the misid
		}
	}
        OneOverNorm += LambdaParam[0];//take care of the misid
	auto Norm = 1/OneOverNorm;
	//Do the normalization
	for( int unsigned iter3=0; iter3<LambdaParam.size(); iter3++) {
		LambdaParam[iter3] *= Norm;
		checksum += LambdaParam[iter3];
		if (verbose) printf("OLDLambdaParam[%u] : %f checksum : %f \n", iter3, LambdaParam[iter3],checksum);
	}
	if ( fabs( (checksum-1.)/(checksum+1.)/2. )<1e-6 ) {//float comparison with presicion of 1e-6
		return LambdaParam;
	} else {
		printf("Check the function the sum of LambdaParam is unequal to 1: %f \nProceed with CAUTION!! \n",checksum);
		return LambdaParam;
	}
}


std::vector<double> CalculateLambdaParam2(double pur1, double pur2, bool SamePart, bool extinput, bool verbose) {

	double checksum = 0;
	//Structure of fractions = (prim, second*[size-2], misid)
	std::vector<double> fractions;
	if(!extinput) {
		fractions.push_back(0.937751*(1-0.122048));	//pri //withomega and everything else
		//fractions.push_back(0.937751*(1-0.075493));	//pri //withoutomega
		fractions.push_back(0.00652577);//Lam
		fractions.push_back(0.014458);	//Kch
		fractions.push_back(0.00321461);//SM
		fractions.push_back(0.0019332); //SP
		fractions.push_back(0.0361176); //KS
                //fractions.push_back(0.937751*0.075493);  //prim-stronglonglived
                fractions.push_back(0.937751*0.122048);  //prim-stronglonglived
		fractions.push_back(-1.03972e-07);//misid
	}
	if(SamePart) pur2=pur1;
	float OneOverNorm=0;
        printf("pur2 %f pur1 %f\n",pur2,pur1);
        //fill misid-misid case actuall also contains all other cases prim/weak-misid and the other way around
	std::vector<double> LambdaParam;
	float Lambdaprim = (pur1)*fractions[0]*(pur2)*fractions[0]; //primaries
	printf("Lambdaprim : %f \n", Lambdaprim);
	float Lambdamisid = (1.-pur1)*((pur2)*2+(1.-pur1)); //misid
	printf("Lambdamisid : %f \n", Lambdamisid);
        float Lambdaother = 1. - Lambdaprim - Lambdamisid;
	printf("Lambdaother : %f \n", Lambdaother);


        //OneOverNorm += Lambdaprim;
        //OneOverNorm += Lambdamisid;
	//OneOverNorm += Lambdaother;
	//auto Norm = 1/OneOverNorm;
	//Lambdaprim *=Norm;
	//Lambdamisid *=Norm;
	//Lambdaother *=Norm;
	LambdaParam.push_back(Lambdaprim);
	LambdaParam.push_back(Lambdaother);
	LambdaParam.push_back(Lambdamisid);
	for( int unsigned iter3=0; iter3<LambdaParam.size(); iter3++) {
		checksum += LambdaParam[iter3];
		if (verbose) printf("NEWLambdaParam[%u] : %f checksum : %f \n", iter3, LambdaParam[iter3],checksum);
	}
	if ( fabs( (checksum-1.)/(checksum+1.)/2. )<1e-6 ) {//float comparison with presicion of 1e-6
		return LambdaParam;
	} else {
		printf("Check the function the sum of LambdaParam is unequal to 1: %f \nProceed with CAUTION!! \n",checksum);
		return LambdaParam;
	}
}



std::vector<TH1D*> ExtractRelevantPlotsGF(string GentleFemto, int k, int m, int iVAR) {
	std::vector<TH1D*> fResults;
	TFile* file = TFile::Open(GentleFemto.c_str());
    TString sVAR = "";
    if(iVAR>=0) sVAR = TString::Format("Var%i",iVAR);

    TH1D* CF;
    if(iVAR>=0) CF = (TH1D*)file->FindObjectAny(TString::Format("hCkTotNormWeightForDimi_pipi%sMeV",sVAR.Data()));
    else CF = (TH1D*)file->FindObjectAny(TString::Format("hCkTotNormWeightpipi%sMeV",sVAR.Data()));

    TH1D* CF_Rebinned;
    if(iVAR>=0) CF_Rebinned = (TH1D*)file->FindObjectAny(TString::Format("hCk_RebinnedForDimi_pipi%sMeV_0",sVAR.Data()));
    else CF_Rebinned = (TH1D*)file->FindObjectAny(TString::Format("hCk_Rebinnedpipi%sMeV_0",sVAR.Data()));
	fResults.emplace_back(CF);
	fResults.emplace_back(CF_Rebinned);
        for ( auto it : fResults ) {
           it->SetDirectory(0);
        }
	file->Close();

	return fResults;
}



void Compare_AOD_nano(){

    int NumMultBins = 3;
    int NumKtBins = 5;

    TString FileNameAOD;
    TString FileNameNano;

    std::vector<TH1D*> fGentleBufferAod;
    std::vector<TH1D*> fGentleBufferNano;

    TFile OutputFile("../Output/Compare_AOD_nano.root","recreate");

    TH1D*** hRatio_AOD_nano = new TH1D** [NumMultBins];
    TH1D*** hAOD = new TH1D** [NumMultBins];
    TH1D*** hNano = new TH1D** [NumMultBins];

    for(unsigned uMult=0; uMult<NumMultBins; uMult++){
        hRatio_AOD_nano[uMult] = new TH1D* [NumKtBins];
        hAOD[uMult] = new TH1D* [NumKtBins];
        hNano[uMult] = new TH1D* [NumKtBins];
        for(unsigned ukT=0; ukT<NumKtBins; ukT++){
            fGentleBufferAod = ExtractRelevantPlotsGF(Form("../Files/AOD/CFOutput_pipi_MB_mult_%i_kT_%i.root",uMult+1,ukT+1),ukT,uMult,-1);
            TH1D* haod = fGentleBufferAod[0];
            fGentleBufferNano = ExtractRelevantPlotsGF(Form("../Files/CFOutput_ForDimi_pipiVar%i_MB_mult_%i_kT_%i.root",0,uMult+1,ukT+1),ukT,uMult,0);
            TH1D* hnano = fGentleBufferNano[0];

            OutputFile.cd();
            hAOD[uMult][ukT] = (TH1D*)haod->Clone(TString::Format("hAod_Mult%u_kT%u",uMult,ukT));
            hNano[uMult][ukT] = (TH1D*)hnano->Clone(TString::Format("hNano_Mult%u_kT%u",uMult,ukT));
            hRatio_AOD_nano[uMult][ukT] = (TH1D*)haod->Clone(TString::Format("hRatio_Mult%u_kT%u",uMult,ukT));
            hRatio_AOD_nano[uMult][ukT]->Divide(hNano[uMult][ukT]);

            OutputFile.cd();
            hAOD[uMult][ukT]->Write();
            hNano[uMult][ukT]->Write();
            hRatio_AOD_nano[uMult][ukT]->Write();
        }
    }


}



