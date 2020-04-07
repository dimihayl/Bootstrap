#include "Basics_boot.h"
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
#include "TMath.h"


//Outside definition in order for the following to work
    DLM_CkDecomposition* CkDe_PiPi_boot;
    CATS* KITTY_CATS_FIT_PiPi_boot;

void SetTheLambdaParamsToDML_CkDecompObj_boot (DLM_CkDecomposition* CkDe_PiPi_boot, bool verbose, int iVar, TRandom3 SEED) {
    //Generate all the Lambdaparams
    //Add here all contributions
    //Exlcude the contribution form the pi_prim---pi_prim case e.g. last the last entry in the vec.
    //The first entry are the misid

    std::vector<double> LambdaHousehold = CalculateLambdaParam_boot();
	
    if ( iVar >= 0 ) {
	printf("Mutate LambdaFeed for iVar: %i.\n",iVar);
  	double mutatedFeed = SEED.Uniform(LambdaHousehold[1]*0.8,LambdaHousehold[1]*1.2);
	LambdaHousehold.at(1) = mutatedFeed;
    }

    CkDe_PiPi_boot->AddContribution(1,LambdaHousehold[2],DLM_CkDecomposition::cFake);
    CkDe_PiPi_boot->AddContribution(0,LambdaHousehold[1],DLM_CkDecomposition::cFeedDown);//Here Lambda from the Feeddown (the pT integrated one)
}


//First try of a TF1 object which gets then passed to CATS obj, for the Fitting procedure...
double MickeyFitter_pipi_boot(double* x, double* par){
    double& MOM = *x;
    KITTY_CATS_FIT_PiPi_boot->SetAnaSource(0,par[0],true);
    if(KITTY_CATS_FIT_PiPi_boot->GetNumSourcePars()>1){
        KITTY_CATS_FIT_PiPi_boot->SetAnaSource(1,par[1],true);
    }
    KITTY_CATS_FIT_PiPi_boot->KillTheCat();
    CkDe_PiPi_boot->Update(true);
    double CkVal;
    CkVal = CkDe_PiPi_boot->EvalCk(MOM);
    double BaseLine = par[2]*(1.+par[3]*MOM+par[4]*MOM*MOM);
    return CkVal*BaseLine;
}
//FakeFitting: All params will be fixed, just to draw quickly the baselines
double MickeyFitter_pipi_boot_Baseline(double* x, double* par){
    double& MOM = *x;
    double BaseLine = par[0]*(1.+par[1]*MOM+par[2]*MOM*MOM);
    return BaseLine;
}

//VariationSettings[0] = 0/1 default data only/all variations
//VariationSettings[1] = 0/1 bootstrap only/also systematics
//VariationSettings[2] = WhichIteration to perform, the 0th is the default
void Basics_PiPiCATS_boot(const bool& Identical, const bool& WithCoulomb, unsigned* VariationSettings){
    //Trace the function calls and construct an abbort condition per the desired iteration maximum. What is here usually used? O(50)???
    double call_count = 0.;
    double call_count_pol1 = 0.; 
    double call_count_pol2 = 0.;
    const int IterMax = 100;

    //Set-up the arrays:
    const double kMin = 0;
    const double kMax = 300;
    const double kWidth = 4;
    const unsigned NumMomBins = TMath::Nint( (kMax-kMin)/kWidth );

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
   
    auto NmultBins = MultBins.size()-1;
    auto NktBins = MultBins.size()-1;
    auto TotalArraySize = NmultBins*NktBins;

    double** fEachBinPol1Mean = new double*[TotalArraySize];
    double** fEachBinPol2Mean = new double*[TotalArraySize];
    double** fEachBinPolMixedMean = new double*[TotalArraySize];
    double** fEachBinPol1BaselineMean = new double*[TotalArraySize];
    double** fEachBinPol2BaselineMean = new double*[TotalArraySize];
    double** fEachBinPolMixedBaselineMean = new double*[TotalArraySize];
    double** fEachBinPol1Stdev = new double*[TotalArraySize];
    double** fEachBinPol2Stdev = new double*[TotalArraySize];
    double** fEachBinPolMixedStdev = new double*[TotalArraySize];
    double** fEachBinPol1BaselineStdev = new double*[TotalArraySize];
    double** fEachBinPol2BaselineStdev = new double*[TotalArraySize];
    double** fEachBinPolMixedBaselineStdev = new double*[TotalArraySize];
    for (int i = 0; i < NumMomBins; ++i){
       fEachBinPol1Mean[i] = new double[NumMomBins];
       fEachBinPol2Mean[i] = new double[NumMomBins];
       fEachBinPolMixedMean[i] = new double[NumMomBins];
       fEachBinPol1BaselineMean[i] = new double[NumMomBins];
       fEachBinPol2BaselineMean[i] = new double[NumMomBins];
       fEachBinPolMixedBaselineMean[i] = new double[NumMomBins];
       fEachBinPol1Stdev[i] = new double[NumMomBins]; //Use for post processing
       fEachBinPol2Stdev[i] = new double[NumMomBins];
       fEachBinPolMixedStdev[i] = new double[NumMomBins];
       fEachBinPol1BaselineStdev[i] = new double[NumMomBins];
       fEachBinPol2BaselineStdev[i] = new double[NumMomBins];
       fEachBinPolMixedBaselineStdev[i] = new double[NumMomBins];
    }

    TGraphErrors* fResultPol1 = new TGraphErrors[TotalArraySize];
    TGraphErrors* fResultPol2 = new TGraphErrors[TotalArraySize];
    TGraphErrors* fResultPolMixed = new TGraphErrors[TotalArraySize];
    TGraphErrors* fResultPol1Baseline = new TGraphErrors[TotalArraySize];
    TGraphErrors* fResultPol2Baseline = new TGraphErrors[TotalArraySize];
    TGraphErrors* fResultPolMixedBaseline = new TGraphErrors[TotalArraySize];

    while(call_count < IterMax) {
    //trace the progess
    float perToBeDone = (IterMax-call_count)/IterMax;
    printf("%.0f%% percent till completion of the task.\n",perToBeDone*100);
    //create the uniqueIdentifier for the Iteration run.
    string uniqueIdent = Form("_Iter_%.0f",call_count);
    //make sure we have different seeds for each iteration. This way we can run the variations in parallel since here TRandom3 is initalized do I need to put here the WHILE or if construction?
    //TRandom3 rangen(VariationSettings[2]+1); // leave in case the function is re-structured to run really in parallel...
    TRandom3 rangen(call_count);
    //const bool DefaultVariation = VariationSettings[2]==0;
    const bool DefaultVariation = call_count==0;
    call_count++;
    //check the configuration:
    printf("VariationSettings[0]: %u.\n",VariationSettings[0]);
    printf("VariationSettings[1]: %u.\n",VariationSettings[1]);
    printf("VariationSettings[2]: %u.\n",VariationSettings[2]);

    // Set-up
    CATS PionKitty;

    const double sourcsize = 1.7;
    bool fUseResoMod = true;
    string SOURCE = "GaussCore";

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
        TFile* F_EposDisto_p_pReso = new TFile("../Files/ForMax_pi_piReso_noOmega.root");
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
        TFile* F_EposDisto_pReso_pReso = new TFile("../Files/ForMax_piReso_piReso_noOmega.root");
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

    if(WithCoulomb) PionKitty.SetRedMass( 0.5*Mass_Pich );
    else if(Identical) PionKitty.SetRedMass( 0.5*Mass_Pi0b );
    else PionKitty.SetRedMass( (Mass_Pi0b*Mass_Pich)/(Mass_Pi0b+Mass_Pich) );

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
    DLM_CkDecomposition CkDe_PiPi_bootRealOne ("PipiDecomp", 2, pipiCk, histoCopy);
    //initialise some counter
    //Put here some if statement or While Condition
    int iVAR = -1;
    if(VariationSettings[0]) {
	iVAR = rangen.Uniform(10); //change asap
	printf("UseOnlyDefault Var is set to false, busy with iVAR: %i.\n",iVAR);
    }
    if(VariationSettings[0]&&VariationSettings[2]) {
	iVAR = 0;
	printf("UseOnlyDefault is set to false, busy with iVAR: %i.\n",iVAR);
    }
    printf("Busy with iVAR: %i.\n",iVAR);
    SetTheLambdaParamsToDML_CkDecompObj_boot(&CkDe_PiPi_bootRealOne, false, iVAR, rangen);
    CkDe_PiPi_bootRealOne.Update(kTRUE);
    KITTY_CATS_FIT_PiPi_boot = &PionKitty;
    CkDe_PiPi_boot= &CkDe_PiPi_bootRealOne;
    KITTY_CATS_FIT_PiPi_boot->KillTheCat();
    CkDe_PiPi_boot->Update(kTRUE);
    std::vector<TH1D*> fGentleCF;
    std::vector<TH1D*> fGentleCF_Rebinned;
    std::vector<TH1D*> fGentleBuffer;
    std::vector<TF1*> fFITTERBuffer;
    std::vector<TF1*> fFITTERBASEBuffer;
    std::vector<TGraph*> fgFITTERBuffer;
    std::vector<TGraph*> fgFITTERBASEBuffer;
    for (unsigned int i = 0; i < MultBins.size() - 1; ++i) {
        TString namePart1 = Form("mult:%i-%i_",MultBins[i],MultBins[i+1]);
        for (unsigned int j = 0; j < kTBins.size() - 1; ++j) {
            //fGentleBuffer = ExtractRelevantPlotsGF(Form("../Files/CFOutput_ForDimi_pipiVar%i_MB_mult_%i_kT_%i.root",iVAR,i+1,j+1),j,i,iVAR);
            if(iVAR>=0) fGentleBuffer = ExtractRelevantPlotsGF_boot(Form("../Files/NanoAOD/CFOutput_pipiVar%i_MB_mult_%i_kT_%i.root",iVAR,i+1,j+1),j,i,iVAR);
            else fGentleBuffer = ExtractRelevantPlotsGF_boot(Form("../Files/AOD/CFOutput_pipi_MB_mult_%i_kT_%i.root",i+1,j+1),j,i,iVAR);
            TH1D* dummy = fGentleBuffer[0];
            dummy->SetTitle(Form("#splitline{%i < mult < %i}{%0.2f < k_{T} #leq %.2f [GeV/#it{c}]}",MultBins[i],MultBins[i+1]+1,kTBins[j],kTBins[j+1]));
            dummy->SetName(TString::Format("hFit_Mult%u_kT%u",i,j));
            fGentleCF.emplace_back(dummy);
            fGentleCF_Rebinned.emplace_back(fGentleBuffer[1]);
        }
    }
    int Index = 0;
    //also int indx3 = rangen.Integer(2); and then if statement to switch bool isPOL1 and take care of the parameters
    bool usePol1 = false;
    if  (DefaultVariation) {
	usePol1 = false;
    } else {
    	int indx3 = rangen.Integer(2);
    	bool usePol1 = (indx3) ? true : false;
    }
    TString AddendumForPolType = (usePol1) ? "pol1" : "pol2";
    printf("AddendumForPolType: %s.\n",AddendumForPolType.Data()); 	
    for ( auto hfit : fGentleCF) {
        //TH1D* hFit = (TH1D*)it->Clone(TString::Format("hFit_Mult%u_kT%u",i,j));
        //bootstrap here
        float Value;
        float Error;
        float RanVal;
        for(unsigned uBin=0; uBin<hfit->GetNbinsX(); uBin++){
            if(DefaultVariation || VariationSettings[1] == 0) break;
            Value = hfit->GetBinContent(uBin+1);
            Error = hfit->GetBinError(uBin+1);
            RanVal = rangen.Gaus(Value,Error);
            hfit->SetBinContent(uBin+1,RanVal);
            hfit->SetBinError(uBin+1,Error);
        }
        Index++;
		std::cout<<"At Index: "<<Index<<"\n";
		string FitterName = Form("FITTER_%i_Baseline_%s",Index,AddendumForPolType.Data());
		//FitterName += uniqueIdent;
		string FitterNameBase = Form("FITTER_%i_Baseline_%s",Index,AddendumForPolType.Data());
		//FitterNameBase += uniqueIdent;
		string gFitterName = Form("gFITTER_%i_Baseline_%s",Index,AddendumForPolType.Data());
  		//gFitterName += uniqueIdent;
		string gFitterNameBase = Form("gFITTER_%i_Baseline_%s",Index,AddendumForPolType.Data());
		//gFitterNameBase += uniqueIdent;
		double kFitMin = 0;
		double kFitMax = 0;
		if (DefaultVariation) {
           	      kFitMin = kMinVars[1];
                      kFitMax = kMaxVars[1];
		} else {
		      int indx1 = rangen.Integer(3);
		      int indx2 = rangen.Integer(3);
		      kFitMin = kMinVars[indx1];
		      kFitMax = kMaxVars[indx2];
		}
		TF1* FITTER = new TF1(FitterName.c_str(),MickeyFitter_pipi_boot,kFitMin,kFitMax,5);
		FITTER->SetParameter(0,sourcsize);FITTER->SetParLimits(0,0.5,3);
		FITTER->FixParameter(1,2);
		FITTER->SetParameter(2,1);
		FITTER->SetParameter(3,0);FITTER->SetParLimits(3,-1e-2,1e-2);
		if (usePol1) {
			FITTER->FixParameter(4,0);
			call_count_pol1++;
		} else {
			FITTER->SetParameter(4,0);FITTER->SetParLimits(4,-1e-3,1e-3);
			call_count_pol2++; // just in order to be sure everything works fine, else IterMax - call_count_pol1 = call_count_pol2
		}
		string ItName =  Form("%s_%i",hfit->GetName(),Index);
		//it->SetName(ItName.c_str());
		hfit->Fit(FITTER,"Q, S, R, M");
        	TF1* FITTERBASE = new TF1(FitterNameBase.c_str(),MickeyFitter_pipi_boot_Baseline,kFitMin,kFitMax,3);
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
			double fitVal = FITTER->Eval(kVal);
			double fitBaseVal = FITTERBASE->Eval(kVal);
            		gFITTER->SetPoint(uPoint,kVal,fitVal);
           		gFITTERBASE->SetPoint(uPoint,kVal,fitBaseVal);
                        //Check pol1or2 and fill arrays, also fill the mixed here.
			if (usePol1) {
				fEachBinPol1Mean[Index-1][uPoint] += fitVal;
				fEachBinPol1BaselineMean[Index-1][uPoint] += fitBaseVal;
				fEachBinPol1Stdev[Index-1][uPoint] += fitVal*fitVal;
				fEachBinPol1BaselineStdev[Index-1][uPoint] += fitBaseVal*fitBaseVal;
			} else {
				fEachBinPol2Mean[Index-1][uPoint] += fitVal;
				fEachBinPol2BaselineMean[Index-1][uPoint] += fitBaseVal;
				fEachBinPol2Stdev[Index-1][uPoint] += fitVal*fitVal;
				fEachBinPol2BaselineStdev[Index-1][uPoint] += fitBaseVal*fitBaseVal;
			}
		        uPoint++;
			fEachBinPolMixedMean[Index-1][uPoint] += fitVal;
			fEachBinPolMixedBaselineMean[Index-1][uPoint] += fitBaseVal*fitBaseVal;
			fEachBinPolMixedStdev[Index-1][uPoint] += fitVal;
			fEachBinPolMixedBaselineStdev[Index-1][uPoint] += fitBaseVal*fitBaseVal;
		}
		fgFITTERBuffer.emplace_back(gFITTER);
		fgFITTERBASEBuffer.emplace_back(gFITTERBASE);
    }
    //Should POL1 and POL2 be separated, also mark whether POL1 or POL2 were used?
    TFile* OutputFile = new TFile(TString::Format("../Output/OUTPUT_%u_%u_%u%s.root",VariationSettings[0],VariationSettings[1],VariationSettings[2],uniqueIdent.c_str()),"recreate");
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
        it->GetXaxis()->SetTitle("k* [MeV/#it{c}]");
        it->GetYaxis()->SetTitle("C(k*)");
        it->GetYaxis()->SetRangeUser(0.6, 1.8);
        it->GetXaxis()->SetRangeUser(kMin, kMax);
        it->SetTitle(Form("#splitline{%s}{<r> = %f #pm %f}",it->GetTitle(),fFITTERBuffer[Index2]->GetParameter(0),fFITTERBuffer[Index2]->GetParError(0)));
	it->Draw("pe");
	TF1* f0 = fFITTERBASEBuffer[Index2];
	f0->SetLineColor(kBlue);
	f0->Draw("same");
        auto legend2 = new TLegend(0.51,0.52,0.89,0.73);
        legend2->AddEntry(it,"CF","lep");
	legend2->AddEntry(fFITTERBuffer[Index2],"Full Fit","l");
	legend2->AddEntry(f0,"Baseline","l");
   	legend2->SetBorderSize(0);
	legend2->Draw("Same");
	string name = (it)->GetName();
        auto savenamepart="Plot_test/"+name;
	auto savename=savenamepart+".pdf";
        if (DefaultVariation) C->Print(savename.c_str());
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
           ParXValueSourceSplit.push_back(j);
           ParYValueSourceSplit.push_back(ParYValueSource[iter5]);
           ParXErrorSourceSplit.push_back(0);
           ParYErrorSourceSplit.push_back(ParYErrorSource[iter5]);
           iter5++;
        }
        auto gr = new TGraphErrors(5,&ParXValueSourceSplit[0],&ParYValueSourceSplit[0],&ParXErrorSourceSplit[0],&ParYErrorSourceSplit[0]);
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
        track++;
    }
    OutputFile->Write();
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
   //Caluculate Stdevs and fill the TGraphErrors
   for (int j = 0; j < TotalArraySize; ++j) {
    for (int i = 0; i < NumMomBins; ++i){
       //take care of meanvalues
       fEachBinPol1Mean[j][i] /= call_count_pol1;
       fEachBinPol1BaselineMean[j][i] /= call_count_pol1;
       fEachBinPol2Mean[j][i] /= call_count_pol2;
       fEachBinPol2BaselineMean[j][i] /= call_count_pol2;
       fEachBinPolMixedMean[j][i] /= IterMax;
       fEachBinPolMixedBaselineMean[j][i] /= IterMax;
       //Calculate Stdev
       fEachBinPol1Stdev[j][i] = sqrt((fEachBinPol1Stdev[j][i]-fEachBinPol1Mean[j][i]*fEachBinPol1Mean[j][i])); //Use for post processing
       fEachBinPol2Stdev[j][i] = sqrt((fEachBinPol2Stdev[j][i]-fEachBinPol2Mean[j][i]*fEachBinPol2Mean[j][i]));
       fEachBinPolMixedStdev[j][i] = sqrt((fEachBinPolMixedStdev[j][i]-fEachBinPolMixedMean[j][i]*fEachBinPolMixedMean[j][i]));
       fEachBinPol1BaselineStdev[j][i] = sqrt((fEachBinPol1BaselineStdev[j][i]-fEachBinPol1BaselineMean[j][i]*fEachBinPol1BaselineMean[j][i]));
       fEachBinPol2BaselineStdev[j][i] = sqrt((fEachBinPol2BaselineStdev[j][i]-fEachBinPol2BaselineMean[j][i]*fEachBinPol2BaselineMean[j][i]));
       fEachBinPolMixedBaselineStdev[j][i] = sqrt((fEachBinPolMixedBaselineStdev[j][i]-fEachBinPolMixedBaselineMean[j][i]*fEachBinPolMixedBaselineMean[j][i]));
    }
   }
   //Fill the TGraphErrors
   for (int j = 0; j < TotalArraySize; ++j) {
      unsigned uPoint=0;
      for (double kVal=kWidth*0.5; kVal<=kMax; kVal+=kWidth){
         //How do I deal with the different fit ranges??
         fResultPol1->SetPoint(uPoint,kVal,fEachBinPol1Mean[j][uPoint]);
         fResultPol1->SetPointError(uPoint,0.,fEachBinPol1Stdev[j][uPoint]);
         fResultPol2->SetPoint(uPoint,kVal,fEachBinPol2Mean[j][uPoint]);
         fResultPol2->SetPointError(uPoint,0.,fEachBinPol2Stdev[j][uPoint]);
         fResultPolMixed->SetPoint(uPoint,kVal,fEachBinPolMixedMean[j][uPoint]);
         fResultPolMixed->SetPointError(uPoint,0.,fEachBinPolMixedStdev[j][uPoint]);
         fResultPol1Baseline->SetPoint(uPoint,kVal,fEachBinPol1BaselineMean[j][uPoint]);
         fResultPol1Baseline->SetPointError(uPoint,0.,fEachBinPol1BaselineStdev[j][uPoint]);
         fResultPol2Baseline->SetPoint(uPoint,kVal,fEachBinPol2BaselineMean[j][uPoint]);
         fResultPol2Baseline->SetPointError(uPoint,0.,fEachBinPol2BaselineStdev[j][uPoint]);
         fResultPolMixedBaseline->SetPoint(uPoint,kVal,fEachBinPolMixedBaselineMean[j][uPoint]);
         fResultPolMixedBaseline->SetPointError(uPoint,0.,fEachBinPolMixedBaselineStdev[j][uPoint]);
	 ++uPoint;
      }
   }
   TFile* OutputFile2 = new TFile(TString::Format("../Output/RESULTS_OUTPUT_%u_%u_%u_IterIs_%i.root",VariationSettings[0],VariationSettings[1],VariationSettings[2],IterMax),"recreate");
   fResultPol1->Write("fResultPol1");
   fResultPol2->Write("fResultPol2");
   fResultPolMixed->Write("fResultPolMixed");
   fResultPol1Baseline->Write("fResultPol1Baseline");
   fResultPol2Baseline->Write("fResultPol2Baseline");
   fResultPolMixedBaseline->Write("fResultPolMixedBaseline");
   OutputFile2->Write();
   OutputFile2->Close();
   //Clean-up of all the arrays...
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


std::vector<double> CalculateLambdaParam_boot(double pur1, double pur2, bool SamePart, bool extinput, bool verbose) {

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



std::vector<TH1D*> ExtractRelevantPlotsGF_boot(string GentleFemto, int k, int m, int iVAR) {
  std::vector<int> variationIDs = {10, 11, 20, 21, 22, 23, 30, 31, 32, 33, 40, 41, 42, 43, 50, 51, 60, 61, 70, 71}; //quick fix maybe just improve the extraction procedure
	std::vector<TH1D*> fResults;
	TFile* file = TFile::Open(GentleFemto.c_str());
    TString sVAR = "";
    if(iVAR>=0) sVAR = TString::Format("Var%i",variationIDs[iVAR]);

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
            fGentleBufferAod = ExtractRelevantPlotsGF_boot(Form("../Files/AOD/CFOutput_pipi_MB_mult_%i_kT_%i.root",uMult+1,ukT+1),ukT,uMult,-1);
            TH1D* haod = fGentleBufferAod[0];
            fGentleBufferNano = ExtractRelevantPlotsGF_boot(Form("../Files/CFOutput_ForDimi_pipiVar%i_MB_mult_%i_kT_%i.root",0,uMult+1,ukT+1),ukT,uMult,0);
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



