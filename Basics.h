#include <vector>
#include "TF1.h"
#include "TH1F.h"
#include "TFile.h"
#include "DLM_CkDecomposition.h"

const double Mass_Pi0 = 134.98;
const double Mass_Pic = 139.57;
const double Mass_p = 938.272;
const double Mass_L = 1115.683;
const double Mass_Xi= 1314.86;

class TGraph;

double Basics_Potential_Usmani(double* pars);
double Basics_Source_Gauss(double* Pars);
//double CATS_FIT_PL(double* x, double* pars);
void SetTheLambdaParamsToDML_DecompObj (DLM_CkDecomposition* CkDe_PiPi, bool verbose=false);
std::vector<TH1D*> ExtractRelevantPlotsGF(string GentleFemto, int k, int m);
std::vector<double> CalculateLambdaParam(double pur1=0.99, double pur2=0, bool SamePart=true, bool extinput=false, bool verbose=true);
std::vector<double> CalculateLambdaParam2(double pur1=0.99, double pur2=0, bool SamePart=true, bool extinput=false, bool verbose=true);

TGraph* Basics_PiPiTheory(const bool& Identical, const bool& WithCoulomb);
TGraph* Basics_PiPiCATS(const bool& Identical, const bool& WithCoulomb);
void Basics_PiPiCATS2(const bool& Identical, const bool& WithCoulomb);
void Basics_ProtonLambda();
