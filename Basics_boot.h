#include <vector>
#include "TF1.h"
#include "TH1F.h"
#include "TFile.h"
#include "DLM_CkDecomposition.h"
#include "TRandom3.h"

const double Mass_Pich = 139.57;
const double Mass_Pi0b = 134.98;

class TGraph;

void SetTheLambdaParamsToDML_CkDecompObj_boot (DLM_CkDecomposition* CkDe_PiPi, bool verbose=false, int iVAR=0, TRandom3 SEED=0);
std::vector<TH1D*> ExtractRelevantPlotsGF_boot(string GentleFemto, int k, int m, int iVAR=-1);
std::vector<double> CalculateLambdaParam_boot(double pur1=0.98, double pur2=0, bool SamePart=true, bool extinput=false, bool verbose=false);

void Basics_PiPiCATS_boot(const bool& Identical, const bool& WithCoulomb, unsigned* VariationSettings);
void Compare_AOD_nano();
