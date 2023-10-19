//2023.2.28
//This could be the final fitting code for the JJ inclusive cross section measurement
//The fitting is applied on four dimensions: J1Mass, J2Mass, ctau1, ctau2 
//The fitting is applied in the narrow mass windows [2.7,3.5]
//No vertex requirement is asked, trigger matching and JJ pair smear have been applied at the pre-selction
//Seven components will be included in the fitting in total:
//
//JJ(Prompt+Prompt): J1Mass -> CBJ1/double CB/float; J2Mass -> CBJ2/double CB/float(share the same set of parameters with the CBJ1);
//					 ctau1 -> ctau1_JJ_Gaus/double gaus/fixed by prompt MC sample;
//					 ctau2 -> ctau2_JJ_Gaus/double gaus/fixed by prompt MC sample(share the same set of parameters with the ctau1_JJ_Gaus)
//JJ(Non-Prompt+Prompt): J1Mass -> CBJ1/double CB/float; J2Mass -> CBJ2/double CB/float
//						 ctau1 -> ctau1_JJ_Convolution/convolution/fixed by non-prompt MC sample;
//						 ctau2 -> ctau2_JJ_Gaus/double gaus/fixed by prompt MC sample
//JJ(Prompt+Non-Prompt): J1Mass -> CBJ1/double CB/float; J2Mass -> CBJ2/double CB/float
//						 ctau1 -> ctau1_JJ_Gaus/double gaus/fixed by prompt MC sample;
//						 ctau2 -> ctau2_JJ_Convolution/convolution/fixed by non-prompt MC sample(share the same set of parameters with the ctau1_JJ_Convolution)
//JJ(Non-Prompt+Non-Prompt): J1Mass -> CBJ1/double CB/float; J2Mass -> CBJ2/double CB/float
//							 ctau1 -> ctau1_JJ_Convolution/convolution/fixed by non-prompt MC sample;
//							 ctau2 -> ctau2_JJ_Convolution/convolution/fixed by non-prompt MC sample
//JMuMu: J1Mass -> CBJ1/double CB/float; J2Mass -> J2_Cheb/second order chebychev/fixed by data sample;
//		 ctau1 -> ctau1_JMuMu_PDF/sum of a convolution and a gaus/fixed by side band sample;
//		 ctau2 -> ctau2_JMuMu_PDF/sum of a convolution and a gaus/fixed by side band sample(share the same set of parameters with the ctau1_JMuMu_PDF)
//MuMuJ: J1Mass -> J1_Cheb/second order chebychev/fixed by data sample(share the same set of parameters with the J2_Cheb);
//		 J2Mass -> CBJ2/double CB/float;
//		 ctau1 -> ctau1_JMuMu_PDF/sum of a convolution and a gaus/fixed by side band sample;
//		 ctau2 -> ctau2_JMuMu_PDF/sum of a convolution and a gaus/fixed by side band sample
//MuMuMuMu: J1Mass -> J1_Cheb/second order chebychev/fixed by data sample; J2Mass -> J2_Cheb/second order chebychev/fixed by data sample;
//			ctau1 -> ctau1_MuMuMuMu_Gaus/gaus/fixed by side band sample;
//			ctau2 -> ctau2_MuMuMuMu_Gaus/gaus/fixed by side band sample(share the same set of parameters with the ctau1_MuMuMuMu_Gaus)
//
//We would like to test the code by the 2016 dataset
//All the paramters can be found in the same directory (2023.2.27)
//Good luck with that!

//2023.4.4
//Huge controversy has been introduced about the time information of the combinatorial background
//It is found that the ctau1 and ctau2 shape of the JMuMu component are different
//It is considered as the shape of Jpsi and MuMu
//The fit for the combinatorial background has been redone and the strategy for the combinatorial bakground has been changed:
// 
//JMuMu: J1Mass -> CBJ1/double CB/float; J2Mass -> J2_Cheb/second order chebychev/fixed by data sample;
//		 ctau1 -> ctau1_JMuMu_PDF/sum of a convolution and a gaus/fixed by side band sample;
//		 ctau2 -> ctau2_JMuMu_Convolution/convolution/fixed by side band sample
//MuMuJ: J1Mass -> J1_Cheb/second order chebychev/fixed by data sample(share the same set of parameters with the J2_Cheb);
//		 J2Mass -> CBJ2/double CB/float;
//		 ctau1 -> ctau1_JMuMu_Convolution/convolution//fixed by side band sample (share the same set of parameters with the ctau2_JMuMu_Convolution);
//		 ctau2 -> ctau2_JMuMu_PDF/sum of a convolution and a gaus/fixed by side band sample (share the same set of parameters with the ctau1_JMuMu_PDF);
//MuMuMuMu: J1Mass -> J1_Cheb/second order chebychev/fixed by data sample; J2Mass -> J2_Cheb/second order chebychev/fixed by data sample;
//			ctau1 -> ctau1_JMuMu_Convolution/convolution/fixed by side band sample;
//			ctau2 -> ctau2_JMuMu_Convolution/convolution/fixed by side band sample;
//
//The new set of parameters can be found in the same directory (2023.4.4)
//Another set of plots can also be found with merged non-prompt components and combinatorial components
//Good luck again!

//2023.4.7
//The previous version works well, we now try to fix the height of four components:
//1. The height of the JJ(Prompt+Non-Prompt) and JJ(Non-Prompt+Prompt) should be the same;
//2. The height of the JMuMu and MuMuJ shoud be the same
//The default will be recovered to the inclusive one to make it calearer

//2023.8.2
//We try the weightted fit now
//Reference: https://root.cern.ch/doc/master/rf403__weightedevts_8C.html
//Only ctau will be fitted
//Parameters will be updated at local directory 8.2/8.5
//Information for binning fit will be introduced

//2023.8.14
// Binned fit

//2023.9.11
//For binned fit, the shape of the Jpsi mass peak will be fixed to the unbinned fit

//We try to write an automatic script

// H1. T Series
#include "TString.h"
#include "TH2.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFile.h"
#include "TROOT.h"
// H2. C Series
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
// H3. RooFit pack
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooCategory.h"
#include "RooPlot.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "RooVoigtian.h"
#include "RooCBShape.h"
#include "RooKeysPdf.h"
#include "RooAddPdf.h"
#include "RooAddition.h"
#include "RooProdPdf.h"
#include "RooProduct.h"
#include "RooFitResult.h"
#include "RooStepFunction.h"
#include "RooGenericPdf.h"
#include "TLatex.h"
#include "RooCrystalBall.h"
#include "RooFFTConvPdf.h"
#include "RooHist.h"

using namespace RooFit;

void Fit_Final(float FourMuonMass_LowCut_Value, float FourMuonMass_HighCut_Value, int MassNumber, TString Output_Name) {

	//Shape Parameters
	//For prompt J/psi on the ctau dimension
	float ctau_JJ_Gaus_Mean_Value = 0.000113;
	float ctau_JJ_Gaus_Sigma1_Value = 0.00210;
	float ctau_JJ_Gaus_Sigma2_Value = 0.00564;
	float ctau_JJ_Gaus_Frac_Value = 0.751;

	//For non-prompt J/psi on the ctau dimension
	float ctau_JJ_Convolution_Tau_Value = -26.75;
	float ctau_JJ_Convolution_Mean_Value = 0.00037;
	float ctau_JJ_Convolution_Sigma_Value = 0.00311;

	//For JMuMu on the ctau(J/psi) dimension 
	float ctau_JMuMu_J_Convolution_Tau_Value = -27.1;
	float ctau_JMuMu_J_Convolution_Mean_Value = -0.0026;
	float ctau_JMuMu_J_Convolution_Sigma_Value = 0.0035;
	float ctau_JMuMu_J_Gaus_Mean_Value = 0.00056;
	float ctau_JMuMu_J_Gaus_Sigma_Value = 0.00116;
	float ctau_JMuMu_J_Frac_Value = 0.090;

	//For JMuMu on the ctau(MuMu) dimension 
	float ctau_JMuMu_MuMu_Convolution_Tau_Value = -23.2;
	float ctau_JMuMu_MuMu_Convolution_Mean_Value = -0.0058;
	float ctau_JMuMu_MuMu_Convolution_Sigma_Value = 0.0110;

	//For Jpsi on the mass dimension
	float J_Mean_Value = 3.09309;
	float J_Sigma1_Value = 0.0234;
	float J_Sigma2_Value = 0.0585;
	float J_Alpha_Value = 2.21;
	float J_N_Value = 1.85;
	float J_Frac_Value = 0.708;

	//Mass Windows
	float LowCut = 2.95;
	float HighCut = 3.25;

	//Bin cut
	float FourMuonMass_LowCut = FourMuonMass_LowCut_Value;
	float FourMuonMass_HighCut = FourMuonMass_HighCut_Value;
	float DeltaY_LowCut = 0;
	float DeltaY_HighCut = 0;

	//Particle Massage
	const float J_Mass = 3.0969;

	//Set the variables
	RooRealVar J1Mass_("J1Mass_", "m(J/#psi_{1}) [GeV]", LowCut, HighCut);
	RooRealVar J2Mass_("J2Mass_", "m(J/#psi_{2}) [GeV]", LowCut, HighCut);
	RooRealVar x1_("x1_", "c#tau^{J/#psi1} [cm]", -0.01, 0.12);
	RooRealVar x2_("x2_", "c#tau^{J/#psi2} [cm]", -0.01, 0.12);

	RooRealVar FourMuonMass_("FourMuonMass_", "m(J/#psi_{1}J/#psi_{2}) [Gev]", 0, 55);
	RooRealVar DeltaY_("DeltaY_", "#Delta(y(J/#psi_{1}),y(J/#psi_{2}))", 0, 5);

	RooRealVar Weight_sum("Weight_sum", "Weight_sum", +0, 5000);

	//1D PDF
	//Mass - Parameters - J/Psi 
	RooRealVar J_Mean("J_Mean", "J_Mean", J_Mean_Value);  //1

	RooRealVar J_Sigma1("J_Sigma1", "J_Sigma1", J_Sigma1_Value);                //2
	RooRealVar J_Sigma2("J_Sigma2", "J_Sigma2", J_Sigma2_Value);                //3

	RooRealVar J_Alpha("J_Alpha", "J_Alpha", J_Alpha_Value);                            //4
	RooRealVar J_N("J_N", "J_N", J_N_Value);                                        //5

	RooRealVar J_Frac("J_Frac", "J_Frac", J_Frac_Value);                             //6

	//Mass1 - Function - CBJ1 - double CB
	RooAbsPdf* CBJ1_1 = new RooCBShape("CBJ1_1", "CBJ1_1", J1Mass_, J_Mean, J_Sigma1, J_Alpha, J_N);
	RooAbsPdf* CBJ1_2 = new RooCBShape("CBJ1_2", "CBJ1_2", J1Mass_, J_Mean, J_Sigma2, J_Alpha, J_N);

	RooAddPdf CBJ1("CBJ1", "CBJ1", RooArgList(*CBJ1_1, *CBJ1_2), J_Frac);

	//Mass2 - Function - CBJ2 - double CB
	RooAbsPdf* CBJ2_1 = new RooCBShape("CBJ2_1", "CBJ2_1", J2Mass_, J_Mean, J_Sigma1, J_Alpha, J_N);
	RooAbsPdf* CBJ2_2 = new RooCBShape("CBJ2_2", "CBJ2_2", J2Mass_, J_Mean, J_Sigma2, J_Alpha, J_N);

	RooAddPdf CBJ2("CBJ2", "CBJ2", RooArgList(*CBJ2_1, *CBJ2_2), J_Frac);

	//Mass - Parameters - MuMu
	//RooRealVar J_Cheb_Co1("J_Cheb_Co1", "J_Cheb_Co1", J_Cheb_Co1_Value);
	//RooRealVar J_Cheb_Co2("J_Cheb_Co2", "J_Cheb_Co2", J_Cheb_Co2_Value);

	RooRealVar J_Cheb_Co1("J_Cheb_Co1", "J_Cheb_Co1", -2.00, -3.00, 0.10);
	RooRealVar J_Cheb_Co2("J_Cheb_Co2", "J_Cheb_Co2", 0.40, -0.20, 1.00);

	//Mass1 - Function - J1_Cheb - second order cheb
	RooAbsPdf* J1_Cheb = new RooChebychev("J1_Cheb", "J1_Cheb", J1Mass_, RooArgList(J_Cheb_Co1, J_Cheb_Co2));

	//Mass2 - Function - J2_Cheb - second order cheb
	RooAbsPdf* J2_Cheb = new RooChebychev("J2_Cheb", "J2_Cheb", J2Mass_, RooArgList(J_Cheb_Co1, J_Cheb_Co2));

	//ctau - Parameters - Prompt JJ
	RooRealVar ctau_JJ_Gaus_Mean("ctau_JJ_Gaus_Mean", "ctau_JJ_Gaus_Mean", ctau_JJ_Gaus_Mean_Value);

	RooRealVar ctau_JJ_Gaus_Sigma1("ctau_JJ_Gaus_Sigma1", "ctau_JJ_Gaus_Sigma1", ctau_JJ_Gaus_Sigma1_Value);
	RooRealVar ctau_JJ_Gaus_Sigma2("ctau_JJ_Gaus_Sigma2", "ctau_JJ_Gaus_Sigma2", ctau_JJ_Gaus_Sigma2_Value);

	RooRealVar ctau_JJ_Gaus_Frac("ctau_JJ_Gaus_Frac", "ctau_JJ_Gaus_Frac", ctau_JJ_Gaus_Frac_Value);

	//ctau1 - Function - ctau1_JJ_Gaus - double gaus
	RooAbsPdf* ctau1_JJ_Gaus1 = new RooGaussian("ctau1_JJ_Gaus1", "ctau1_JJ_Gaus1", x1_, ctau_JJ_Gaus_Mean, ctau_JJ_Gaus_Sigma1);
	RooAbsPdf* ctau1_JJ_Gaus2 = new RooGaussian("ctau1_JJ_Gaus2", "ctau1_JJ_Gaus2", x1_, ctau_JJ_Gaus_Mean, ctau_JJ_Gaus_Sigma2);

	RooAddPdf ctau1_JJ_Gaus("ctau1_JJ_Gaus", "ctau1_JJ_Gaus", RooArgList(*ctau1_JJ_Gaus1, *ctau1_JJ_Gaus2), ctau_JJ_Gaus_Frac);

	//ctau2 - Function - ctau1_JJ_Gaus - double gaus
	RooAbsPdf* ctau2_JJ_Gaus1 = new RooGaussian("ctau2_JJ_Gaus1", "ctau2_JJ_Gaus1", x2_, ctau_JJ_Gaus_Mean, ctau_JJ_Gaus_Sigma1);
	RooAbsPdf* ctau2_JJ_Gaus2 = new RooGaussian("ctau2_JJ_Gaus2", "ctau2_JJ_Gaus2", x2_, ctau_JJ_Gaus_Mean, ctau_JJ_Gaus_Sigma2);

	RooAddPdf ctau2_JJ_Gaus("ctau2_JJ_Gaus", "ctau2_JJ_Gaus", RooArgList(*ctau2_JJ_Gaus1, *ctau2_JJ_Gaus2), ctau_JJ_Gaus_Frac);

	//ctau - Parameters - Non-prompt JJ
	RooRealVar ctau_JJ_Convolution_Tau("ctau_JJ_Convolution_Tau", "ctau_JJ_Convolution_Tau", ctau_JJ_Convolution_Tau_Value);
	
	RooRealVar ctau_JJ_Convolution_Mean("ctau_JJ_Convolution_Mean", "ctau_JJ_Convolution_Mean", ctau_JJ_Convolution_Mean_Value);
	RooRealVar ctau_JJ_Convolution_Sigma("ctau_JJ_Convolution_Sigma", "ctau_JJ_Convolution_Sigma", ctau_JJ_Convolution_Sigma_Value);

	RooRealVar ctau_StepValue("ctau_StepValue", "ctau_StepValue", 1);
	RooRealVar ctau_LowLimit("ctau_LowLimit", "ctau_LowLimit", 0);
	RooRealVar ctau_HighLimit("ctau_HighLimit", "ctau_HighLimit", 0.12);

	RooAbsReal* ctau1_Step = new RooStepFunction("ctau1_Step", "ctau1_Step", x1_, RooArgList(ctau_StepValue), RooArgList(ctau_LowLimit, ctau_HighLimit));
	RooAbsReal* ctau2_Step = new RooStepFunction("ctau2_Step", "ctau2_Step", x2_, RooArgList(ctau_StepValue), RooArgList(ctau_LowLimit, ctau_HighLimit));

	//ctau1 - Function - ctau1_JJ_Convolution - Convolution
	RooAbsPdf* ctau1_JJ_Exp = new RooGenericPdf("ctau1_JJ_Exp", "ctau1_JJ_Exp", "ctau1_Step * exp(x1_ * ctau_JJ_Convolution_Tau)", RooArgSet(*ctau1_Step, x1_, ctau_JJ_Convolution_Tau));
	RooAbsPdf* ctau1_JJ_Gaus_ForConvolution = new RooGaussian("ctau1_JJ_Gaus_ForConvolution", "ctau1_JJ_Gaus_ForConvolution", x1_, ctau_JJ_Convolution_Mean, ctau_JJ_Convolution_Sigma);

	RooFFTConvPdf ctau1_JJ_Convolution("ctau1_JJ_Convolution", "ctau1_JJ_Convolution", x1_, *ctau1_JJ_Exp, *ctau1_JJ_Gaus_ForConvolution);

	//ctau2 - Function - ctau2_JJ_Convolution - Convolution
	RooAbsPdf* ctau2_JJ_Exp = new RooGenericPdf("ctau2_JJ_Exp", "ctau2_JJ_Exp", "ctau2_Step * exp(x2_ * ctau_JJ_Convolution_Tau)", RooArgSet(*ctau2_Step, x2_, ctau_JJ_Convolution_Tau));
	RooAbsPdf* ctau2_JJ_Gaus_ForConvolution = new RooGaussian("ctau2_JJ_Gaus_ForConvolution", "ctau2_JJ_Gaus_ForConvolution", x2_, ctau_JJ_Convolution_Mean, ctau_JJ_Convolution_Sigma);

	RooFFTConvPdf ctau2_JJ_Convolution("ctau2_JJ_Convolution", "ctau2_JJ_Convolution", x2_, *ctau2_JJ_Exp, *ctau2_JJ_Gaus_ForConvolution);

	//ctau - Parameters - JMuMu - J/psi
	RooRealVar ctau_JMuMu_J_Convolution_Tau("ctau_JMuMu_J_Convolution_Tau", "ctau_JMuMu_J_Convolution_Tau", ctau_JMuMu_J_Convolution_Tau_Value);

	RooRealVar ctau_JMuMu_J_Convolution_Mean("ctau_JMuMu_J_Convolution_Mean", "ctau_JMuMu_J_Convolution_Mean", ctau_JMuMu_J_Convolution_Mean_Value);
	RooRealVar ctau_JMuMu_J_Convolution_Sigma("ctau_JMuMu_J_Convolution_Sigma", "ctau_JMuMu_J_Convolution_Sigma", ctau_JMuMu_J_Convolution_Sigma_Value);

	RooRealVar ctau_JMuMu_J_Gaus_Mean("ctau_JMuMu_J_Gaus_Mean", "ctau_JMuMu_J_Gaus_Mean", ctau_JMuMu_J_Gaus_Mean_Value);
	RooRealVar ctau_JMuMu_J_Gaus_Sigma("ctau_JMuMu_J_Gaus_Sigma", "ctau_JMuMu_J_Gaus_Sigma", ctau_JMuMu_J_Gaus_Sigma_Value);

	RooRealVar ctau_JMuMu_J_Frac("ctau_JMuMu_J_Frac", "ctau_JMuMu_J_Frac", ctau_JMuMu_J_Frac_Value);

	//ctau1 - Function - ctau1_JMuMu_J_PDF - Convolution+Gaus
	RooAbsPdf* ctau1_JMuMu_J_Exp = new RooGenericPdf("ctau1_JMuMu_J_Exp", "ctau1_JMuMu_J_Exp", "ctau1_Step * exp(x1_ * ctau_JMuMu_J_Convolution_Tau)", RooArgSet(*ctau1_Step, x1_, ctau_JMuMu_J_Convolution_Tau));
	RooAbsPdf* ctau1_JMuMu_J_Gaus_ForConvolution = new RooGaussian("ctau1_JMuMu_J_Gaus_ForConvolution", "ctau1_JMuMu_J_Gaus_ForConvolution", x1_, ctau_JMuMu_J_Convolution_Mean, ctau_JMuMu_J_Convolution_Sigma);

	RooFFTConvPdf ctau1_JMuMu_J_Convolution("ctau1_JMuMu_J_Convolution", "ctau1_JMuMu_J_Convolution", x1_, *ctau1_JMuMu_J_Exp, *ctau1_JMuMu_J_Gaus_ForConvolution);

	RooAbsPdf* ctau1_JMuMu_J_Gaus = new RooGaussian("ctau1_JMuMu_J_Gaus", "ctau1_JMuMu_J_Gaus2", x1_, ctau_JMuMu_J_Gaus_Mean, ctau_JMuMu_J_Gaus_Sigma);
	RooAddPdf ctau1_JMuMu_PDF("ctau1_JMuMu_PDF", "ctau1_JMuMu_PDF", RooArgList(*ctau1_JMuMu_J_Gaus, ctau1_JMuMu_J_Convolution), ctau_JMuMu_J_Frac);

	//ctau2 - Function - ctau2_JMuMu_J_PDF - Convolution+Gaus
	RooAbsPdf* ctau2_JMuMu_J_Exp = new RooGenericPdf("ctau2_JMuMu_J_Exp", "ctau2_JMuMu_J_Exp", "ctau2_Step * exp(x2_ * ctau_JMuMu_J_Convolution_Tau)", RooArgSet(*ctau2_Step, x2_, ctau_JMuMu_J_Convolution_Tau));
	RooAbsPdf* ctau2_JMuMu_J_Gaus_ForConvolution = new RooGaussian("ctau2_JMuMu_J_Gaus_ForConvolution", "ctau2_JMuMu_J_Gaus_ForConvolution", x2_, ctau_JMuMu_J_Convolution_Mean, ctau_JMuMu_J_Convolution_Sigma);

	RooFFTConvPdf ctau2_JMuMu_J_Convolution("ctau2_JMuMu_J_Convolution", "ctau2_JMuMu_J_Convolution", x2_, *ctau2_JMuMu_J_Exp, *ctau2_JMuMu_J_Gaus_ForConvolution);

	RooAbsPdf* ctau2_JMuMu_J_Gaus = new RooGaussian("ctau2_JMuMu_J_Gaus", "ctau2_JMuMu_J_Gaus2", x2_, ctau_JMuMu_J_Gaus_Mean, ctau_JMuMu_J_Gaus_Sigma);
	RooAddPdf ctau2_JMuMu_PDF("ctau2_JMuMu_PDF", "ctau2_JMuMu_PDF", RooArgList(*ctau2_JMuMu_J_Gaus, ctau2_JMuMu_J_Convolution), ctau_JMuMu_J_Frac);

	//ctau - Parameters - JMuMu - MuMu
	RooRealVar ctau_JMuMu_MuMu_Convolution_Tau("ctau_JMuMu_MuMu_Convolution_Tau", "ctau_JMuMu_MuMu_Convolution_Tau", ctau_JMuMu_MuMu_Convolution_Tau_Value);

	RooRealVar ctau_JMuMu_MuMu_Convolution_Mean("ctau_JMuMu_MuMu_Convolution_Mean", "ctau_JMuMu_MuMu_Convolution_Mean", ctau_JMuMu_MuMu_Convolution_Mean_Value);
	RooRealVar ctau_JMuMu_MuMu_Convolution_Sigma("ctau_JMuMu_MuMu_Convolution_Sigma", "ctau_JMuMu_MuMu_Convolution_Sigma", ctau_JMuMu_MuMu_Convolution_Sigma_Value);

	//ctau1 - Function - ctau1_JMuMu_Convolution - Convolution
	RooAbsPdf* ctau1_JMuMu_MuMu_Exp = new RooGenericPdf("ctau1_JMuMu_MuMu_Exp", "ctau1_JMuMu_MuMu_Exp", "ctau1_Step * exp(x1_ * ctau_JMuMu_MuMu_Convolution_Tau)", RooArgSet(*ctau1_Step, x1_, ctau_JMuMu_MuMu_Convolution_Tau));
	RooAbsPdf* ctau1_JMuMu_MuMu_Gaus_ForConvolution = new RooGaussian("ctau1_JMuMu_MuMu_Gaus_ForConvolution", "ctau1_JMuMu_MuMu_Gaus_ForConvolution", x1_, ctau_JMuMu_MuMu_Convolution_Mean, ctau_JMuMu_MuMu_Convolution_Sigma);

	RooFFTConvPdf ctau1_JMuMu_Convolution("ctau1_JMuMu_Convolution", "ctau1_JMuMu_Convolution", x1_, *ctau1_JMuMu_MuMu_Exp, *ctau1_JMuMu_MuMu_Gaus_ForConvolution);

	//ctau2 - Function - ctau2_JMuMu_Convolution - Convolution
	RooAbsPdf* ctau2_JMuMu_MuMu_Exp = new RooGenericPdf("ctau2_JMuMu_MuMu_Exp", "ctau2_JMuMu_MuMu_Exp", "ctau2_Step * exp(x2_ * ctau_JMuMu_MuMu_Convolution_Tau)", RooArgSet(*ctau2_Step, x2_, ctau_JMuMu_MuMu_Convolution_Tau));
	RooAbsPdf* ctau2_JMuMu_MuMu_Gaus_ForConvolution = new RooGaussian("ctau2_JMuMu_MuMu_Gaus_ForConvolution", "ctau2_JMuMu_MuMu_Gaus_ForConvolution", x2_, ctau_JMuMu_MuMu_Convolution_Mean, ctau_JMuMu_MuMu_Convolution_Sigma);

	RooFFTConvPdf ctau2_JMuMu_Convolution("ctau2_JMuMu_Convolution", "ctau2_JMuMu_Convolution", x2_, *ctau2_JMuMu_MuMu_Exp, *ctau2_JMuMu_MuMu_Gaus_ForConvolution);

	//4D PDF
	//J1Mass: CBJ1, *J1_Cheb
	//J2Mass: CBJ2, *J2_Cheb
	//ctau1: ctau1_JJ_Gaus, ctau1_JJ_Convolution, ctau1_JMuMu_PDF, ctau1_JMuMu_Convolution
	//ctau2: ctau2_JJ_Gaus, ctau2_JJ_Convolution, ctau2_JMuMu_PDF, ctau2_JMuMu_Convolution
	RooProdPdf PDF_JJ_PP("PDF_JJ_PP", "PDF_JJ_PP", RooArgList(CBJ1, CBJ2, ctau1_JJ_Gaus, ctau2_JJ_Gaus));

	RooProdPdf PDF_JJ_NPP("PDF_JJ_NPP", "PDF_JJ_NPP", RooArgList(CBJ1, CBJ2, ctau1_JJ_Convolution, ctau2_JJ_Gaus));
	RooProdPdf PDF_JJ_PNP("PDF_JJ_PNP", "PDF_JJ_PNP", RooArgList(CBJ1, CBJ2, ctau1_JJ_Gaus, ctau2_JJ_Convolution));
	RooProdPdf PDF_JJ_NPNP("PDF_JJ_NPNP", "PDF_JJ_NPNP", RooArgList(CBJ1, CBJ2, ctau1_JJ_Convolution, ctau2_JJ_Convolution));

	RooProdPdf PDF_J1MuMu("PDF_J1MuMu", "PDF_J1MuMu", RooArgList(CBJ1, *J2_Cheb, ctau1_JMuMu_PDF, ctau2_JMuMu_Convolution));
	RooProdPdf PDF_MuMuJ2("PDF_MuMuJ2", "PDF_MuMuJ2", RooArgList(*J1_Cheb, CBJ2, ctau1_JMuMu_Convolution, ctau2_JMuMu_PDF));
	RooProdPdf PDF_MuMuMuMu("PDF_MuMuMuMu", "PDF_MuMuMuMu", RooArgList(*J1_Cheb, *J2_Cheb, ctau1_JMuMu_Convolution, ctau2_JMuMu_Convolution));

	//Yield
	RooRealVar N_JJ_PP("N_JJ_PP", "N_JJ_PP", 2000, 0, 1E6);

	RooRealVar N_JJ_NPP("N_JJ_NPP", "N_JJ_NPP", 500, 0, 1E4);
	RooRealVar N_JJ_NPNP("N_JJ_NPNP", "N_JJ_NPNP", 1000, 0, 1E6);

	RooRealVar N_JMuMu("N_JMuMu", "N_JMuMu", 0, 0, 5000);
	RooRealVar N_MuMuMuMu("N_MuMuMuMu", "N_MuMuMuMu", 0, 0, 300);

	//Fitting function
	RooArgList PDF_List(PDF_JJ_PP, PDF_JJ_NPP, PDF_JJ_PNP, PDF_JJ_NPNP, PDF_J1MuMu, PDF_MuMuJ2, PDF_MuMuMuMu);
	RooArgList Yield_List(N_JJ_PP, N_JJ_NPP, N_JJ_NPP, N_JJ_NPNP, N_JMuMu, N_JMuMu, N_MuMuMuMu);

	RooAbsPdf* FittingFunction = new RooAddPdf("FittingFunction", "FittingFunction", PDF_List, Yield_List);

	//Mass windows cut, for insurance
	TString MassCut = "";
	MassCut = MassCut + "(J1Mass_ >" + LowCut + ") && (J1Mass_ < " + HighCut + ") && (J2Mass_ >" + LowCut + ") && (J2Mass_ < " + HighCut + ")";
	MassCut = MassCut + "&& (FourMuonMass_ > " + FourMuonMass_LowCut + ") && (FourMuonMass_ < " + FourMuonMass_HighCut + ")";
	//MassCut = MassCut + "&& (DeltaY_ > " + DeltaY_LowCut + ") && (DeltaY_ < " + DeltaY_HighCut + ")";

	//Data acquire
	TFile F_input("~/work1/jinfeng/CrossSection/CMSSW_10_2_5/src/PreSelection/Feb/PreSelection/Data/2016/Weighted_2016.root", "read");
	TTree* Tree_ = (TTree*)gROOT->FindObject("WeightedTree");

	RooArgList Parameter_List(J1Mass_, J2Mass_, x1_, x2_, FourMuonMass_, DeltaY_, Weight_sum);
	RooDataSet* Set = new RooDataSet("Set", "Set", Tree_, Parameter_List, MassCut, Weight_sum.GetName());

	//Fit
	RooFitResult* FitResult = FittingFunction->fitTo(*Set, Extended(kTRUE), Hesse(kTRUE), Strategy(1), NumCPU(4), Save(kTRUE), Minos(kFALSE));

	//Plot
	//Plot J1
	RooPlot* J1Frame = J1Mass_.frame(Title("M_{J/#psi1}"));

	Set->plotOn(J1Frame, MarkerColor(kBlack), LineColor(kBlack), Name("Dat"));
	FittingFunction->plotOn(J1Frame, LineColor(kBlue), Name("Total"));

	FittingFunction->plotOn(J1Frame, Components(PDF_JJ_PP), LineColor(kRed), Name("PP"));
	
	FittingFunction->plotOn(J1Frame, Components(PDF_JJ_NPP), LineColor(kMagenta + 4), LineStyle(kDashed), Name("NPP"));
	FittingFunction->plotOn(J1Frame, Components(PDF_JJ_PNP), LineColor(kMagenta - 1), LineStyle(kDashed), Name("PNP"));
	FittingFunction->plotOn(J1Frame, Components(PDF_JJ_NPNP), LineColor(kMagenta - 5), LineStyle(kDashed), Name("NPNP"));

	FittingFunction->plotOn(J1Frame, Components(PDF_J1MuMu), LineColor(kCyan + 3), LineStyle(kDashed), Name("J1MuMu"));
	FittingFunction->plotOn(J1Frame, Components(PDF_MuMuJ2), LineColor(kCyan - 2), LineStyle(kDashed), Name("MuMuJ2"));
	FittingFunction->plotOn(J1Frame, Components(PDF_MuMuMuMu), LineColor(kCyan - 6), LineStyle(kDashed), Name("MuMuMuMu"));
	/*
	FittingFunction->plotOn(J1Frame, Components(RooArgSet(PDF_JJ_NPP, PDF_JJ_PNP, PDF_JJ_NPNP)), LineColor(kMagenta + 4), LineStyle(kDashed), Name("NP"));
	FittingFunction->plotOn(J1Frame, Components(RooArgSet(PDF_J1MuMu, PDF_MuMuJ2, PDF_MuMuMuMu)), LineColor(kOrange), LineStyle(kDashed), Name("Combi"));
	*/
	//Pull
	RooHist* J1_Pull = J1Frame->pullHist("Dat", "Total");
	RooPlot* J1_Pull_Frame = J1Mass_.frame(Title(" "));
	J1_Pull_Frame->addPlotable(J1_Pull, "P");

	//Legend
	TLegend legend(0.7, 0.5, 0.9, 0.9);
	//TLegend legend(0.7, 0.65, 0.9, 0.9);

	legend.AddEntry(J1Frame->findObject("Dat"), "Data(2016)", "LP");
	legend.AddEntry(J1Frame->findObject("Total"), "Total PDF", "L");
	
	legend.AddEntry(J1Frame->findObject("PP"), "J/#psi_{1}(P)J/#psi_{2}(P)", "L");
	
	legend.AddEntry(J1Frame->findObject("NPP"), "J/#psi_{1}(NP)J/#psi_{2}(P)", "L");
	legend.AddEntry(J1Frame->findObject("PNP"), "J/#psi_{1}(P)J/#psi_{2}(NP)", "L");
	legend.AddEntry(J1Frame->findObject("NPNP"), "J/#psi_{1}(NP)J/#psi_{2}(NP)", "L");

	legend.AddEntry(J1Frame->findObject("J1MuMu"), "J/#psi_{1}#mu^{+}#mu^{-}", "L");
	legend.AddEntry(J1Frame->findObject("MuMuJ2"), "#mu^{+}#mu^{-}J/#psi_{2}", "L");
	legend.AddEntry(J1Frame->findObject("MuMuMuMu"), "#mu^{+}#mu^{-}#mu^{+}#mu^{-}", "L");
	/*
	legend.AddEntry(J1Frame->findObject("PP"), "J/#psi_{1}J/#psi_{2}(P)", "L");
	legend.AddEntry(J1Frame->findObject("NP"), "J/#psi_{1}J/#psi_{2}(NP)", "L");
	legend.AddEntry(J1Frame->findObject("Combi"), "Combi", "L");
	*/
	//Draw
	TString OutputNumber = "";
	OutputNumber = OutputNumber + MassNumber;

	TCanvas* J1_Canvas = new TCanvas("J1_Canvas", "J1_Canvas", 3000, 3000);
	J1_Canvas->Divide(1, 2);
	J1_Canvas->cd(1)->SetPad(0.01, 0.20, 0.95, 0.98);
	gPad->SetLeftMargin(0.15);
	J1Frame->Draw();
	legend.DrawClone();

	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.5 * J1_Canvas->GetTopMargin());
	latex.SetTextFont(62);
	latex.SetTextAlign(11);
	latex.DrawLatex(0.18, 0.85, "CMS");

	J1_Canvas->cd(2)->SetPad(0.01, 0.03, 0.95, 0.25);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.25);
	J1_Pull_Frame->GetYaxis()->SetTitle("Pull");
	J1_Pull_Frame->GetYaxis()->SetTitleSize(0.1);
	J1_Pull_Frame->GetYaxis()->SetLabelSize(0.1);
	J1_Pull_Frame->GetXaxis()->SetTitleSize(0.1);
	J1_Pull_Frame->GetXaxis()->SetLabelSize(0.1);
	J1_Pull_Frame->Draw();

	J1_Canvas->SaveAs("./Plot/J1_Result" + OutputNumber + ".pdf");
	delete J1_Canvas;

	//Plot J2
	RooPlot* J2Frame = J2Mass_.frame(Title("M_{J/#psi2}"));

	Set->plotOn(J2Frame, MarkerColor(kBlack), LineColor(kBlack), Name("Dat"));
	FittingFunction->plotOn(J2Frame, LineColor(kBlue), Name("Total"));

	FittingFunction->plotOn(J2Frame, Components(PDF_JJ_PP), LineColor(kRed), Name("PP"));
	
	FittingFunction->plotOn(J2Frame, Components(PDF_JJ_NPP), LineColor(kMagenta + 4), LineStyle(kDashed), Name("NPP"));
	FittingFunction->plotOn(J2Frame, Components(PDF_JJ_PNP), LineColor(kMagenta - 1), LineStyle(kDashed), Name("PNP"));
	FittingFunction->plotOn(J2Frame, Components(PDF_JJ_NPNP), LineColor(kMagenta - 5), LineStyle(kDashed), Name("NPNP"));

	FittingFunction->plotOn(J2Frame, Components(PDF_J1MuMu), LineColor(kCyan + 3), LineStyle(kDashed), Name("J1MuMu"));
	FittingFunction->plotOn(J2Frame, Components(PDF_MuMuJ2), LineColor(kCyan - 2), LineStyle(kDashed), Name("MuMuJ2"));
	FittingFunction->plotOn(J2Frame, Components(PDF_MuMuMuMu), LineColor(kCyan - 6), LineStyle(kDashed), Name("MuMuMuMu"));
	/*
	FittingFunction->plotOn(J2Frame, Components(RooArgSet(PDF_JJ_NPP, PDF_JJ_PNP, PDF_JJ_NPNP)), LineColor(kMagenta + 4), LineStyle(kDashed), Name("NP"));
	FittingFunction->plotOn(J2Frame, Components(RooArgSet(PDF_J1MuMu, PDF_MuMuJ2, PDF_MuMuMuMu)), LineColor(kOrange), LineStyle(kDashed), Name("Combi"));
	*/
	//Pull
	RooHist* J2_Pull = J2Frame->pullHist("Dat", "Total");
	RooPlot* J2_Pull_Frame = J2Mass_.frame(Title(" "));
	J2_Pull_Frame->addPlotable(J2_Pull, "P");

	//Draw
	TCanvas* J2_Canvas = new TCanvas("J2_Canvas", "J2_Canvas", 3000, 3000);
	J2_Canvas->Divide(1, 2);
	J2_Canvas->cd(1)->SetPad(0.01, 0.20, 0.95, 0.98);
	gPad->SetLeftMargin(0.15);
	J2Frame->Draw();
	legend.DrawClone();

	latex.SetTextSize(0.5 * J2_Canvas->GetTopMargin());
	latex.SetTextFont(62);
	latex.SetTextAlign(11);
	latex.DrawLatex(0.18, 0.85, "CMS");

	J2_Canvas->cd(2)->SetPad(0.01, 0.03, 0.95, 0.25);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.25);
	J2_Pull_Frame->GetYaxis()->SetTitle("Pull");
	J2_Pull_Frame->GetYaxis()->SetTitleSize(0.1);
	J2_Pull_Frame->GetYaxis()->SetLabelSize(0.1);
	J2_Pull_Frame->GetXaxis()->SetTitleSize(0.1);
	J2_Pull_Frame->GetXaxis()->SetLabelSize(0.1);
	J2_Pull_Frame->Draw();

	J2_Canvas->SaveAs("./Plot/J2_Result" + OutputNumber + ".pdf");
	delete J2_Canvas;

	//Plot ctau1
	RooPlot* ctau1Frame = x1_.frame(Title("c#tau^{J/#psi1}"), Range(-0.01, 0.10));

	Set->plotOn(ctau1Frame, MarkerColor(kBlack), LineColor(kBlack), Name("Dat"));
	FittingFunction->plotOn(ctau1Frame, LineColor(kBlue), Name("Total"));

	FittingFunction->plotOn(ctau1Frame, Components(PDF_JJ_PP), LineColor(kRed), Name("PP"));
	
	FittingFunction->plotOn(ctau1Frame, Components(PDF_JJ_NPP), LineColor(kMagenta + 4), LineStyle(kDashed), Name("NPP"));
	FittingFunction->plotOn(ctau1Frame, Components(PDF_JJ_PNP), LineColor(kMagenta - 1), LineStyle(kDashed), Name("PNP"));
	FittingFunction->plotOn(ctau1Frame, Components(PDF_JJ_NPNP), LineColor(kMagenta - 5), LineStyle(kDashed), Name("NPNP"));

	FittingFunction->plotOn(ctau1Frame, Components(PDF_J1MuMu), LineColor(kCyan + 3), LineStyle(kDashed), Name("J1MuMu"));
	FittingFunction->plotOn(ctau1Frame, Components(PDF_MuMuJ2), LineColor(kCyan - 2), LineStyle(kDashed), Name("MuMuJ2"));
	FittingFunction->plotOn(ctau1Frame, Components(PDF_MuMuMuMu), LineColor(kCyan - 6), LineStyle(kDashed), Name("MuMuMuMu"));
	/*
	FittingFunction->plotOn(ctau1Frame, Components(RooArgSet(PDF_JJ_NPP, PDF_JJ_PNP, PDF_JJ_NPNP)), LineColor(kMagenta + 4), LineStyle(kDashed), Name("NP"));
	FittingFunction->plotOn(ctau1Frame, Components(RooArgSet(PDF_J1MuMu, PDF_MuMuJ2, PDF_MuMuMuMu)), LineColor(kOrange), LineStyle(kDashed), Name("Combi"));
	*/
	//Pull
	RooHist* ctau1_Pull = ctau1Frame->pullHist("Dat", "Total");
	RooPlot* ctau1_Pull_Frame = x1_.frame(Title(" "), Range(-0.01, 0.10));
	ctau1_Pull_Frame->addPlotable(ctau1_Pull, "P");

	//Draw
	TCanvas* ctau1_Canvas = new TCanvas("ctau1_Canvas", "ctau1_Canvas", 3000, 3000);
	ctau1_Canvas->Divide(1, 2);
	ctau1_Canvas->cd(1)->SetPad(0.01, 0.20, 0.95, 0.98);
	gPad->SetLeftMargin(0.15);
	gPad->SetLogy();
	ctau1Frame->Draw();
	legend.DrawClone();

	latex.SetTextSize(0.5 * ctau1_Canvas->GetTopMargin());
	latex.SetTextFont(62);
	latex.SetTextAlign(11);
	latex.DrawLatex(0.60, 0.85, "CMS");

	ctau1_Canvas->cd(2)->SetPad(0.01, 0.03, 0.95, 0.25);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.25);
	ctau1_Pull_Frame->GetYaxis()->SetTitle("Pull");
	ctau1_Pull_Frame->GetYaxis()->SetTitleSize(0.1);
	ctau1_Pull_Frame->GetYaxis()->SetLabelSize(0.1);
	ctau1_Pull_Frame->GetXaxis()->SetTitleSize(0.1);
	ctau1_Pull_Frame->GetXaxis()->SetLabelSize(0.1);
	ctau1_Pull_Frame->Draw();

	ctau1_Canvas->SaveAs("./Plot/ctau1_Result" + OutputNumber + ".pdf");
	delete ctau1_Canvas;

	//Plot ctau2
	RooPlot* ctau2Frame = x2_.frame(Title("c#tau^{J/#psi2}"), Range(-0.01, 0.10));

	Set->plotOn(ctau2Frame, MarkerColor(kBlack), LineColor(kBlack), Name("Dat"));
	FittingFunction->plotOn(ctau2Frame, LineColor(kBlue), Name("Total"));

	FittingFunction->plotOn(ctau2Frame, Components(PDF_JJ_PP), LineColor(kRed), Name("PP"));
	
	FittingFunction->plotOn(ctau2Frame, Components(PDF_JJ_NPP), LineColor(kMagenta + 4), LineStyle(kDashed), Name("NPP"));
	FittingFunction->plotOn(ctau2Frame, Components(PDF_JJ_PNP), LineColor(kMagenta - 1), LineStyle(kDashed), Name("PNP"));
	FittingFunction->plotOn(ctau2Frame, Components(PDF_JJ_NPNP), LineColor(kMagenta - 5), LineStyle(kDashed), Name("NPNP"));

	FittingFunction->plotOn(ctau2Frame, Components(PDF_J1MuMu), LineColor(kCyan + 3), LineStyle(kDashed), Name("J1MuMu"));
	FittingFunction->plotOn(ctau2Frame, Components(PDF_MuMuJ2), LineColor(kCyan - 2), LineStyle(kDashed), Name("MuMuJ2"));
	FittingFunction->plotOn(ctau2Frame, Components(PDF_MuMuMuMu), LineColor(kCyan - 6), LineStyle(kDashed), Name("MuMuMuMu"));
	/*
	FittingFunction->plotOn(ctau2Frame, Components(RooArgSet(PDF_JJ_NPP, PDF_JJ_PNP, PDF_JJ_NPNP)), LineColor(kMagenta + 4), LineStyle(kDashed), Name("NP"));
	FittingFunction->plotOn(ctau2Frame, Components(RooArgSet(PDF_J1MuMu, PDF_MuMuJ2, PDF_MuMuMuMu)), LineColor(kOrange), LineStyle(kDashed), Name("Combi"));
	*/
	//Pull
	RooHist* ctau2_Pull = ctau2Frame->pullHist("Dat", "Total");
	RooPlot* ctau2_Pull_Frame = x2_.frame(Title(" "), Range(-0.01, 0.10));
	ctau2_Pull_Frame->addPlotable(ctau2_Pull, "P");

	//Draw
	TCanvas* ctau2_Canvas = new TCanvas("ctau2_Canvas", "ctau2_Canvas", 3000, 3000);
	ctau2_Canvas->Divide(1, 2);
	ctau2_Canvas->cd(1)->SetPad(0.01, 0.20, 0.95, 0.98);
	gPad->SetLeftMargin(0.15);
	gPad->SetLogy();
	ctau2Frame->Draw();
	legend.DrawClone();

	latex.SetTextSize(0.5 * ctau2_Canvas->GetTopMargin());
	latex.SetTextFont(62);
	latex.SetTextAlign(11);
	latex.DrawLatex(0.60, 0.85, "CMS");

	ctau2_Canvas->cd(2)->SetPad(0.01, 0.03, 0.95, 0.25);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.25);
	ctau2_Pull_Frame->GetYaxis()->SetTitle("Pull");
	ctau2_Pull_Frame->GetYaxis()->SetTitleSize(0.1);
	ctau2_Pull_Frame->GetYaxis()->SetLabelSize(0.1);
	ctau2_Pull_Frame->GetXaxis()->SetTitleSize(0.1);
	ctau2_Pull_Frame->GetXaxis()->SetLabelSize(0.1);
	ctau2_Pull_Frame->Draw();

	ctau2_Canvas->SaveAs("./Plot/ctau2_Result" + OutputNumber + ".pdf");
	delete ctau2_Canvas;

	/*
	//Print
	cout << "//////////////////" << endl;
	cout << "Status: " << FitResult->status() << endl;
	cout << "JJ(Prompt+Prompt): " << N_JJ_PP.getValV() << " +/- " << N_JJ_PP.getError() << endl;
	
	cout << "JJ(NonPrompt+Prompt): " << N_JJ_NPP.getValV() << " +/- " << N_JJ_NPP.getError() << endl;
	//cout << "JJ(Prompt+NonPrompt): " << N_JJ_PNP.getValV() << " +/- " << N_JJ_PNP.getError() << endl;
	cout << "JJ(NonPrompt+NonPrompt): " << N_JJ_NPNP.getValV() << " +/- " << N_JJ_NPNP.getError() << endl;

	cout << "JMuMu: " << N_JMuMu.getValV() << " +/- " << N_JMuMu.getError() << endl;
	//cout << "MuMuJ: " << N_MuMuJ2.getValV() << " +/- " << N_MuMuJ2.getError() << endl;
	cout << "MuMuMuMu: " << N_MuMuMuMu.getValV() << " +/- " << N_MuMuMuMu.getError() << endl;
	*/

	//Print the result to the target file
	FILE* OutputFile = NULL;
	OutputFile = fopen(Output_Name + ".txt", "w");

	fprintf(OutputFile, "%i: %i\n", MassNumber, FitResult->status());
	fprintf(OutputFile, "PP: %f +/- %f\n", N_JJ_PP.getValV(), N_JJ_PP.getError());
	fprintf(OutputFile, "PNP: %f +/- %f\n", N_JJ_NPP.getValV(), N_JJ_NPP.getError());
	fprintf(OutputFile, "NPNP: %f +/- %f\n", N_JJ_NPNP.getValV(), N_JJ_NPNP.getError());
	fprintf(OutputFile, "JMuMu: %f +/- %f\n", N_JMuMu.getValV(), N_JMuMu.getError());
	fprintf(OutputFile, "MuMuMuMu: %f +/- %f\n", N_MuMuMuMu.getValV(), N_MuMuMuMu.getError());

	fclose(OutputFile);
}
