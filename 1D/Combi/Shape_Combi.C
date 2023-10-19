//2023.2.27
//This is the code to obtain the shape of the combinatorial background
//The background is acquired by the side band cut [2.7-2.95], [3.25-3.5]
//First check the shape then decide what function to be used
//As for now, the time dimension is the ctau
//Since we have smeared the two Jpsis in the pre-selction part, there are only two kind of combinatorial background
//
//  | |\| |   |\| |\|
//  |\| |\|   | | | |
//  | |\| |   |\| |\|
//   JMuMu    MuMuMuMu
// 
//And, for the same reason, we only need the shape of one dimension

//The shape for the JMuMu component on the lifetime dimension should be a convolution

//2023.4.4
//It is obvious that we misunderstand the principle of the time shape of the combinatorial background
//For J1MuMu, the ctau2 dimension keeps the MuMu time information although the smearing has been carried out
//Similar for the MuMuJ2 on the ctau1 dimension
//The strategy now is:
//  \ J2  
//  | |1|2|3|   
//  | |4|5|6|   
//  | |7|8|9|   
//    ----->J1
//2 and 8 to represent the J1MuMu component and fit the two dimension
//The ctau1 would be kept as a merging of a convolution and a gaussian (shape5)
//The ctau2 would be a convolution (shape6)
//In the final fitting, the arrangement of the combinatorial background would be:
//				ctau1	ctau2
//  J1MuMu		 S5		 S6
//  MuMuJ2		 S6		 S5
//MuMuMuMuMu	 S6		 S6

//2023.8.2
//We try the weightted fit now
//Reference: https://root.cern.ch/doc/master/rf403__weightedevts_8C.html
//Only ctau will be fitted
//ctau range has been changed to [-0.01, 0.12/0.10]
//Fittings on mass dimensions will be removed

//2023.10.18
//We will try a four muon mass/Delta rapidity windows cut here
//To check the shift of the shape

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
#include "RooCrystalBall.h"
#include "RooFFTConvPdf.h"
#include "RooHist.h"

using namespace RooFit;

void Shape_Combi() {

	//Mass windows
	float LowCut = 2.95;
	float HighCut = 3.25;

	//Bin cut
	float FourMuonMass_LowCut = 7.5;
	float FourMuonMass_HighCut = 12.5;
	float DeltaY_LowCut = 0;
	float DeltaY_HighCut = 1;

	//Particle Massage
	const float J_Mass = 3.0969;

	//Set the variable 
	RooRealVar J1Mass_("J1Mass_", "m(J/#psi_{1}) [GeV]", 2.7, 3.5);
	RooRealVar J2Mass_("J2Mass_", "m(J/#psi_{2}) [GeV]", 2.7, 3.5);
	RooRealVar x1_("x1_", "c#tau^{J/#psi1} [cm]", -0.01, 0.12);
	RooRealVar x2_("x2_", "c#tau^{J/#psi2} [cm]", -0.01, 0.12);
	RooRealVar Weight_sum("Weight_sum", "Weight_sum", +0, 5000);

	RooRealVar FourMuonMass_("FourMuonMass_", "m(J/#psi_{1}J/#psi_{2}) [Gev]", 0, 55);
	RooRealVar DeltaY_("DeltaY_", "#Delta(y(J/#psi_{1}),y(J/#psi_{2}))", 0, 5);

	//PDF - J1MuMu - ctau1
	RooRealVar ctau1_JMuMu_Tau("ctau1_JMuMu_Tau", "ctau1_JMuMu_Tau", -20, -40, -0);

	RooRealVar ctau_JMuMu_StepValue("ctau_JMuMu_StepValue", "ctau_JMuMu_StepValue", 1);
	RooRealVar ctau_JMuMu_LowerLimit("ctau_JMuMu_LowerLimit", "ctau_JMuMu_LowerLimit", 0);
	RooRealVar ctau_JMuMu_HigherLimit("ctau_JMuMu_HigherLimit", "ctau_JMuMu_HigherLimit", 0.12);

	RooRealVar ctau1_JMuMu_Mean1("ctau1_JMuMu_Mean1", "ctau1_JMuMu_Mean1", -0.01, -0.015, 0.005);
	RooRealVar ctau1_JMuMu_Sigma1("ctau1_JMuMu_Sigma1", "ctau1_JMuMu_Sigma1", 0.01, 0.001, 0.02);

	RooAbsReal* ctau1_JMuMu_Step = new RooStepFunction("ctau1_JMuMu_Step", "ctau1_JMuMu_Step", x1_, RooArgList(ctau_JMuMu_StepValue), RooArgList(ctau_JMuMu_LowerLimit, ctau_JMuMu_HigherLimit));
	RooAbsPdf* ctau1_JMuMu_Exp = new RooGenericPdf("ctau1_JMuMu_Exp", "ctau1_JMuMu_Exp", "ctau1_JMuMu_Step * exp(x1_ * ctau1_JMuMu_Tau)", RooArgSet(*ctau1_JMuMu_Step, x1_, ctau1_JMuMu_Tau));
	RooAbsPdf* ctau1_JMuMu_Gaus1 = new RooGaussian("ctau1_JMuMu_Gaus1", "ctau1_JMuMu_Gaus1", x1_, ctau1_JMuMu_Mean1, ctau1_JMuMu_Sigma1);

	RooRealVar ctau1_JMuMu_Mean2("ctau1_JMuMu_Mean2", "ctau1_JMuMu_Mean2", -0.01, -0.015, 0.02);
	RooRealVar ctau1_JMuMu_Sigma2("ctau1_JMuMu_Sigma2", "ctau1_JMuMu_Sigma2", 0.001, 0.0005, 0.005);
	RooAbsPdf* ctau1_JMuMu_Gaus2 = new RooGaussian("ctau1_JMuMu_Gaus2", "ctau1_JMuMu_Gaus2", x1_, ctau1_JMuMu_Mean2, ctau1_JMuMu_Sigma2);

	RooFFTConvPdf ctau1_JMuMu_Convolution("ctau1_JMuMu_Convolution", "ctau1_JMuMu_Convolution", x1_, *ctau1_JMuMu_Exp, *ctau1_JMuMu_Gaus1);

	RooRealVar ctau1_JMuMu_Frac("ctau1_JMuMu_Frac", "ctau1_JMuMu_Frac", 0.4, 0, 1);
	RooAddPdf ctau1_JMuMu_PDF("ctau1_JMuMu_PDF", "ctau1_JMuMu_PDF", RooArgList(*ctau1_JMuMu_Gaus2, ctau1_JMuMu_Convolution), ctau1_JMuMu_Frac);

	//PDF - J1MuMu - ctau2
	RooRealVar ctau2_JMuMu_Tau("ctau2_JMuMu_Tau", "ctau2_JMuMu_Tau", -20, -40, -0);

	RooAbsReal* ctau2_JMuMu_Step = new RooStepFunction("ctau2_JMuMu_Step", "ctau2_JMuMu_Step", x2_, RooArgList(ctau_JMuMu_StepValue), RooArgList(ctau_JMuMu_LowerLimit, ctau_JMuMu_HigherLimit));
	
	RooRealVar ctau2_JMuMu_Mean("ctau2_JMuMu_Mean", "ctau2_JMuMu_Mean", -0.01, -0.015, 0.02);
	RooRealVar ctau2_JMuMu_Sigma("ctau2_JMuMu_Sigma", "ctau2_JMuMu_Sigma", 0.01, 0.001, 0.02);
	
	RooAbsPdf* ctau2_JMuMu_Exp = new RooGenericPdf("ctau2_JMuMu_Exp", "ctau2_JMuMu_Exp", "ctau2_JMuMu_Step * exp(x2_ * ctau2_JMuMu_Tau)", RooArgSet(*ctau2_JMuMu_Step, x2_, ctau2_JMuMu_Tau));
	RooAbsPdf* ctau2_JMuMu_Gaus = new RooGaussian("ctau2_JMuMu_Gaus", "ctau2_JMuMu_Gaus", x2_, ctau2_JMuMu_Mean, ctau2_JMuMu_Sigma);

	RooFFTConvPdf ctau2_JMuMu_Convolution("ctau2_JMuMu_Convolution", "ctau2_JMuMu_Convolution", x2_, *ctau2_JMuMu_Exp, *ctau2_JMuMu_Gaus);

	//Input file
	TFile F_input("~/work1/jinfeng/CrossSection/CMSSW_10_2_5/src/PreSelection/Feb/PreSelection/Data/2016/Weighted_2016.root", "read");
	TTree* Tree_ = (TTree*)gROOT->FindObject("WeightedTree");

	RooArgList ParameterList(J1Mass_, J2Mass_, x1_, x2_, FourMuonMass_, DeltaY_, Weight_sum);

	//The side band mass windows 
	TString Cut_JMuMu = "";
	Cut_JMuMu = Cut_JMuMu + "((J2Mass_ <" + LowCut + ") || (J2Mass_ >" + HighCut + ")) && (J1Mass_ >" + LowCut + ") && (J1Mass_ <" + HighCut + ")";
	//Cut_JMuMu = Cut_JMuMu + "&& (FourMuonMass_ > " + FourMuonMass_LowCut + ") && (FourMuonMass_ < " + FourMuonMass_HighCut + ")";
	//Cut_JMuMu = Cut_JMuMu + "&& (DeltaY_ > " + DeltaY_LowCut + ") && (DeltaY_ < " + DeltaY_HighCut + ")";

	RooDataSet* JMuMu_Set = new RooDataSet("JMuMu_Set", "JMuMu_Set", Tree_, ParameterList, Cut_JMuMu, Weight_sum.GetName());

	//Fit
	RooFitResult* ctau1_JMuMu_FitResult = ctau1_JMuMu_PDF.fitTo(*JMuMu_Set, Hesse(kTRUE), Strategy(1), NumCPU(4), Save(kTRUE), Minos(kFALSE));
	RooFitResult* ctau2_JMuMu_FitResult = ctau2_JMuMu_Convolution.fitTo(*JMuMu_Set, Hesse(kTRUE), Strategy(1), NumCPU(4), Save(kTRUE), Minos(kFALSE));

	//Plot
	//JMuMu - ctau1
	RooPlot* x1Frame_JMuMu = x1_.frame(Title("c#tau^{J/#psi1}"), Range(-0.01, 0.10));
	JMuMu_Set->plotOn(x1Frame_JMuMu, MarkerColor(kBlack), Name("Dat"));

	ctau1_JMuMu_PDF.plotOn(x1Frame_JMuMu, LineColor(kBlue), Name("Total"));

	//Pull
	RooHist* ctau1_JMuMu_Pull = x1Frame_JMuMu->pullHist();
	RooPlot* ctau1_JMuMu_Pull_Frame = x1_.frame(Title(" "), Range(-0.01, 0.10));
	ctau1_JMuMu_Pull_Frame->addPlotable(ctau1_JMuMu_Pull, "P");

	TLegend legend(0.7, 0.8, 0.9, 0.9);
	legend.AddEntry(x1Frame_JMuMu->findObject("Dat"), "Data", "LP");
	legend.AddEntry(x1Frame_JMuMu->findObject("Total"), "Fitting", "L");

	//Draw
	TCanvas* x1Canvas_JMuMu = new TCanvas("x1Canvas_JMuMu", "x1Canvas_JMuMu", 3000, 2000);
	x1Canvas_JMuMu->Divide(1, 2);
	x1Canvas_JMuMu->cd(1)->SetPad(0.01, 0.20, 0.95, 0.98);
	gPad->SetLeftMargin(0.15);
	gPad->SetLogy();
	x1Frame_JMuMu->Draw();
	legend.DrawClone();

	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.5 * x1Canvas_JMuMu->GetTopMargin());
	latex.SetTextFont(62);
	latex.SetTextAlign(11);
	latex.DrawLatex(0.60, 0.85, "CMS");

	x1Canvas_JMuMu->cd(2)->SetPad(0.01, 0.03, 0.95, 0.25);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.25);
	ctau1_JMuMu_Pull_Frame->GetYaxis()->SetTitle("Pull");
	ctau1_JMuMu_Pull_Frame->GetYaxis()->SetTitleSize(0.1);
	ctau1_JMuMu_Pull_Frame->GetYaxis()->SetLabelSize(0.1);
	ctau1_JMuMu_Pull_Frame->GetXaxis()->SetTitleSize(0.1);
	ctau1_JMuMu_Pull_Frame->GetXaxis()->SetLabelSize(0.1);
	ctau1_JMuMu_Pull_Frame->Draw();

	x1Canvas_JMuMu->SaveAs("JMuMu_x1.pdf");

	//JMuMu - ctau2
	RooPlot* x2Frame_JMuMu = x2_.frame(Title("c#tau^{J/#psi2}"), Range(-0.01, 0.10));
	JMuMu_Set->plotOn(x2Frame_JMuMu, MarkerColor(kBlack), Name("Dat"));

	ctau2_JMuMu_Convolution.plotOn(x2Frame_JMuMu, LineColor(kBlue), Name("Total"));

	//Pull
	RooHist* ctau2_JMuMu_Pull = x2Frame_JMuMu->pullHist();
	RooPlot* ctau2_JMuMu_Pull_Frame = x2_.frame(Title(" "), Range(-0.01, 0.10));
	ctau2_JMuMu_Pull_Frame->addPlotable(ctau2_JMuMu_Pull, "P");

	//Draw
	TCanvas* x2Canvas_JMuMu = new TCanvas("x2Canvas_JMuMu", "x2Canvas_JMuMu", 3000, 2000);
	x2Canvas_JMuMu->Divide(1, 2);
	x2Canvas_JMuMu->cd(1)->SetPad(0.01, 0.20, 0.95, 0.98);
	gPad->SetLeftMargin(0.15);
	gPad->SetLogy();
	x2Frame_JMuMu->Draw();
	legend.DrawClone();

	latex.SetTextSize(0.5 * x2Canvas_JMuMu->GetTopMargin());
	latex.SetTextFont(62);
	latex.SetTextAlign(11);
	latex.DrawLatex(0.60, 0.85, "CMS");

	x2Canvas_JMuMu->cd(2)->SetPad(0.01, 0.03, 0.95, 0.25);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.25);
	ctau2_JMuMu_Pull_Frame->GetYaxis()->SetTitle("Pull");
	ctau2_JMuMu_Pull_Frame->GetYaxis()->SetTitleSize(0.1);
	ctau2_JMuMu_Pull_Frame->GetYaxis()->SetLabelSize(0.1);
	ctau2_JMuMu_Pull_Frame->GetXaxis()->SetTitleSize(0.1);
	ctau2_JMuMu_Pull_Frame->GetXaxis()->SetLabelSize(0.1);
	ctau2_JMuMu_Pull_Frame->Draw();

	x2Canvas_JMuMu->SaveAs("JMuMu_x2.pdf");


	//Print
	cout << "//////////////////" << endl;
	cout << "JMuMu - ctau1: " << ctau1_JMuMu_FitResult->status() << endl;
	cout << "NLL: " << ctau1_JMuMu_FitResult->minNll() << endl;
	cout << "Chi2: " << x1Frame_JMuMu->chiSquare() << endl;
	cout << "Fraction(Gaussia/Convolution): " << ctau1_JMuMu_Frac.getValV() << " +/- " << ctau1_JMuMu_Frac.getError() << endl;
	cout << "Gaussian: " << endl;
	cout << "Mean: " << ctau1_JMuMu_Mean2.getValV() << " +/- " << ctau1_JMuMu_Mean2.getError() << endl;
	cout << "Sigma: " << ctau1_JMuMu_Sigma2.getValV() << " +/- " << ctau1_JMuMu_Sigma2.getError() << endl;
	cout << "Convolution: " << endl;
	cout << "Tau: " << ctau1_JMuMu_Tau.getValV() << " +/- " << ctau1_JMuMu_Tau.getError() << endl;
	cout << "Mean: " << ctau1_JMuMu_Mean1.getValV() << " +/- " << ctau1_JMuMu_Mean1.getError() << endl;
	cout << "Sigma: " << ctau1_JMuMu_Sigma1.getValV() << " +/- " << ctau1_JMuMu_Sigma1.getError() << endl;

	cout << "//////////////////" << endl;
	cout << "JMuMu - ctau2: " << ctau2_JMuMu_FitResult->status() << endl;
	cout << "NLL: " << ctau2_JMuMu_FitResult->minNll() << endl;
	cout << "Chi2: " << x2Frame_JMuMu->chiSquare() << endl;
	cout << "Convolution: " << endl;
	cout << "Tau: " << ctau2_JMuMu_Tau.getValV() << " +/- " << ctau2_JMuMu_Tau.getError() << endl;
	cout << "Mean: " << ctau2_JMuMu_Mean.getValV() << " +/- " << ctau2_JMuMu_Mean.getError() << endl;
	cout << "Sigma: " << ctau2_JMuMu_Sigma.getValV() << " +/- " << ctau2_JMuMu_Sigma.getError() << endl;

}
