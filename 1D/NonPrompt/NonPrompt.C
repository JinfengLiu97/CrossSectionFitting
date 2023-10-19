
//10.5
//To acquire the shape parameters for the (possible) forth diemnsion and fifth dimension
//The forth dimension is the significance of the vertex distance
//The fifth dimension is the significance of the LxyPV
//The convolution of a gaussian and an exponent will be used for these two dimensions for all the components

//10.6
//Add another two dimensions, lxyPV and ctau, we will determine the shapes for all these dimensions all at once
//For prompt sample, the double gaussian will be applied
//For b decay sample. the double-side crystal ball will be applied
//It has to worked in 12_4_0 then

//10.24
//This version is for the prompt sample without four muon vertex cut
//Convolution of gaus and exponent will be applied on all the dimensions

//2023.2.16
// Plotting updated and trigger match added

//2023.8.2
//We try the weightted fit now
//Reference: https://root.cern.ch/doc/master/rf403__weightedevts_8C.html
//Only ctau will be fitted
//ctau range has been changed to [-0.01, 0.12/0.10]

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
#include "TPad.h"
#include "TVirtualPad.h"

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

void NonPrompt() {

	//Mass windows
	float LowCut = 2.95;
	float HighCut = 3.25;

	//Bin cut
	float FourMuonMass_LowCut = 7.5;
	float FourMuonMass_HighCut = 12.5;
	float DeltaY_LowCut = 0;
	float DeltaY_HighCut = 1;

	RooRealVar J1Mass_("J1Mass_", "m(J/#psi_{1}) [GeV]", 2.95, 3.25);
	RooRealVar J2Mass_("J2Mass_", "m(J/#psi_{2}) [GeV]", 2.95, 3.25);

	RooRealVar x1_("x1_", "c#tau^{J/#psi1} [cm]", -0.01, 0.12);

	RooRealVar TriggerMatch_("TriggerMatch_", "TriggerMatch_", -1, 2);
	RooRealVar Weight_sum("Weight_sum", "Weight_sum", +0, 5000);

	RooRealVar FourMuonMass_("FourMuonMass_", "m(J/#psi_{1}J/#psi_{2}) [Gev]", 0, 55);
	RooRealVar DeltaY_("DeltaY_", "#Delta(y(J/#psi_{1}),y(J/#psi_{2}))", 0, 5);

	//ctau
	//Non-prompt (convolution)
	RooRealVar ctau_Tau("ctau_Tau", "ctau_Tau", -20, -40, -0);

	RooRealVar ctau_StepValue("ctau_StepValue", "ctau_StepValue", 1);
	RooRealVar ctau_LowerLimit("ctau_LowerLimit", "ctau_LowerLimit", 0);
	RooRealVar ctau_HigherLimit("ctau_HigherLimit", "ctau_HigherLimit", 0.12);

	RooRealVar ctau_Mean("ctau_Mean", "ctau_Mean", 0.01, 0, 0.02);
	RooRealVar ctau_Sigma("ctau_Sigma", "ctau_Sigma", 0.01, 0.001, 0.02);

	RooAbsReal* ctau_Step = new RooStepFunction("ctau_Step", "ctau_Step", x1_, RooArgList(ctau_StepValue), RooArgList(ctau_LowerLimit, ctau_HigherLimit));
	RooAbsPdf* ctau_Exp = new RooGenericPdf("ctau_Exp", "ctau_Exp", "ctau_Step * exp(x1_ * ctau_Tau)", RooArgSet(*ctau_Step, x1_, ctau_Tau));
	RooAbsPdf* ctau_Gaus = new RooGaussian("ctau_Gaus", "ctau_Gaus", x1_, ctau_Mean, ctau_Sigma);

	RooFFTConvPdf ctau_Convolution("ctau_Convolution", "ctau_Convolution", x1_, *ctau_Exp, *ctau_Gaus);

	//For insurance, we add a mass windows cut here:
	TString MassCut = "";
	MassCut = MassCut + "(J1Mass_ >" + LowCut + ") && (J1Mass_ < " + HighCut + ") && (J2Mass_ >" + LowCut + ") && (J2Mass_ < " + HighCut + ")";
	MassCut = MassCut + "&& (TriggerMatch_ == 1)";
	//MassCut = MassCut + "&& (FourMuonMass_ > " + FourMuonMass_LowCut + ") && (FourMuonMass_ < " + FourMuonMass_HighCut + ")";
	//MassCut = MassCut + "&& (DeltaY_ > " + DeltaY_LowCut + ") && (DeltaY_ < " + DeltaY_HighCut + ")";

	//Data acquire
	TFile F_input("~/work1/jinfeng/CrossSection/CMSSW_10_2_5/src/PreSelection/Feb/PreSelection/MC/NoVertex/Weighted_B.root", "read");
	TTree* Tree_ = (TTree*)gROOT->FindObject("WeightedTree");
	RooDataSet* Set = new RooDataSet("Set", "Set", Tree_, RooArgList(J1Mass_, J2Mass_, x1_, TriggerMatch_, FourMuonMass_, DeltaY_, Weight_sum), MassCut, Weight_sum.GetName());

	//Fit
	RooFitResult* ctau_FitResult = ctau_Convolution.fitTo(*Set, Hesse(kTRUE), Strategy(1), NumCPU(4), Save(kTRUE), Minos(kFALSE));
	
	//Plot
	RooPlot* ctau_Frame = x1_.frame(Title("c#tau^{J/#psi1}"), Range(-0.01, 0.10));

	//Line
	Set->plotOn(ctau_Frame, MarkerColor(kBlack), LineColor(kBlack), Name("Dat"));

	ctau_Convolution.plotOn(ctau_Frame, LineColor(kBlue), Name("Total"));
	
	//Pull
	RooHist* ctau_Pull = ctau_Frame->pullHist();
	RooPlot* ctau_Pull_Frame = x1_.frame(Range(-0.01, 0.10), Title(" "));
	ctau_Pull_Frame->addPlotable(ctau_Pull, "P");

	//Legend
	TLegend legend(0.7, 0.8, 0.9, 0.9);
	legend.AddEntry(ctau_Frame->findObject("Dat"), "MC", "LP");
	legend.AddEntry(ctau_Frame->findObject("Total"), "Fitting", "L");

	//Draw
	//ctau
	TCanvas* ctau_Canvas = new TCanvas("ctau_Canvas", "ctau_Canvas", 2000, 2000);
	ctau_Canvas->Divide(1, 2);
	ctau_Canvas->cd(1)->SetPad(0.01, 0.20, 0.95, 0.98);
	gPad->SetLeftMargin(0.15);
	gPad->SetLogy();
	ctau_Frame->Draw();
	legend.DrawClone();

	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.5 * ctau_Canvas->GetTopMargin());
	latex.SetTextFont(62);
	latex.SetTextAlign(11);
	latex.DrawLatex(0.25, 0.70, "CMS");
	/*
	latex.SetTextSize(0.4 * ctau_Canvas->GetTopMargin());
	latex.SetTextFont(52);
	latex.DrawLatex(0.30, 0.85, "Preliminary");
	*/

	ctau_Canvas->cd(2)->SetPad(0.01, 0.03, 0.95, 0.25);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.25);
	ctau_Pull_Frame->GetYaxis()->SetTitle("Pull");
	ctau_Pull_Frame->GetYaxis()->SetTitleSize(0.1);
	ctau_Pull_Frame->GetYaxis()->SetLabelSize(0.1);
	ctau_Pull_Frame->GetXaxis()->SetTitleSize(0.1);
	ctau_Pull_Frame->GetXaxis()->SetLabelSize(0.1);
	ctau_Pull_Frame->Draw();
	ctau_Canvas->SaveAs("NonPrompt_ctau_Result.pdf");


	//Print the result
	cout << "ctau: " << ctau_FitResult->status() << endl;
	cout << "NLL: " << ctau_FitResult->minNll() << endl;
	cout << "Chi2: " << ctau_Frame->chiSquare() << endl;
	cout << "Tau: " << ctau_Tau.getValV() << " +/- " << ctau_Tau.getError() << endl;
	cout << "Mean: " << ctau_Mean.getValV() << " +/- " << ctau_Mean.getError() << endl;
	cout << "Sigma: " << ctau_Sigma.getValV() << " +/- " << ctau_Sigma.getError() << endl;
	
}
