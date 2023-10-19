//2023.9.5
//2D template fit code to extract the DPS/SPS fraction

//2023.10.18
//Using the DPS fraction to construct the function instead of the height of SPS/DPS

#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TFrame.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TTree.h"
#include "TMath.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFFTConvPdf.h"
using namespace RooFit;

void Template() {

	//Data histogram
	int DataBinContent[5][6] = {
		{2372,318,394,428,273,143},
		{261,210,217,330,241,107},
		{0,78,150,163,291,87},
		{0,0,3,97,208,106},
		{0,0,0,0,18,172}
	};

	int DataBinError[5][6] = {
		{50,17,20,30,19,14},
		{17,17,17,20,18,13},
		{0,11,14,15,19,12},
		{0,0,0,10,17,14},
		{0,0,0,0,5,16}
	};

	auto DataHistogram = new TH2F("DataHistogram", "", 5, 0, 2.5, 6, 7.5, 37.5);
	DataHistogram->SetLineColor(kBlack);
	DataHistogram->SetStats(kFALSE);

	for (int MassNumber = 1; MassNumber < 7; MassNumber++) {
		for (int YNumber = 1; YNumber < 6; YNumber++) {

			DataHistogram->SetBinContent(YNumber, MassNumber, DataBinContent[YNumber - 1][MassNumber - 1]);
			DataHistogram->SetBinError(YNumber, MassNumber, DataBinError[YNumber - 1][MassNumber - 1]);
		}
	}

	//SPS&DPS histograms
	float SPS_J1Mass;
	float SPS_J2Mass;
	float SPS_FourMuonMass;
	float SPS_DeltaY;
	float SPS_Weight;

	float DPS_J1Mass;
	float DPS_J2Mass;
	float DPS_FourMuonMass;
	float DPS_DeltaY;
	float DPS_Weight;

	TFile* SPS_File = TFile::Open("~/work1/jinfeng/CrossSection/CMSSW_10_2_5/src/PreSelection/Feb/PreSelection//MC/NoVertex/Weighted_SPS.root");
	TTree* SPS_Tree = new TTree("SPS_Tree", "SPS_Tree");
	SPS_File->GetObject("WeightedTree", SPS_Tree);
	SPS_Tree->SetBranchAddress("J1Mass_", &SPS_J1Mass);
	SPS_Tree->SetBranchAddress("J2Mass_", &SPS_J2Mass);
	SPS_Tree->SetBranchAddress("FourMuonMass_", &SPS_FourMuonMass);
	SPS_Tree->SetBranchAddress("DeltaY_", &SPS_DeltaY);
	SPS_Tree->SetBranchAddress("Weight_sum", &SPS_Weight);

	TFile* DPS_File = TFile::Open("~/work1/jinfeng/CrossSection/CMSSW_10_2_5/src/PreSelection/Feb/PreSelection//MC/NoVertex/Weighted_DPS.root");
	TTree* DPS_Tree = new TTree("DPS_Tree", "DPS_Tree");
	DPS_File->GetObject("WeightedTree", DPS_Tree);
	DPS_Tree->SetBranchAddress("J1Mass_", &DPS_J1Mass);
	DPS_Tree->SetBranchAddress("J2Mass_", &DPS_J2Mass);
	DPS_Tree->SetBranchAddress("FourMuonMass_", &DPS_FourMuonMass);
	DPS_Tree->SetBranchAddress("DeltaY_", &DPS_DeltaY);
	DPS_Tree->SetBranchAddress("Weight_sum", &DPS_Weight);

	int SPS_Length = SPS_Tree->GetEntries();
	int DPS_Length = DPS_Tree->GetEntries();

	auto SPSHistogram = new TH2F("SPSHistogram", "", 5, 0, 2.5, 6, 7.5, 37.5);
	SPSHistogram->SetStats(kFALSE);
	auto DPSHistogram = new TH2F("DPSHistogram", "", 5, 0, 2.5, 6, 7.5, 37.5);
	DPSHistogram->SetStats(kFALSE);

	for (int i = 0; i < SPS_Length; i++) {
		SPS_Tree->GetEntry(i);

		if ((SPS_J1Mass > 2.95) && (SPS_J1Mass < 3.25) && (SPS_J2Mass > 2.95) && (SPS_J2Mass < 3.25) && (SPS_Weight < 500)) {

			SPSHistogram->Fill(SPS_DeltaY, SPS_FourMuonMass, SPS_Weight);
		}
	}

	for (int i = 0; i < DPS_Length; i++) {
		DPS_Tree->GetEntry(i);

		if ((DPS_J1Mass > 2.95) && (DPS_J1Mass < 3.25) && (DPS_J2Mass > 2.95) && (DPS_J2Mass < 3.25) && (DPS_Weight < 500)) {

			DPSHistogram->Fill(DPS_DeltaY, DPS_FourMuonMass, DPS_Weight);
		}
	}

	//Prepare the template
	RooRealVar FourMuonMass("FourMuonMass", "M(J/#psi_{1}J/#psi_{2}) [Gev]", 0, 37.5);
	RooRealVar DeltaY("DeltaY", "#Delta(y(J/#psi_{1}),y(J/#psi_{2}))", 0, 2.5);

	RooDataHist SPS_Template("SPS_Template", "SPS_Template", RooArgList(DeltaY, FourMuonMass), SPSHistogram);
	RooDataHist DPS_Template("DPS_Template", "DPS_Template", RooArgList(DeltaY, FourMuonMass), DPSHistogram);
	RooDataHist Data_Template("Data_Template", "Data_Template", RooArgList(DeltaY, FourMuonMass), DataHistogram);

	RooHistPdf SPS_PDF("SPS_PDF", "SPS_PDF", RooArgSet(DeltaY, FourMuonMass), SPS_Template);
	RooHistPdf DPS_PDF("DPS_PDF", "DPS_PDF", RooArgSet(DeltaY, FourMuonMass), DPS_Template);

	//Fit
	RooRealVar DPS_Fraction("DPS_Fraction", "DPS_Fraction", 0.5, 0, 1);

	RooAddPdf PDF("PDF", "PDF", RooArgList(DPS_PDF, SPS_PDF), RooArgList(DPS_Fraction));

	RooFitResult* FitResult = PDF.fitTo(Data_Template, Extended(kTRUE), Hesse(kTRUE), Strategy(1), NumCPU(4), Save(kTRUE), Minos(kFALSE));

	//Plot
	//Mass
	RooPlot* MassFrame = FourMuonMass.frame(Title("M(J/#psi_{1}J/#psi_{2}) [Gev]"));

	Data_Template.plotOn(MassFrame, LineColor(kBlack), Name("Data"));

	PDF.plotOn(MassFrame, LineColor(kRed), LineWidth(4), Name("Total"));
	PDF.plotOn(MassFrame, Components(SPS_PDF), LineColor(kMagenta + 2), Name("SPS"));
	PDF.plotOn(MassFrame, Components(DPS_PDF), LineColor(kBlue), Name("DPS"));

	//Pull
	RooHist* MassPull = MassFrame->pullHist("Data", "Total");
	RooPlot* MassPullFrame = FourMuonMass.frame(Title(" "));
	MassPullFrame->addPlotable(MassPull, "P");

	//Legend
	TLegend legend(0.7, 0.65, 0.9, 0.9);
	legend.AddEntry(MassFrame->findObject("Data"), "Data(2016)", "L");

	legend.AddEntry(MassFrame->findObject("Total"), "Total PDF", "L");
	legend.AddEntry(MassFrame->findObject("SPS"), "SPS", "L");
	legend.AddEntry(MassFrame->findObject("DPS"), "DPS", "L");

	//Y
	RooPlot* YFrame = DeltaY.frame(Title("#Delta(y(J/#psi_{1}),y(J/#psi_{2}))"));

	Data_Template.plotOn(YFrame, LineColor(kBlack), Name("Data"));

	PDF.plotOn(YFrame, LineColor(kRed), LineWidth(4), Name("Total"));
	PDF.plotOn(YFrame, Components(SPS_PDF), LineColor(kMagenta + 2), Name("SPS"));
	PDF.plotOn(YFrame, Components(DPS_PDF), LineColor(kBlue), Name("DPS"));

	//Pull
	RooHist* YPull = YFrame->pullHist("Data", "Total");
	RooPlot* YPullFrame = DeltaY.frame(Title(" "));
	YPullFrame->addPlotable(YPull, "P");

	//2D
	//auto Frame2D = new RooPlot(FourMuonMass, DeltaY);

	TH1* PDF_Histogram = PDF.createHistogram("PDF_Histogram", DeltaY, YVar(FourMuonMass));
	
	TH1* SPS_PDF_Histogram = SPS_PDF.createHistogram("SPS_PDF_Histogram", DeltaY, YVar(FourMuonMass));
	TH1* DPS_PDF_Histogram = DPS_PDF.createHistogram("DPS_PDF_Histogram", DeltaY, YVar(FourMuonMass));

	PDF_Histogram->SetLineColor(kRed);
	PDF_Histogram->SetLineWidth(4);
	PDF_Histogram->SetStats(kFALSE);

	SPS_PDF_Histogram->Scale(1/SPS_Height.getValV());
	SPS_PDF_Histogram->SetLineColor(kMagenta + 2);
	SPS_PDF_Histogram->SetFillColor(kMagenta + 2);
	SPS_PDF_Histogram->SetStats(kFALSE);
	DPS_PDF_Histogram->Scale(1/DPS_Height.getValV());
	DPS_PDF_Histogram->SetLineColor(kBlue);
	DPS_PDF_Histogram->SetFillColor(kBlue);
	DPS_PDF_Histogram->SetStats(kFALSE);

	//Draw
	//Mass
	TCanvas* MassCanvas = new TCanvas("MassCanvas", "MassCanvas", 3000, 3000);
	MassCanvas->Divide(1, 2);
	MassCanvas->cd(1)->SetPad(0.01, 0.20, 0.95, 0.98);
	gPad->SetLeftMargin(0.15);
	MassFrame->Draw();
	legend.DrawClone();

	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.5 * MassCanvas->GetTopMargin());
	latex.SetTextFont(62);
	latex.SetTextAlign(11);
	latex.DrawLatex(0.55, 0.85, "CMS");

	MassCanvas->cd(2)->SetPad(0.01, 0.03, 0.95, 0.25);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.25);
	MassPullFrame->GetYaxis()->SetTitle("Pull");
	MassPullFrame->GetYaxis()->SetTitleSize(0.1);
	MassPullFrame->GetYaxis()->SetLabelSize(0.1);
	MassPullFrame->GetXaxis()->SetTitleSize(0.1);
	MassPullFrame->GetXaxis()->SetLabelSize(0.1);
	MassPullFrame->Draw();

	MassCanvas->SaveAs("Mass_Result.pdf");

	//Y
	TCanvas* YCanvas = new TCanvas("YCanvas", "YCanvas", 3000, 3000);
	YCanvas->Divide(1, 2);
	YCanvas->cd(1)->SetPad(0.01, 0.20, 0.95, 0.98);
	gPad->SetLeftMargin(0.15);
	YFrame->Draw();
	legend.DrawClone();

	latex.SetTextSize(0.5 * YCanvas->GetTopMargin());
	latex.SetTextFont(62);
	latex.SetTextAlign(11);
	latex.DrawLatex(0.55, 0.85, "CMS");

	YCanvas->cd(2)->SetPad(0.01, 0.03, 0.95, 0.25);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.25);
	YPullFrame->GetYaxis()->SetTitle("Pull");
	YPullFrame->GetYaxis()->SetTitleSize(0.1);
	YPullFrame->GetYaxis()->SetLabelSize(0.1);
	YPullFrame->GetXaxis()->SetTitleSize(0.1);
	YPullFrame->GetXaxis()->SetLabelSize(0.1);
	YPullFrame->Draw();

	YCanvas->SaveAs("Y_Result.pdf");

	//2D
	TCanvas* Canvas2D = new TCanvas("Canvas2D", "Canvas2D", 3000, 3000);
	gPad->SetLeftMargin(0.15);
	DataHistogram->Draw("lego");
	PDF_Histogram->Draw("same, lego");
	//SPS_PDF_Histogram->Draw("same, lego");
	//DPS_PDF_Histogram->Draw("same, lego");
	legend.DrawClone();

	latex.SetTextSize(0.5 * Canvas2D->GetTopMargin());
	latex.SetTextFont(62);
	latex.SetTextAlign(11);
	latex.DrawLatex(0.55, 0.85, "CMS");

	Canvas2D->SaveAs("2D_Result.png");

	//Print
	cout << "////////////////" << endl;
	cout << "Status: " << FitResult->status() << endl;
	cout << "DPS fraction: " << DPS_Fraction.getValV() << " +/- " << DPS_Fraction.getError() << endl;
}