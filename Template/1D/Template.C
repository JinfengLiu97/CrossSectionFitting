//2023.8.23
//Template fit to acquire the DPS/SPS fraction

//2023.9.5
//Rapidity will be added

//2023.9.20
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

	//Scale factor (16)
	double ScaleFactor = 0.001 / (37.9 * 0.05961 * 0.05961);

	//Data histogram
	float DataMassBinContent[] = { 2577,635,852,1208,1158,609,477 };
	float DataMassBinError[] = { 60,30,30,40,40,30,30 };
	auto DataMassHistogram = new TH1F("DataMassHistogram", "", 7, 7.5, 42.5);
	DataMassHistogram->SetStats(kFALSE);

	float DataYBinContent[] = { 4079,1585,1039,768,480 };
	float DataYBinError[] = { 70,50,40,30,30};
	auto DataYHistogram = new TH1F("DataYHistogram", "", 5, 0.0, 2.5);
	DataYHistogram->SetStats(kFALSE);

	for (int i = 1; i < 8; i++) {
		DataMassHistogram->SetBinContent(i, DataMassBinContent[i - 1]);
		DataMassHistogram->SetBinError(i, DataMassBinError[i - 1]);

		DataYHistogram->SetBinContent(i, DataYBinContent[i - 1]);
		DataYHistogram->SetBinError(i, DataYBinError[i - 1]);

	}

	DataMassHistogram->Scale(ScaleFactor);
	DataYHistogram->Scale(ScaleFactor);

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

	TFile* SPS_File = TFile::Open("~/work1/jinfeng/CrossSection/CMSSW_10_2_5/src/PreSelection/Feb/PreSelection/MC/NoVertex/Weighted_SPS.root");
	TTree* SPS_Tree = new TTree("SPS_Tree", "SPS_Tree");
	SPS_File->GetObject("WeightedTree", SPS_Tree);
	SPS_Tree->SetBranchAddress("J1Mass_", &SPS_J1Mass);
	SPS_Tree->SetBranchAddress("J2Mass_", &SPS_J2Mass);
	SPS_Tree->SetBranchAddress("FourMuonMass_", &SPS_FourMuonMass);
	SPS_Tree->SetBranchAddress("DeltaY_", &SPS_DeltaY);
	SPS_Tree->SetBranchAddress("Weight_sum", &SPS_Weight);

	TFile* DPS_File = TFile::Open("~/work1/jinfeng/CrossSection/CMSSW_10_2_5/src/PreSelection/Feb/PreSelection/MC/NoVertex/Weighted_DPS.root");
	TTree* DPS_Tree = new TTree("DPS_Tree", "DPS_Tree");
	DPS_File->GetObject("WeightedTree", DPS_Tree);
	DPS_Tree->SetBranchAddress("J1Mass_", &DPS_J1Mass);
	DPS_Tree->SetBranchAddress("J2Mass_", &DPS_J2Mass);
	DPS_Tree->SetBranchAddress("FourMuonMass_", &DPS_FourMuonMass);
	DPS_Tree->SetBranchAddress("DeltaY_", &DPS_DeltaY);
	DPS_Tree->SetBranchAddress("Weight_sum", &DPS_Weight);

	int SPS_Length = SPS_Tree->GetEntries();
	int DPS_Length = DPS_Tree->GetEntries();

	auto SPSMassHistogram = new TH1F("SPSMassHistogram", "", 7, 7.5, 42.5);
	auto DPSMassHistogram = new TH1F("DPSMassHistogram", "", 7, 7.5, 42.5);

	auto SPSYHistogram = new TH1F("SPSYHistogram", "", 5, 0.0, 2.5);
	auto DPSYHistogram = new TH1F("DPSYHistogram", "", 5, 0.0, 2.5);

	SPSMassHistogram->SetStats(kFALSE);
	DPSMassHistogram->SetStats(kFALSE);
	SPSYHistogram->SetStats(kFALSE);
	DPSYHistogram->SetStats(kFALSE);

	for (int i = 0; i < SPS_Length; i++) {
		SPS_Tree->GetEntry(i);

		if ((SPS_J1Mass > 2.95) && (SPS_J1Mass < 3.25) && (SPS_J2Mass > 2.95) && (SPS_J2Mass < 3.25) && (SPS_Weight < 500)) {

			SPSMassHistogram->Fill(SPS_FourMuonMass, SPS_Weight);
			SPSYHistogram->Fill(SPS_DeltaY, SPS_Weight);
		}
	}

	SPSMassHistogram->Scale(ScaleFactor);
	SPSYHistogram->Scale(ScaleFactor);

	for (int i = 0; i < DPS_Length; i++) {
		DPS_Tree->GetEntry(i);

		if ((DPS_J1Mass > 2.95) && (DPS_J1Mass < 3.25) && (DPS_J2Mass > 2.95) && (DPS_J2Mass < 3.25) && (DPS_Weight<500)) {

			DPSMassHistogram->Fill(DPS_FourMuonMass, DPS_Weight);
			DPSYHistogram->Fill(DPS_DeltaY, DPS_Weight);
		}
	}

	DPSMassHistogram->Scale(ScaleFactor);
	DPSYHistogram->Scale(ScaleFactor);

	//Prepare the template
	RooRealVar FourMuonMass("FourMuonMass", "M(J/#psi_{1}J/#psi_{2}) [Gev]", 0, 45);

	RooDataHist SPS_MassTemplate("SPS_MassTemplate", "SPS_MassTemplate", FourMuonMass, SPSMassHistogram);
	RooDataHist DPS_MassTemplate("DPS_MassTemplate", "DPS_MassTemplate", FourMuonMass, DPSMassHistogram);
	RooDataHist Data_MassTemplate("Data_MassTemplate", "Data_MassTemplate", FourMuonMass, DataMassHistogram);

	RooHistPdf SPS_MassPDF("SPS_MassPDF", "SPS_MassPDF", RooArgSet(FourMuonMass), SPS_MassTemplate);
	RooHistPdf DPS_MassPDF("DPS_MassPDF", "DPS_MassPDF", RooArgSet(FourMuonMass), DPS_MassTemplate);

	RooRealVar DeltaY("DeltaY", "#Delta(y(J/#psi_{1}),y(J/#psi_{2}))", 0, 3);

	RooDataHist SPS_YTemplate("SPS_YTemplate", "SPS_YTemplate", DeltaY, SPSYHistogram);
	RooDataHist DPS_YTemplate("DPS_YTemplate", "DPS_YTemplate", DeltaY, DPSYHistogram);
	RooDataHist Data_YTemplate("Data_YTemplate", "Data_YTemplate", DeltaY, DataYHistogram);

	RooHistPdf SPS_YPDF("SPS_YPDF", "SPS_YPDF", RooArgSet(DeltaY), SPS_YTemplate);
	RooHistPdf DPS_YPDF("DPS_YPDF", "DPS_YPDF", RooArgSet(DeltaY), DPS_YTemplate);

	//Fit
	RooRealVar DPS_MassFraction("DPS_MassFraction", "DPS_MassFraction", 0.5, 0, 1);

	RooAddPdf MassPDF("MassPDF", "MassPDF", RooArgList(DPS_MassPDF, SPS_MassPDF), RooArgList(DPS_MassFraction));

	RooRealVar DPS_YFraction("DPS_YFraction", "DPS_YFraction", 0.5, 0, 1);

	RooAddPdf YPDF("YPDF", "YPDF", RooArgList(DPS_YPDF, SPS_YPDF), RooArgList(DPS_YFraction));

	RooFitResult* MassFitResult = MassPDF.fitTo(Data_MassTemplate, Extended(kTRUE), Hesse(kTRUE), Strategy(1), NumCPU(4), Save(kTRUE), Minos(kFALSE));
	RooFitResult* YFitResult = YPDF.fitTo(Data_YTemplate, Extended(kTRUE), Hesse(kTRUE), Strategy(1), NumCPU(4), Save(kTRUE), Minos(kFALSE));

	//Plot
	//M
	RooPlot* MassFrame = FourMuonMass.frame(Title("M(J/#psi_{1}J/#psi_{2}) [Gev]"));

	Data_MassTemplate.plotOn(MassFrame, LineColor(kBlack), Name("Data"));

	MassPDF.plotOn(MassFrame, LineColor(kRed), LineWidth(4), Name("Total"));
	MassPDF.plotOn(MassFrame, Components(SPS_MassPDF), LineColor(kMagenta + 2), Name("SPS"));
	MassPDF.plotOn(MassFrame, Components(DPS_MassPDF), LineColor(kBlue), Name("DPS"));

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

	Data_YTemplate.plotOn(YFrame, LineColor(kBlack), Name("Data"));

	YPDF.plotOn(YFrame, LineColor(kRed), LineWidth(4), Name("Total"));
	YPDF.plotOn(YFrame, Components(SPS_YPDF), LineColor(kMagenta + 2), Name("SPS"));
	YPDF.plotOn(YFrame, Components(DPS_YPDF), LineColor(kBlue), Name("DPS"));

	//Pull
	RooHist* YPull = YFrame->pullHist("Data", "Total");
	RooPlot* YPullFrame = DeltaY.frame(Title(" "));
	YPullFrame->addPlotable(YPull, "P");

	//Draw
	//Mass
	TCanvas* MassCanvas = new TCanvas("MassCanvas", "MassCanvas", 3000, 3000);
	MassCanvas->Divide(1, 2);
	MassCanvas->cd(1)->SetPad(0.01, 0.20, 0.95, 0.98);
	gPad->SetLeftMargin(0.15);
	MassFrame->GetYaxis()->SetTitle("#sigma per bin [pb]");
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
	YFrame->GetYaxis()->SetTitle("#sigma per bin [pb]");
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

	//Print
	cout << "//////////////" << endl;
	cout << "Mass: " << endl;
	cout << "Status: " << MassFitResult->status() << endl;
	cout << "DPS fraction: " << DPS_MassFraction.getValV() << " +/- " << DPS_MassFraction.getError() << endl;

	cout << "Rapidity: " << endl;
	cout << "Status: " << YFitResult->status() << endl;
	cout << "DPS fraction: " << DPS_YFraction.getValV() << " +/- " << DPS_YFraction.getError() << endl;

}

