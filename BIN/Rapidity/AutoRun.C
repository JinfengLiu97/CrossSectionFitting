#include "Fit_Final.C"
#include "TString.h"

void AutoRun() {

	TString Output_Name = "Record";

	float YBin[7] = { 0.0,0.5,1.0,1.5,2.0,2.5,3.0 };

	for (int YNumber = 1; YNumber < 7; YNumber++) {

			Fit_Final(YBin[YNumber - 1], YBin[YNumber], YNumber, Output_Name);
	}

}