#include "Fit_Final.C"
#include "TString.h"

void AutoRun() {

	TString Output_Name = "Record";

	float YBin[7] = { 0.0,0.5,1.0,1.5,2.0,2.5,3.0 };
	float MassBin[8] = { 7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5 };

	for (int YNumber = 1; YNumber < 7; YNumber++) {
		for (int MassNumber = 1; MassNumber < 8; MassNumber++) {

			Fit_Final(MassBin[MassNumber - 1], MassBin[MassNumber], YBin[YNumber - 1], YBin[YNumber], YNumber, MassNumber, Output_Name);
		}
	}

}