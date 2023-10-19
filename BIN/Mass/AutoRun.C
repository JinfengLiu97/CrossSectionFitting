#include "Fit_Final.C"
#include "TString.h"

void AutoRun() {

	TString Output_Name = "Record";

	float MassBin[8] = { 7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5 };

	for (int MassNumber = 1; MassNumber < 8; MassNumber++) {
		Fit_Final(MassBin[MassNumber - 1], MassBin[MassNumber], MassNumber, Output_Name);
	}

}