/*!
 * \file main.cpp
 * \brief  this program aims to randomly test and debugging the TD Library for SU2 .
 * \author S.Vitale.
 * \version 1.0
 */

#include <iostream>
#include <stdio.h>
#include "../include/fluid_model.hpp"
using namespace std;

int main(int argc, char *argv[]) {
	double gamma = 1.4;
	double r = 287.058;
	double pcr = 3588550.0;
	double tcr = 131.00;
	double P = 101325.0;
	double T = 288.15;

	CFluidModel *vw;
	vw= new CVanDerWaalsGas(gamma, r, pcr, tcr);
	CFluidModel *id;
	id= new CIdealGas(gamma, r);

	vw->SetTDState_PT(P,T);
	id->SetTDState_PT(P,T);


	double h_vw = vw->GetStaticEnergy() + vw->GetPressure()/vw->GetDensity();
	double h_id = id->GetStaticEnergy() + id->GetPressure()/id->GetDensity();

	cout << "h, s, rho for Van der Waals " << h_vw << " " << vw->GetEntropy()<< " " << vw->GetDensity() <<endl;
	cout << "h, s, rho  for Ideal Gas " << h_id << " " << id->GetEntropy()<< " " << id->GetDensity()<< endl;

	vw->SetTDState_hs(h_vw, vw->GetEntropy());
	id->SetTDState_hs(h_id,id->GetEntropy() );

	cout << "P, T for Van der Waals " << vw->GetPressure() << " " << vw->GetTemperature()<< endl;
	cout << "P, T for Ideal Gas " << id->GetPressure() << " " << id->GetTemperature()<< endl;
	cout <<"-------------------------Hello World------------------------" << endl;
	getchar();

}

