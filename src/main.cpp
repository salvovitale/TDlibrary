/*!
 * \file main.cpp
 * \brief  this program aims to randomly test and debugging the TD Library for SU2 .
 * \author S.Vitale.
 * \version 1.0
 */

#include <iostream>
#include <stdio.h>
#include "math.h"
#include <cmath>
#include <stdlib.h>
#include "../include/fluid_model.hpp"
using namespace std;

void randInRange(int min, int max, int np, double* number);

int main(int argc, char *argv[]) {
//	double gamma = 1.4;
//	double r = 287.058;
//	double pcr = 3588550.0;
//	double tcr = 131.00;
//	double w= 0.035;
//	double P = 1013250.0;
//	double T = 588.15;

	double gamma = 1.05;
	double r = 35.149;
	double pcr = 1415000.0;
	double tcr = 564.09;
	double w= 0.529;
	double P = 500000.0;
	double T = 550.00;

	string  thlib = "RefProp";
	string  fluid = "MDM";
	int ncomp = 1;
	double* conc = new double [20];
	conc[0] = 1.0;

	CFluidModel *vw;
	vw= new CVanDerWaalsGas(gamma, r, pcr, tcr);
	CFluidModel *id;
	id= new CIdealGas(gamma, r);
	CFluidModel *pr;
	pr= new CPengRobinson(gamma, r, pcr, tcr, w);
	CFluidModel *flp;
	flp= new CFluidProp( thlib, fluid, ncomp, conc );

	vw->SetTDState_PT(P,T);
	id->SetTDState_PT(P,T);
	pr->SetTDState_PT(P,T);
	flp->SetTDState_PT(P/pow(10.0,5),T-273.15);


	double h_vw = vw->GetStaticEnergy() + vw->GetPressure()/vw->GetDensity();
	double h_id = id->GetStaticEnergy() + id->GetPressure()/id->GetDensity();
	double h_pr = pr->GetStaticEnergy() + pr->GetPressure()/pr->GetDensity();
	double h_flp = flp->GetStaticEnergy() + flp->GetPressure()/flp->GetDensity();

	cout << "h, s, rho, P, T for Van der Waals " << h_vw << " " << vw->GetEntropy()<< " " << vw->GetDensity()<< " " << vw->GetPressure()<< " " << vw->GetTemperature()<<endl;
	cout << "h, s, rho,  P, T  for Ideal Gas " << h_id << " " << id->GetEntropy()<< " " << id->GetDensity()<<" " << id->GetPressure()<<" " << id->GetTemperature()<< endl;
	cout << "h, s, rho,  P, T  for Peng Robinson " << h_pr << " " << pr->GetEntropy() << " " << pr->GetDensity()<<" " << pr->GetPressure()<<" " << pr->GetTemperature()<< endl;
	cout << "h, s, rho,  P, T  for FluidProp " << h_flp << " " << flp->GetEntropy() << " " << flp->GetDensity()<<" " << flp->GetPressure()<<" " << flp->GetTemperature()<< endl;

	vw->SetTDState_hs(3*h_vw, 1.0*vw->GetEntropy());
	id->SetTDState_hs(3*h_id, 1.0*id->GetEntropy());
	pr->SetTDState_hs(3*h_pr, 1.0*pr->GetEntropy());

	cout << "P, T for Van der Waals " << vw->GetPressure() << " " << vw->GetTemperature()<< endl;
	cout << "P, T for Ideal Gas " << id->GetPressure() << " " << id->GetTemperature()<< endl;
	cout << "P, T for for Peng Robinson " << pr->GetPressure() << " " << pr->GetTemperature()<< endl;
	cout <<"-------------------------Hello World------------------------" << endl;

	int np = 10;
	double* Pvector = new double [np];
	double* Tvector = new double [np];

	randInRange(P, 2*P, np, Pvector);
	randInRange(T, T+100, np, Tvector);

    for( int i=0; i<np; i++ ){
    	for( int j=0; j<np; j++ ){
    		vw->SetTDState_PT(Pvector[i],Tvector[j]);
    		id->SetTDState_PT(Pvector[i],Tvector[j]);
    		pr->SetTDState_PT(Pvector[i],Tvector[j]);
    		flp->SetTDState_PT(Pvector[i]/pow(10.0,5),Tvector[j]-273.15);
    		//cout << "rho " << id->GetDensity() << " " << vw->GetDensity() << " " << pr->GetDensity() << " " << flp->GetDensity() << endl;
    		cout << "c   " << id->GetSoundSpeed() << " " << vw->GetSoundSpeed() << " " << pr->GetSoundSpeed() << " " << flp->GetSoundSpeed() << endl;

    	}
    }

}


void randInRange(int min, int max, int np, double* number)
{
  for( int i=0; i<np; i++ )
  {
   number[i] = min + ( rand() / (double) RAND_MAX ) * (max - min);
   //cout << "random number " << number[i] << endl;
  }
}





