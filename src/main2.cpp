/*!
 * \file main2.cpp
 * \brief  this program is to compute total/static conditions for SU2 config file.
 * \author M. Pini, S. Vitale
 * \version 1.0
 */

#include <iostream>
#include <stdio.h>
#include "math.h"
#include <cmath>
#include <stdlib.h>
#include "../include/fluid_model.hpp"
using namespace std;

double Energy_hs ( double h, CFluidModel *FluidModel, double ht, double s, double Ma );
void SetTotalTDState_prho ( CFluidModel *FluidModel, double P, double rho, double Mach );

typedef double (*fct)(...);  // Pointer to a Function returning a Float.
int zbrac( fct func, double *x1, double *x2, CFluidModel *FluidModel, double a, double b, double c);

double rtbis( fct func, double x1, double x2, double xacc, CFluidModel *FluidModel, double a, double b, double c );

int main(int argc, char *argv[]) {

	enum eos_kind {IDEAL_GAS, VW_GAS, PR_GAS, FLP} eos;

    // Open file for reading

	double gamma = 1.12;
	double r = 35.149;
	double pcr = 1415000.0;
	double tcr = 564.09;
	double w= 0.529;
	double P = 1000000;
	double T = 550.00;

	double Mach = 1.5;

//    FILE * inputfile = fopen( "input.inp", "r" );
//    fscanf (inputfile, "%s", eos, "\n");
//    fscanf (inputfile, "%f", &gamma);
//    fscanf (inputfile, "%f", &Mach);

    cout << eos << endl;
    cout << gamma << endl;
    cout << Mach << endl;

	string  thlib = "RefProp";
	string  fluid = "MDM";
	int ncomp = 1;
	double* conc = new double [20];
	conc[0] = 1.0;

	CFluidModel *FluidModel;

	//eos = IDEAL_GAS;
	//eos = VW_GAS;
	//eos = PR_GAS;
	eos = FLP;

	switch ( eos ) {
	case IDEAL_GAS:
		FluidModel = new CIdealGas(gamma, r);
		cout << "FLUID MODEL: IDEAL_GAS" << endl;
		break;
	case VW_GAS:
		FluidModel = new CVanDerWaalsGas(gamma, r, pcr, tcr);
		cout << "FLUID MODEL: VW_GAS" << endl;
		break;
	case PR_GAS:
		FluidModel = new CPengRobinson(gamma, r, pcr, tcr, w);
		cout << "FLUID MODEL: PR_GAS" << endl;
		break;
	case FLP:
		FluidModel = new CFluidProp( thlib, fluid, ncomp, conc );
		cout << "FLUID MODEL: FLUIDPROP LIBRARY" << endl;
		break;
	}

	FluidModel->SetTDState_PT(P,T);
	double ht = FluidModel->GetStaticEnergy() + FluidModel->GetPressure()/FluidModel->GetDensity();
    double s  = FluidModel->GetEntropy();

    double x1 = ht;
    double x2 = ht*0.9;
    cout << "INITIAL GUESS " << endl;
    cout << "x1: " << x1 << endl;
    cout << "x2: " << x2 << endl;

    zbrac( (fct) Energy_hs, &x1, &x2, FluidModel, ht, s, Mach);
    cout << "BRACKETING COMPLETED " << endl;
    cout << "x1: " << x1 << endl;
    cout << "x2: " << x2 << endl;

    double h = rtbis( (fct) Energy_hs, x1, x2, 0.00001, FluidModel, ht, s, Mach );
    cout << "ROOT FOUND " << endl;
    cout << "root: " << h << endl;

    FluidModel->SetTDState_hs(h,s);
	double Pressure = FluidModel->GetPressure();
	double Temperature = FluidModel->GetTemperature();
    cout << "Pressure " << Pressure << endl;
    cout << "Temperature " << Temperature << endl;
    cout << endl;

    SetTotalTDState_prho ( FluidModel, 1542350.0, 225.9312, 1.1 );

}

double Energy_hs ( double h, CFluidModel *FluidModel, double ht, double s, double Mach )
{
	FluidModel->SetTDState_hs(h,s);
	double c   = FluidModel->GetSoundSpeed();

	double f = ht -h -( pow(Mach,2)*pow(c,2) )/2;
	return f;
}

void SetTotalTDState_prho ( CFluidModel *FluidModel, double P, double rho, double Mach )
{
	FluidModel->SetTDState_Prho(P,rho);
	double h = FluidModel->GetStaticEnergy() + FluidModel->GetPressure()/FluidModel->GetDensity();
	double c = FluidModel->GetSoundSpeed();
	double s = FluidModel->GetEntropy();
	double ht = h + ( pow(Mach,2)*pow(c,2) )/2;

	FluidModel->SetTDState_hs(ht,s);
    cout << "Total Pressure " <<  FluidModel->GetPressure() << endl;
    cout << "Total Temperature " << FluidModel->GetTemperature() << endl;
    cout << "Total Density " << FluidModel->GetDensity() << endl;

}


int zbrac( fct func, double *x1, double *x2, CFluidModel *FluidModel, double a, double b, double c)
{
	double FACTOR = 0.2;
	int NTRY = 50;
	int j;
	double f1,f2;

	if (*x1 == *x2) cout << "Bad initial range in zbrac " << endl;
	f1=(*func)(*x1,FluidModel,a,b,c);
	f2=(*func)(*x2,FluidModel,a,b,c);
	for (j=1;j<=NTRY;j++) {
		if (f1*f2 < 0.0) return 1;
		if (fabs(f1) < fabs(f2))
			f1=(*func)(*x1 += FACTOR*(*x1-*x2),FluidModel,a,b,c);
		else
			f2=(*func)(*x2 += FACTOR*(*x2-*x1),FluidModel,a,b,c);
	}
	return 0;
}

double rtbis( fct func, double x1, double x2, double xacc, CFluidModel *FluidModel, double a, double b, double c )
{
	int JMAX = 40;
	int j;
	double dx,f,fmid,xmid,rtb;
	f=(*func)(x1,FluidModel,a,b,c);
	fmid=(*func)(x2,FluidModel,a,b,c);
	if (f*fmid >= 0.0) cout << "Root must be bracketed for bisection in rtbis " << endl;
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	for (j=1;j<=JMAX;j++) {
		fmid=(*func)(xmid=rtb+(dx *= 0.5),FluidModel,a,b,c);
		if (fmid <= 0.0) rtb=xmid;
		if (fabs(dx) < xacc || fmid == 0.0) return rtb;
	}
	cout << "Too many bisections in rtbis " << endl;
	return 0.0;
}










