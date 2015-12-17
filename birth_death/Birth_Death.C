//
//  Simple birth_death process simulation 
//
//  Created by Peng,Cheng-Chieh on 12/7/15.
//  Copyright © 2015 Peng,Cheng-Chieh. All rights reserved.
//
// 1. X-> X+1 ,   a[0]=kb ,     v[0]=1
// 2. X-> X-1 ,   a[1]=kd*X,    v[1]=-1


#include "TGraph.h"  //root library for plot
#include "TCanvas.h" 
#include "TAxis.h"
#include "TLatex.h"
#include "TTree.h"
#include "TH1D.h"
#include "TFile.h"
#include "TBranch.h"
#include "TMath.h"
#include "TF1.h"

#include <iostream>
#include <fstream>
#include <time.h>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

//#include <random> //c++11 random


using namespace std;

int factorial(int x, int result = 1) {
	if (x == 1) return result; else return factorial(x - 1, x * result);
}

void Birth_Death(){
	TFile *fout = new TFile("myOutput.root","recreate"); fout->cd();
  TTree *t0 = new TTree("t0","t0");

	
	// create two different random value with two different seed.
	srand( (unsigned) time(NULL) ); // use time to set rand seed.
	cout<<"RAND_MAX="<< RAND_MAX<<endl;
	//    default_random_engine r1((unsigned)time(NULL));
	//    default_random_engine r2(rand());
	//    cout<<"r1max = "<<r1.max()<<endl;
	double r2t=0;
	double ajsum=0;

	//
	const int N_system = 2000; //final result use 2000+
	const double endtime = 200;
	const int N_reaction = 2;

	double propensity[N_reaction];
	double propensity0; // Sum over all prepensity, a0
	double CSV[N_reaction]={1,-1};   // change of state vector

	double kb=2;
	double kd=0.05;

	// initail parameter
	int x_ini=0;
	double t=0;


	int timepart=1; // the time step for recording data.
	int nparts=endtime/timepart+1;

	double x_temp[N_system][nparts];
	double x2_temp[N_system][nparts];
	double x_mean[nparts];
	double x2_mean[nparts];
	double variance[nparts];
	double variance2[nparts];
	double skewness[nparts];

	int fpflag=0; //first-passage time flag
	const int fpthreshold=0.5*kb/kd; // set to half of the eq # with enough time, all sys should pass threshold
	double fptime[N_system]; 

	int x[N_system];
	for (int isys = 0; isys<N_system ; isys++){
		x[isys]=x_ini;
		fptime[isys]=-1; 
		for (int time=0; time<nparts;time++){
			x_temp[isys][time]=0;
			x2_temp[isys][time]=0;
		}

	}
	for (int iparts=0; iparts<nparts;iparts++){
		x_mean[iparts]=0;
		x2_mean[iparts]=0;
		variance[iparts]=0;
		variance2[iparts]=0;
		skewness[iparts]=0;
	}



	propensity[0]=kb;
	propensity[1]=kd*x[0];
	propensity0=propensity[0]+propensity[1];

	double dt=0;
	int i_re=0; //index to sum over ai, determine which reaction.

	for (int isys=0; isys<N_system ;isys++){
		t=0;
		x[isys]=x_ini;
		x_temp[isys][0]=0; 
		fpflag=0; //reset fpflag to 0;
		while (t<endtime){

			propensity[0]=kb;
			propensity[1]=kd*x[isys];
			propensity0=propensity[0]+propensity[1];
			dt=(1/propensity0)*log((double)RAND_MAX/(double)rand()); 
			r2t=(double)rand()/(double)RAND_MAX;

			ajsum=0;
			i_re=0;
			for (i_re=0;ajsum<r2t*propensity0;i_re++){
				ajsum+=propensity[i_re];
				if(i_re>=N_reaction){
					cout<<"error 1 in determine the reaction, terminate"<<endl;
				}
			}
			i_re=i_re-1; // the for loop would add extra one to i_re
			x[isys]=x[isys]+CSV[i_re];
			t=t+dt;
			if(x[isys]<0){
				cout<<"x<0"<<endl;

			}

			int iparts=t/timepart;
			if ((int)floor(t)%timepart==0 && x_temp[isys][iparts]==0 &&iparts>0){
				x_temp[isys][iparts]=x[isys];
				x2_temp[isys][iparts]=pow(x[isys],2);
				x_mean[iparts]=x_mean[iparts]+x_temp[isys][iparts]/N_system;
				x2_mean[iparts]=x2_mean[iparts]+x2_temp[isys][iparts]/N_system;

			}
			if (x[isys]>fpthreshold && fpflag==0){
			fpflag=1;
			fptime[isys]=t;			
			}


		} // end while time
	} // end for isys


/////////////////////////////////
// First-Passage time
/////////////////////////////
	double tfp=0;
  t0->Branch("tfp",&tfp,"tfp/D");
	for (int isys=0; isys<N_system;isys++){
	tfp=fptime[isys];
	t0->Fill();
	}
	t0->Write();

	// analytical calculation for first passage time pdf
double fp_pdf[100]; 
double fp_t[100];
for (int i=0;i<100;i++){
fp_pdf[i]=0;
fp_t[i]=i;
}

for(int it=0 ; it<100 ; it++){
double lambda_temp=kb/kd*(1-exp(-kd*it));
for(int ix=0 ; ix<=fpthreshold ; ix++){ 
	if (ix ==0) {fp_pdf[it]-=-exp(-lambda_temp)*kb*exp(-kd*it);   }
	if (ix !=0) {
	fp_pdf[it]-= 1/TMath::Factorial(ix) * exp(-lambda_temp) * (ix*pow(lambda_temp,ix-1)-pow(lambda_temp,ix)) * kb*exp(-kd*it);
	}

	}	// end for ix
} // end for it

  TCanvas *cfp= new TCanvas("cfp","fp");
  TGraph *grfp= new TGraph(50,fp_t,fp_pdf);
  grfp->GetXaxis()->SetTitle("t");
  grfp->GetYaxis()->SetTitle("Probability");
  grfp->GetXaxis()->CenterTitle();
  grfp->GetXaxis()->CenterTitle();
  grfp->SetLineColor(1);
  grfp->SetMarkerColor(1);
  grfp->Draw("AC*");
  grfp->SetTitle("fp");
	cfp->SaveAs("cfp.pdf");



	// calculate variance directly from (X-mean)^2
	for(int iparts=0; iparts<nparts; iparts++){
		for (int isys=0;isys<N_system;isys++){
			variance2[iparts]=variance2[iparts]+pow((x_temp[isys][iparts]-x_mean[iparts]),2)/N_system;
		}
	}
	// calculate skewness
	for(int iparts=0; iparts<nparts; iparts++){
		for (int isys=0;isys<N_system;isys++){
			skewness[iparts]=skewness[iparts] +	pow(((x_temp[isys][iparts]-x_mean[iparts])/sqrt(variance2[iparts])),3)/N_system;	

		}	
	}    
	for (int iparts=1; iparts<nparts;iparts++){
		variance[iparts]=x2_mean[iparts]-pow((x_mean[iparts]),2); // another way to calculate variance, double check
	}


	// analytical calculation
	// mean x= kb/kd[1-e^-kdt]

	double xmean_ana[nparts];
	for (int iparts=0;iparts<nparts;iparts++){
		xmean_ana[iparts]=kb/kd*(1-exp(-kd*iparts*timepart));
	}



	////////////////////////////////////////////////
	// -- two-time auto-correlation function c(t,t1)=<x(t+t1)x(t)>-<x(t+t1)><x(t)>
	// -- analytical c(t,t1)=variance(t)*exp((-kd*t1)
	////////////////////////////////////////////////

	int timestep_cxx=1;
	int nstep_cxx=40;
	int ini_time_cxx=10;
	double Cxx_MC[nstep_cxx];
	double Cxx_ana[nstep_cxx];
	double timearray_cxx[nstep_cxx];
	for (int istep=0 ; istep<nstep_cxx; istep++){
		Cxx_MC[istep]=0;

		for (int isys=0 ; isys<N_system; isys++){
			Cxx_MC[istep]+=x_temp[isys][ini_time_cxx]*x_temp[isys][ini_time_cxx+timestep_cxx*istep]/N_system;

		}//end for isys
		Cxx_MC[istep]=Cxx_MC[istep]-x_mean[ini_time_cxx]*x_mean[ini_time_cxx+timestep_cxx*istep];
		timearray_cxx[istep]=timestep_cxx*istep;
		Cxx_ana[istep]=variance2[ini_time_cxx]*exp(-kd*timestep_cxx*istep);
	}//end for istep
	TCanvas *cxx= new TCanvas("cx","two time auto correlation");
	cxx->Divide(1,2);
	TGraph *grcxx= new TGraph(nstep_cxx,timearray_cxx,Cxx_MC);
	TGraph *grcxxa= new TGraph(nstep_cxx,timearray_cxx,Cxx_ana);

	cxx->cd(1);
	grcxxa->GetXaxis()->SetTitle("delta Time from t=10");
	grcxxa->GetYaxis()->SetTitle("Cxx(t,t1)");
	grcxxa->GetXaxis()->CenterTitle();
	grcxxa->GetXaxis()->CenterTitle();
	grcxxa->SetLineColor(4);
	grcxxa->SetMarkerColor(4);
	grcxxa->Draw("AC*");
	grcxxa->SetTitle("two time auto correlation");

	cxx->cd(2);
	grcxx->SetMarkerStyle(21);
	grcxx->SetMarkerColor(1);
	grcxx->SetMarkerSize(0.8);
	grcxx->Draw("AP*");
	cxx->SaveAs("cxx.pdf");


	////////////////////////////////////////////////
	//	 	for plot, using root library to plot 		//
	////////////////////////////////////////////////

	double timearray[nparts];
	for (int iparts=0; iparts<nparts; iparts++){
		timearray[iparts]=iparts*timepart;
	}

	TCanvas *c1= new TCanvas("c1","X_mean");
	TGraph *gr1= new TGraph(nparts,timearray,xmean_ana);
	TGraph *gr2= new TGraph(nparts,timearray,x_mean);
	gr1->GetXaxis()->SetTitle("Time");
	gr1->GetYaxis()->SetTitle("X");
	gr1->GetXaxis()->CenterTitle();
	gr1->GetXaxis()->CenterTitle();

	gr1->SetLineColor(4);
	gr1->SetMarkerColor(4);
	gr1->Draw("AC*");
	gr1->SetTitle("X-Mean");

	gr2->SetMarkerStyle(21);
	gr2->SetMarkerColor(1);
	gr2->SetMarkerSize(0.8);
	gr2->Draw("P");

	c1->SaveAs("X_mean.pdf");

	TCanvas *c2= new TCanvas("c2","Variance");
	TGraph *gr3= new TGraph(nparts,timearray,variance2);
	gr3->GetXaxis()->SetTitle("Time");
	gr3->GetYaxis()->SetTitle("Variance");
	gr3->GetXaxis()->CenterTitle();
	gr3->GetXaxis()->CenterTitle();
	gr3->SetLineColor(1);
	gr3->SetMarkerColor(1);
	gr3->Draw("AC*");
	gr3->SetTitle("Variance");
  c2->SaveAs("Variance.pdf");


	
		 TCanvas *c3= new TCanvas("c3","Skewness");
		 TGraph *gr4= new TGraph(nparts,timearray,skewness);
		 gr4->GetXaxis()->SetTitle("Time");
		 gr4->GetYaxis()->SetTitle("Skewness");
		 gr4->GetXaxis()->CenterTitle();
		 gr4->GetXaxis()->CenterTitle();
	gr4->SetLineColor(1);
	gr4->SetMarkerColor(1);
	gr4->Draw("AC*");
	gr4->SetTitle("Skewness");
  c3->SaveAs("Skewness.pdf");

	int xmin = 0;
	int xmax = (int)(kb/kd*2-1);
	int xbins = kb/kd;
	int xpdf_temp1=0;
	int xpdf_temp2=0;
	int xpdf_temp3=0;
	int xpdf_temp4=0;
	int xpdf_temp5=0;
	int pdf_parts=5;


//	TFile *fout = new TFile("myOutput.root","recreate"); fout->cd();
	TTree *t1 = new TTree("t1","t1");
	TH1D *h_xdata1= new TH1D("h_xdata1","h_xdata1",xbins,xmin,xmax);
	TH1D *h_xdata2= new TH1D("h_xdata2","h_xdata2",xbins,xmin,xmax);
	TH1D *h_xdata3= new TH1D("h_xdata3","h_xdata3",xbins,xmin,xmax);
	TH1D *h_xdata4= new TH1D("h_xdata4","h_xdata4",xbins,xmin,xmax);
	TH1D *h_xdata5= new TH1D("h_xdata5","h_xdata5",xbins,xmin,xmax);

	t1->Branch("xpdf1_data",&xpdf_temp1,"xpdf_temp1/I");
	t1->Branch("xpdf2_data",&xpdf_temp2,"xpdf_temp2/I");
	t1->Branch("xpdf3_data",&xpdf_temp3,"xpdf_temp3/I");
	t1->Branch("xpdf4_data",&xpdf_temp4,"xpdf_temp4/I");
	t1->Branch("xpdf5_data",&xpdf_temp5,"xpdf_temp5/I");

	int itimepart=10;
	int timesteps_pdf=timepart*itimepart;
	for (int isys=0; isys<N_system ; isys++){
		xpdf_temp1=x_temp[isys][itimepart*1];
		xpdf_temp2=x_temp[isys][itimepart*2];
		xpdf_temp3=x_temp[isys][itimepart*3];
		xpdf_temp4=x_temp[isys][itimepart*4];
		xpdf_temp5=x_temp[isys][itimepart*5];
		t1->Fill();
	}
	t1->Write();

	//analytical prediction for pdf P(x,t)=exp(-lambda(t))* lambda(t)^x / x! , lambda(t)= kb/kd(1-exp(-kdt))
	double lambda_t[5];
	double pdf1[xmax];
  double pdf2[xmax];
  double pdf3[xmax];
  double pdf4[xmax];
  double pdf5[xmax];

	double xinput[xmax];
	double pdf1sum=0;
  double pdf2sum=0;
  double pdf3sum=0;
  double pdf4sum=0;
  double pdf5sum=0;


	for (int il=0;il<5;il++){
		lambda_t[il]=kb/kd*(1-exp(-kd*timesteps_pdf*(il+1)));
	}
	for (int xin=0; xin<xmax ; xin++){ //calculate analytical value manually 
		xinput[xin]=xin;
		pdf1[xin]=(exp(-lambda_t[0])*pow(lambda_t[0],xin))/(TMath::Factorial(xin)); 
    pdf2[xin]=(exp(-lambda_t[1])*pow(lambda_t[1],xin))/(TMath::Factorial(xin));
    pdf3[xin]=(exp(-lambda_t[2])*pow(lambda_t[2],xin))/(TMath::Factorial(xin));
    pdf4[xin]=(exp(-lambda_t[3])*pow(lambda_t[3],xin))/(TMath::Factorial(xin));
    pdf5[xin]=(exp(-lambda_t[4])*pow(lambda_t[4],xin))/(TMath::Factorial(xin));

		pdf1sum+=pdf1[xin];
    pdf2sum+=pdf2[xin];
    pdf3sum+=pdf3[xin];
    pdf4sum+=pdf4[xin];
    pdf5sum+=pdf5[xin];

	}
	cout<<"pdf1sum = "<<pdf1sum<<endl;
	for (int xin=0; xin<xmax ; xin++){
		pdf1[xin]=pdf1[xin]/pdf1sum;
    pdf2[xin]=pdf2[xin]/pdf2sum;
    pdf3[xin]=pdf3[xin]/pdf3sum;
    pdf4[xin]=pdf4[xin]/pdf4sum;
    pdf5[xin]=pdf5[xin]/pdf5sum;

	}



	TCanvas *cpdf1= new TCanvas("cpdf1","pdf1");
	TGraph *grpdf1= new TGraph(xmax,xinput,pdf1);
	grpdf1->GetXaxis()->SetTitle("X");
	grpdf1->GetYaxis()->SetTitle("Probability");
	grpdf1->GetXaxis()->CenterTitle();
	grpdf1->GetXaxis()->CenterTitle();
	grpdf1->SetLineColor(2);
	grpdf1->SetMarkerColor(2);
	grpdf1->Draw("AC*");
	grpdf1->SetTitle("Pdf1 (time=10)");
	cpdf1->SaveAs("cpdf1.pdf");

  TCanvas *cpdf2= new TCanvas("cpdf2","pdf2");
  TGraph *grpdf2= new TGraph(xmax,xinput,pdf2);
  grpdf2->GetXaxis()->SetTitle("X");
  grpdf2->GetYaxis()->SetTitle("Probability");
  grpdf2->GetXaxis()->CenterTitle();
  grpdf2->GetXaxis()->CenterTitle();
  grpdf2->SetLineColor(2);
  grpdf2->SetMarkerColor(2);
  grpdf2->Draw("AC*");
  grpdf2->SetTitle("Pdf2 (time=20)");
  cpdf2->SaveAs("cpdf2.pdf");

  TCanvas *cpdf3= new TCanvas("cpdf3","pdf3");
  TGraph *grpdf3= new TGraph(xmax,xinput,pdf3);
  grpdf3->GetXaxis()->SetTitle("X");
  grpdf3->GetYaxis()->SetTitle("Probability");
  grpdf3->GetXaxis()->CenterTitle();
  grpdf3->GetXaxis()->CenterTitle();
  grpdf3->SetLineColor(2);
  grpdf3->SetMarkerColor(2);
  grpdf3->Draw("AC*");
  grpdf3->SetTitle("Pdf3 (time=30)");
  cpdf3->SaveAs("cpdf3.pdf");

  TCanvas *cpdf4= new TCanvas("cpdf4","pdf4");
  TGraph *grpdf4= new TGraph(xmax,xinput,pdf4);
  grpdf4->GetXaxis()->SetTitle("X");
  grpdf4->GetYaxis()->SetTitle("Probability");
  grpdf4->GetXaxis()->CenterTitle();
  grpdf4->GetXaxis()->CenterTitle();
  grpdf4->SetLineColor(2);
  grpdf4->SetMarkerColor(2);
  grpdf4->Draw("AC*");
  grpdf4->SetTitle("Pdf4 (time=40)");
  cpdf4->SaveAs("cpdf4.pdf");

  TCanvas *cpdf5= new TCanvas("cpdf5","pdf5");
  TGraph *grpdf5= new TGraph(xmax,xinput,pdf5);
  grpdf5->GetXaxis()->SetTitle("X");
  grpdf5->GetYaxis()->SetTitle("Probability");
  grpdf5->GetXaxis()->CenterTitle();
  grpdf5->GetXaxis()->CenterTitle();
  grpdf5->SetLineColor(2);
  grpdf5->SetMarkerColor(2);
  grpdf5->Draw("AC*");
  grpdf5->SetTitle("Pdf5 (time=50)");
  cpdf5->SaveAs("cpdf5.pdf");

/*	TF1 *f1 = new TF1("f1",Form("[0]*(exp(-%f)*pow(%f,x))/TMath::Factorial(x)",lambda_t[0],lambda_t[0]),0,100);
	f1->SetParameter(0.001,1000);
	TCanvas *ctest = new TCanvas("ctest","ctest");
	f1->Draw();
*/

	TCanvas *c4= new TCanvas("c4","xpdf");
	c4->Divide(3,2);
	c4->cd(1);
	t1->Draw("xpdf1_data>>h_xdata1"); // could plot with weight = 1/N_system to normalize
	c4->cd(2);
	t1->Draw("xpdf2_data>>h_xdata2");
	c4->cd(3);
	t1->Draw("xpdf3_data>>h_xdata3");
	c4->cd(4);
	t1->Draw("xpdf4_data>>h_xdata4");
	c4->cd(5);
	t1->Draw("xpdf5_data>>h_xdata5");

	TCanvas *c5= new TCanvas("c5","xpdf1");
	h_xdata1->Fit("f1");
	h_xdata1->Draw();


	////////////////////////////////////////////////
	//  Time average v.s ensemble average
	//
	////////////////////////////////////////////////

	double endtime1=5000;
	double xmean_tave=0;
	double x2mean_tave=0;
	double xtime_ave[(int)endtime1],x2time_ave[(int)endtime1];
	for(int it=0; it<endtime ; it++){
	xtime_ave[it]=0;
	x2time_ave[it]=0;
	} 
	int tave_count=0; // count how many entry used.
		t=0;
		x[0]=x_ini;
	timepart=1/kd;
	cout<<"timepart = "<<timepart<<endl;
		while (t<endtime1){
			propensity[0]=kb;
			propensity[1]=kd*x[0];
			propensity0=propensity[0]+propensity[1];
			dt=(1/propensity0)*log((double)RAND_MAX/(double)rand());
			r2t=(double)rand()/(double)RAND_MAX;
			ajsum=0;
			i_re=0;
			for (i_re=0;ajsum<r2t*propensity0;i_re++){
				ajsum+=propensity[i_re];
				if(i_re>=N_reaction){
					cout<<"error 1 in determine the reaction, terminate"<<endl;
				}
			}
			i_re=i_re-1; // the for loop would add extra one to i_re
			x[0]=x[0]+CSV[i_re];
			t=t+dt;
			if(x[0]<0){
				cout<<"x<0"<<endl;
			}
			int iparts=t/timepart;
			if ((int)floor(t)%timepart==0 && xtime_ave[iparts]==0 &&iparts>(10/kd)){ 
				xtime_ave[iparts]=x[0];
				x2time_ave[iparts]=pow(x[0],2);
				xmean_tave=xmean_tave+xtime_ave[iparts];
        x2mean_tave=x2mean_tave+x2time_ave[iparts];
				tave_count++;
			}
		}// end of while time

	xmean_tave=xmean_tave/tave_count;
	x2mean_tave=x2mean_tave/tave_count;
	double variance_tave=x2mean_tave-xmean_tave*xmean_tave;
	cout<<"xmean_tave= "<< xmean_tave<<" ,x2mean_tave = "<<x2mean_tave<<" , variance_tave = "<<variance_tave<<endl;
	cout<<"xmean_ens= "<< x_mean[nparts-1]<<" ,x2mean_ens = "<<x2_mean[nparts-1]<<" , variance_ens = "<<variance[nparts-1]<<endl;




	std::cout << "Hello, World!\n";
	fout->Close();
}

