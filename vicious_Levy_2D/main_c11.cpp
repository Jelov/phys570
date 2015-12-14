//
//  main.cpp
//  vicious Levy flight 2D
//
//  Created by Peng,Cheng-Chieh on 11/25/15.
//  Copyright © 2015 Peng,Cheng-Chieh. All rights reserved.
//


#include <iostream>
#include <fstream>
#include <time.h>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <random>

using namespace std;


long double levyP(long double dx, long double exponent){
	long double answer= powl(dx, -1-exponent);
	return answer;// return pdf value
}

int main(int argc, const char * argv[]) {

default_random_engine r1((unsigned)time(NULL));
default_random_engine r2(r1()+(unsigned)time(NULL));
default_random_engine r3(r2()+(unsigned)time(NULL));
default_random_engine r4(r3()+(unsigned)time(NULL));

	const int dimension=2;
	long double LevyExponent=2;
	//long double exponent=(long double)dimension+LevyExponent;
	long double exponent=LevyExponent;

	int Totaltime= 1000; // total simulate time
	int timeS=0; //current time  in simulation
	int timestep=1;

	long double Xmax=2*powl((long double)Totaltime,1/LevyExponent); // the longest single step

	int tries=2000;
	int diecount=0;
	double dierate=0;
	long double collision_d=0.5;

	// divide tries into part to see fluctuation.
	int part=5;
	int diecountpart[part];
	memset(diecountpart, 0, sizeof(diecountpart));
	double dieratepart[part];
	memset(dieratepart, 0 , sizeof(dieratepart));
	double surratepart[part];
	memset(surratepart, 0 , sizeof(surratepart));
	double deviation[part];
	memset(deviation, 0 , sizeof(deviation));
	double deviation2[part];
	memset(deviation2, 0 , sizeof(deviation2));
	double variance=0;

	int triespart=tries/part;
	int partcount=0;


	// initial position setup
	long double ini_dis=5.0;
	long double ini_angle=1;
	long double x1i=0.0;
	long double x2i=x1i+ini_dis*cos(ini_angle);

	long double y1i=0.0;
	long double y2i=y1i+ini_dis*sin(ini_angle);


	// initialize position for object, p#- predator #kind, n # - #th in this kind, pos[dimension]
	//long double p1n1_pos[dimension];

	// initialize postion for object, p1x[n], p1- creature first kind , x1, x component of n-1th creature
	int np1=1;
	long double p1x[np1];
	memset(p1x, 0, sizeof(p1x));
	long double p1y[np1];
	memset(p1y, 0, sizeof(p1y));

	int np2=1;
	long double p2x[np2];
	memset(p2x, 0, sizeof(p2x));
	long double p2y[np2];
	memset(p2y, 0, sizeof(p2y));




	long double deltar=0;
	long double deltad=0;

	cout<< RAND_MAX<<endl;
	cout<< r1.max()<<endl;
	srand( (unsigned) time(NULL) ); // use time to set rand seed.


	ini_dis=3.0; //survive ~0.75
	ini_angle=1;
	for (LevyExponent=1.25;LevyExponent<3.2;LevyExponent=LevyExponent+0.25){
        cout<<"check 1 , inside Levyexponent loop, LevyExponent = "<<LevyExponent<<endl;

		for(Totaltime=400;Totaltime<7000;Totaltime=Totaltime*2) {
          cout<<"check 2 , inside Total time loop"<<"LevyExponent = "<<LevyExponent<<" Totaltime = "<<Totaltime<<endl;

			ini_angle=1;
			x1i=0.0;
			x2i=x1i+ini_dis;

			//reset counter
			for(int ipart=0;ipart<part;ipart++){
				diecountpart[ipart]=0;
				surratepart[ipart]=0;
				deviation[ipart]=0;
				deviation2[ipart]=0;
			}
			diecount=0;
			variance=0;
			partcount=0;

			//x2i=x1i+ini_dis*cos(ini_angle);


			//y1i=0.0;
			//y2i=y1i+ini_dis*sin(ini_angle);

			//Totaltime= 2000;
			exponent=LevyExponent;
			Xmax=powl((long double)2000,1/LevyExponent); // the longest single step
//			Xmax=powl((long double)Totaltime,1/LevyExponent); // the longest single step

			// for (ini_dis=2.5;ini_dis<3.5;ini_dis=ini_dis+0.5) { // test for different ini_dis
			/*
				 x1i=0.0;
				 x2i=x1i+ini_dis*cos(ini_angle);

				 y1i=0.0;
				 y2i=y1i+ini_dis*sin(ini_angle);

			//reset counter
			for(int ipart=0;ipart<part;ipart++){
			diecountpart[ipart]=0;
			surratepart[ipart]=0;
			deviation[ipart]=0;
			deviation2[ipart]=0;
			}
			diecount=0;
			variance=0;
			partcount=0;


			Totaltime= 2000;
			Xmax=2*powl((long double)Totaltime,1/LevyExponent); // the longest single step
			*/    

			for (int trycount=1;trycount<(tries+1);trycount++){
				// for each try ,reset to the initial condition (position and time)
				long double distance_temp=0; // the closet distance between different group members.
				// initialize all members of group to the same initial position
				for (int ip1=0;ip1<np1;ip1++){
					p1x[ip1]=x1i;
					p1y[ip1]=y1i;
				}
				for (int ip2=0;ip2<np2;ip2++){
					p2x[ip2]=x2i;
					p2y[ip2]=y2i;
				}
				timeS=0;

				for  (timeS=0;timeS<Totaltime;timeS+=timestep){

					// update p1, first kind
					deltar=(Xmax-1)*(long double)r1()/(long double)r1.max()+1;
					while ( levyP(deltar,exponent)<((long double)r1()/r1.max())){
						deltar=(Xmax-1)*(long double)r1()/(long double)r1.max()+1;
					}
					deltad=M_PI*(long double)r2()/(long double)r2.max();

					p1x[0]+=deltar*cos(deltad);
					p1y[0]+=deltar*sin(deltad);

					//calculate distance
					distance_temp=sqrt(pow((p1x[0]-p2x[0]),2)+pow((p1y[0]-p2y[0]),2));
					//cout<<distance_temp;
					if (distance_temp<collision_d){
						diecountpart[partcount]++;
						break;
					}

					// update p2, second kind
					deltar=(Xmax-1)*(long double)r3()/(long double)r3.max()+1;
					while ( levyP(deltar,exponent)<((long double)r3()/r3.max())){
						deltar=(Xmax-1)*(long double)r3()/(long double)r3.max()+1;
					}
					deltad=M_PI*(long double)r4()/(long double)r4();

					p2x[0]+=deltar*cos(deltad);
					p2y[0]+=deltar*sin(deltad);

					//calculate distance
					distance_temp=sqrt(pow((p1x[0]-p2x[0]),2)+pow((p1y[0]-p2y[0]),2));
					if (distance_temp<collision_d){
						diecountpart[partcount]++;
						break;
					}




				}

				if ((trycount%triespart)==0){
					diecount+=diecountpart[partcount];
					//cout<<"diecountpart = "<<diecountpart[partcount]<<endl;
					dieratepart[partcount]=(double)diecountpart[partcount]/(double)triespart;
					surratepart[partcount]=1-dieratepart[partcount];
					cout<<"trycount = "<<trycount<<"  surive rate = "<< surratepart[partcount]<<endl;
					partcount++;

				}
			}

			// calculate and print out

			dierate=(double)diecount/(double)tries;
			cout<<"die rate= "<< dierate<<"  diecount = "<<diecount<<endl;
			double sur_rate=1-dierate;
			cout<<"sur_rate = "<<sur_rate<<endl;
			variance=0;
			for(int icount=0; icount<part; icount++)
			{
				deviation[icount]=fabs( (double)surratepart[icount]-(double)sur_rate);///(double)sur_rate);
				deviation2[icount]=pow(deviation[icount],2);
				//  cout<<"error "<<fabs( (double)surratepart[icount]-(double)sur_rate)<< "  deviation = "<<deviation[icount]<<endl;
				//  cout<<"survive rate "<<surratepart[icount]<<" part = "<<icount<<"  deviation "<<deviation[icount]<<endl;
				if (deviation[icount]>0.1) cout<<"warning deviation too big increase time or MC tries"<<endl;
				variance+=deviation2[icount]/part;
			}

			double Std=sqrt(variance);
			double StdoverMean=Std/sur_rate;
			if (Std>0.04||StdoverMean>0.05) {cout<<"warning:: std (relative error) too large"<<endl; ;}
			//cout<<"diecount = "<<diecount<<endl;
			//cout<<"die rate = "<<dierate<<endl;
			//cout<<"surive rate = "<< 1-dierate<<endl;


			//writing result in a file
			cout<<"dimesion ; "<<"levy exponent ; "<<"simulations ; "<<"Timesteps ; "<<" ini_distance;"<<"survive rate ; "<<"Std;"<<"Std over mean;"<<endl;
			cout<<dimension<<";"<<LevyExponent<<";"<<tries<<";"<<Totaltime<<";"<<ini_dis<<";"<<sur_rate<<";"<<Std<<";"<<StdoverMean<<";\n"<<endl;

			char filename[]="viciouswalk2D.txt";
			fstream myfile;
			myfile.open (filename,ios::out|ios::app);
			//myfile << "vicious walk 2D"<<endl;;
			//myfile<<"dimesion ; "<<"levy exponent ; "<<"simulations ;"<<"Timesteps ; "<<"inidistance ;"<<"survive rate ; "<<"Std;"<<"Std over mean;"<<endl;
			myfile<<dimension<<";"<<LevyExponent<<";"<<tries<<";"<<Totaltime<<";"<<ini_dis<<";"<<sur_rate<<";"<<Std<<";"<<StdoverMean<<";"<<endl;
			myfile.close();
		//}
		}
		}
		cout <<"runtime = "<< (double)clock() / CLOCKS_PER_SEC << " S"<<endl;
		return 0;
	}
