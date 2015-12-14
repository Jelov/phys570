//
//  main.cpp
//  vicious Levy flight 1D
//
//  Created by Peng,Cheng-Chieh on 11/18/15.
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

using namespace std;


long double levyP(long double dx, long double exponent){
	long double answer= powl(dx, -1-exponent);
	return answer;// return pdf value
}

int main(int argc, const char * argv[]) {

	const int dimension=1;
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
	long double x2i=0.0; //x1i+ini_dis*cos(ini_angle);

	//long double y1i=0.0;
	//long double y2i=0.0; //y1i+ini_dis*sin(ini_angle);


	// initialize position for object, p#- predator #kind, n # - #th in this kind, pos[dimension]
	//long double p1n1_pos[dimension];

	// initialize postion for object, p1x[n], p1- creature first kind , x1, x component of n-1th creature
	int np1=1;
	long double p1x[np1];
	memset(p1x, 0, sizeof(p1x));
	//   long double p1y[np1];
	//   memset(p1y, y1i, sizeof(p1y));

	int np2=1;
	long double p2x[np2];
	memset(p2x, 0, sizeof(p2x));
	//   long double p2y[np2];
	//   memset(p2y, y2i, sizeof(p2y));


	long double deltar=0;
	// long double deltad=0;

	cout<< RAND_MAX<<endl;;

	srand( (unsigned) time(NULL) ); // use time to set rand seed.


	ini_dis=40.0;
    //check time at 200,500,1000,1500,2000,2500,3000,3500,4000

	// for (ini_dis=16;ini_dis<40;ini_dis=ini_dis+8){ // test for different ini_dis
	for (LevyExponent=0.75;LevyExponent<0.9;LevyExponent=LevyExponent+0.25){
        cout<<"check 1 , inside Levyexponent loop"<<endl;
		for(Totaltime=100;Totaltime<4000;Totaltime=Totaltime*2) {
            cout<<"check 2 , inside Total time loop"<<endl;
            cout<<"LevyExponent = "<<LevyExponent<<" Totaltime = "<<Totaltime<<endl;
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
			Xmax=powl((long double)2000,1/LevyExponent);
//			Xmax=2*powl((long double)Totaltime,1/LevyExponent); // the longest single step


			for (int trycount=1;trycount<(tries+1);trycount++){
                //cout<<"chekc 3, inside trycount loop"<<endl;
				// for each try ,reset to the initial condition (position and time)
				long double distance_temp=0; // the closet distance between different group members.
				for (int ip1=0;ip1<np1;ip1++){
					p1x[ip1]=x1i;
				}
				for (int ip2=0;ip2<np2;ip2++){
					p2x[ip2]=x2i;
				}
				timeS=0;

				for  (timeS=0;timeS<Totaltime;timeS+=timestep){

					// update p1, first kind
					deltar=(Xmax-1)*(long double)rand()/(long double)RAND_MAX+1;
					while ( levyP(deltar,exponent)<((long double)random()/RAND_MAX)){
						deltar=(Xmax-1)*(long double)rand()/(long double)RAND_MAX+1;
					}
					//deltad=M_PI*(long double)rand()/(long double)RAND_MAX;

					p1x[0]+=deltar;
					//p1x[0]+=deltar*cos(deltad);
					//p1y[0]+=deltar*sin(deltad);

					//calculate distance
					distance_temp=fabsl(p1x[0]-p2x[0]);
					//cout<<distance_temp;
					if (distance_temp<collision_d){
						diecountpart[partcount]++;
						break;
					}


					// update p2, second kind
					deltar=(Xmax-1)*(long double)rand()/(long double)RAND_MAX+1;
					while ( levyP(deltar,exponent)<((long double)random()/RAND_MAX)){
						deltar=(Xmax-1)*(long double)rand()/(long double)RAND_MAX+1;
					}
					//deltad=M_PI*(long double)rand()/(long double)RAND_MAX;

					p2x[0]+=deltar;
					//p1x[0]+=deltar*cos(deltad);
					//p1y[0]+=deltar*sin(deltad);

					//calculate distance
					distance_temp=fabsl(p1x[0]-p2x[0]);
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

			char filename[]="viciouswalk1D.txt";
			fstream myfile;
			myfile.open (filename,ios::out|ios::app);
			//myfile << "vicious walk 1D"<<endl;;
			//myfile<<"dimesion ; "<<"levy exponent ; "<<"simulations ;"<<"Timesteps ; "<<"inidistance ;"<<"survive rate ; "<<"Std;"<<"Std over mean;"<<endl;
			myfile<<dimension<<";"<<LevyExponent<<";"<<tries<<";"<<Totaltime<<";"<<ini_dis<<";"<<sur_rate<<";"<<Std<<";"<<StdoverMean<<";"<<endl;
			myfile.close();
	//}
		}
	}

	cout <<"runtime = "<< (double)clock() / CLOCKS_PER_SEC << " S"<<endl;
	return 0;
}

