//
//  random_2D.cpp
//  random walk 2D
//
//  Created by Peng,Cheng-Chieh on 12/3/15.
//  Copyright (c) 2015 Peng,Cheng-Chieh. All rights reserved.
//
//  Create a number based simple 2-D random walk

#include <iostream>
#include <fstream>
#include <time.h>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

using namespace std;

const long double PI = 3.141592653589793238462643383279502884;
long double pi=atan(1)*4;
int random_2D(int argc, const char * argv[]) {
    cout<<"pi = "<<pi<<endl;
    // creat particle position
    long double ini_dis=3.0;
    long double ini_angle=1;
    long double x1i=0.0;
    long double x2i=x1i+ini_dis*cos(ini_angle);
    long double x1=0.0;
    long double x2=0.0;
    
    long double y1i=0.0;
    long double y2i=y1i+ini_dis*sin(ini_angle);
    long double y1=0.0;
    long double y2=0.0;

    
    int tries=40;
    int diecount=0;
    long double deltar1=0;
    long double deltar2=0;
    long double deltad1=0;
    long double deltad2=0;

    cout<< RAND_MAX<<endl;;
    
    srand( (unsigned) time(NULL) ); // use time to set rand seed.
    
    // set
    long double Totaltime= 500; // total simulate time
    long double timeS=0; //simulation current time
    long double timestep=1;
    /*
    for (int trycount=0;trycount<tries;trycount++){
        //reset the initial position
        x1=x1i;
        x2=x2i;
        y1=y1i;
        y2=y2i;
        timeS=0;
        while (timeS<Totaltime){
            
            deltar1=(long double)rand()/(long double)RAND_MAX;
            deltar2=(long double)rand()/(long double)RAND_MAX;
            deltad1=(long double)rand()/(long double)RAND_MAX;
            deltad2=(long double)rand()/(long double)RAND_MAX;
            if (deltax1D<=0.5)  x1=x1+1;
            if (deltax1D>0.5)   x1=x1-1;
            if (x1==x2) {    // calculate if they meet, check before update x2 for if they cross
                cout<<"meet"<<endl;
                diecount++;
                break;
            }
            deltax2D=(long double)rand()/(long double)RAND_MAX;
            if (deltax2D<=0.5) x2=x2+1;
            if (deltax2D>0.5) x2=x2-1;
            
            if (x1==x2) {    // calculate if they meet
                cout<<"meet"<<endl;
                diecount++;
                break;
            }
            
            timeS=timeS+timestep;
        }
        cout<<"trycount= "<<trycount<<endl;
        cout<<"distance= "<<x2-x1<<endl;
    }*/
    cout<<"diecount = "<<diecount<<endl;
    double die_rate=(double)diecount/(double)tries;
    cout<<"die rate = "<<die_rate<<endl;
    return 0;
}
