//#define CLUSTER
#define PI 3.14159
#include <cstdlib>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <sstream>
#include "Langevin_dynamics.h"
#ifdef CLUSTER
#include "/home/svaikunt/local/lib/include/fftw3.h"
#endif 
#ifndef CLUSTER
#include <fftw3.h>
#endif
#include <assert.h>
using namespace std;

Langevin_dynamics::~Langevin_dynamics(){
	gsl_rng_free(r);
}

Langevin_dynamics::Langevin_dynamics(){}



void Langevin_dynamics::initialize(int N_max1,int N_max2,long int randomseed1,double S1,int inputsoft,double ksolute_solvent){
    RCUT=1.0;
    ksoft=1.0;
    k12=ksolute_solvent;
    epsilon=ksolute_solvent;//WARNING, REUSED PARAMETER.
    soft=inputsoft;
    asoft=1.0;
    Lx=4.0; //size of box
    Ly=4.0; //size of box
    N_type=2;// number of particle types
    Ntype[0]=N_max2; //There are two partivles
    Ntype[1]=0;//Number of particles of type 1 input by user set to zero.
	S=S1; //value of S for biasing.
	N_max=N_max1;
    randomseed=randomseed1;
	r=gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(r,randomseed);
    for (int loopi=0;loopi<N_type;loopi++){
        vector< vector<double> > pos2;
        vector< vector<double> >fpos2;
        for (int loopj=0;loopj<Ntype[loopi];loopj++){
            vector<double> pos1;
            vector<double> fpos1;
            for(int loopk=0;loopk<2;loopk++){
                pos1.push_back(gsl_rng_uniform(r)*Lx);
                fpos1.push_back(0);
            }
            pos2.push_back(pos1);
            fpos2.push_back(fpos1);
        }
        pos.push_back(pos2);
        fpos.push_back(fpos2);
        fpostype.push_back(fpos2);
        upostype.push_back(fpos2);
        zerovec.push_back(fpos2);
    }
    for (int loopi=0;loopi<N_type;loopi++){
            vector<double> pos1;
            vector<double> fpos1;
            for(int loopj=0;loopj<Ntype[loopi];loopj++){
                pos1.push_back(gsl_rng_uniform(r)*Lx);
                fpos1.push_back(0);
            }
            ftheta.push_back(fpos1);
            theta.push_back(pos1);
	    zerovectheta.push_back(fpos1);
        }

    //equilibrate();//initial equilibration run.
    //pos[type][Number][dimension];
}

void Langevin_dynamics::equilibrate(double Pe, double tau){
    double fd_term,noise_0,noise_1;
    double del1,del2,deltheta1;
    gamma_i=2;
    int typebin=0;
    //Compute forces, evolve dynamics, return y factor
    double time=0;
    double dt=0.001;
    int totaltimeint=dt*10000;
    int count=0;
    char outputfile[100];
    sprintf(outputfile,"EquilibriumN%d.S%.3f.XYZ",N_max,S);
    ofstream fileout;
    fileout.open(outputfile);
    for (time=0;time<totaltimeint;time=time+dt){
        compute_forces_theta();
        count+=1;
       
        if (count%1000==0){
           fileout<<N_max<<"\n";
           fileout<<"Next\n";
            for (int i=0;i<N_max;i++){
                fileout<<"O"<<"\t"<<pos[0][i][0]<<"\t"<<pos[0][i][1]<<"\t"<<0.1<<"\n";
            }
	}

        for (int loopi=0;loopi<N_type;loopi++){
            for (int loopj=0;loopj<Ntype[loopi];loopj++){
                fd_term = sqrt( 2 * epsilon*dt);
                noise_0 = fd_term*gsl_ran_gaussian(r,1);
                noise_1 = fd_term * gsl_ran_gaussian(r,1);
                deltheta1=dt*gamma_i* (ftheta[loopi][loopj])+ noise_0;
		del1=Pe*dt*cos(theta[loopi][loopj]);
		del2=Pe*dt*sin(theta[loopi][loopj]);
                if (fabs(del1)>1| fabs(del2)>1|| pos[loopi][loopj][0]>Lx || pos[loopi][loopj][0]<0||pos[loopi][loopj][1]>Ly||pos[loopi][loopj][1]<0){
                    cout<<del1<<"\t"<<noise_0<<"\t"<<dt * fpos[loopi][loopj][0]<<"Kill program\t"<<pos[loopi][loopj][0]<<"\t"<<pos[loopi][loopj][1]<<"\n";
                    cout<<del2<<"\t"<<noise_1<<"\t"<<dt * fpos[loopi][loopj][1]<<"Kill program\t"<<pos[loopi][loopj][0]<<"\t"<<pos[loopi][loopj][1]<<"\n";
                    cout<<"Time"<<"\t"<<time<<"\t"<<"Count\t"<<count<<"\n";
                    cout<<"Atom"<<"\t"<<loopj<<"\n\n";
                    cout.flush();
                    //fileout<<N_max+4<<"\n";
                    //fileout<<"Next\n";
                    //for (int i=0;i<N_max;i++){
                    //    fileout<<"H"<<"\t"<<pos[1][i][0]<<"\t"<<pos[1][i][1]<<"\t"<<0.1<<"\n";
                    // }
                    //fileout<<"O"<<"\t"<<pos[0][0][0]<<"\t"<<pos[0][0][1]<<"\t"<<0.1<<"\n";
                    //fileout<<"O"<<"\t"<<pos[0][1][0]<<"\t"<<pos[0][1][1]<<"\t"<<0.1<<"\n";
                    //fileout<<"O"<<"\t"<<pos[0][2][0]<<"\t"<<pos[0][2][1]<<"\t"<<0.1<<"\n";
                    //fileout<<"O"<<"\t"<<pos[0][3][0]<<"\t"<<pos[0][3][1]<<"\t"<<0.1<<"\n";
                    //fileout.flush();

                    //Check for big overlaps during equilibration and handle them.
                    if (del1>1)
                        del1=1;
                    if (del1<-1)
                        del1=-1;
                    if (del2>1)
                        del2=1;
                    if (del2<-1)
                        del2=-1;
                    
		    theta[loopi][loopj]+=deltheta1;
                    pos[loopi][loopj][0] += del1;
                    pos[loopi][loopj][1] += del2;
                }
                else{
                    pos[loopi][loopj][0] += del1;
                    pos[loopi][loopj][1] += del2;
		    theta[loopi][loopj]+=deltheta1;
                }
                if (pos[loopi][loopj][0]>Lx)
                    pos[loopi][loopj][0]=pos[loopi][loopj][0]-Lx;
                if (pos[loopi][loopj][0]<0)
                    pos[loopi][loopj][0]+=Lx;
                if (pos[loopi][loopj][1]>Ly)
                    pos[loopi][loopj][1]=pos[loopi][loopj][1]-Ly;
                if (pos[loopi][loopj][1]<0)
                    pos[loopi][loopj][1]+=Ly;
            }
        }
    }

}

double Langevin_dynamics::propogate_dynamics(double dt,double Pe, double t_c, double tau){
    double fd_term,noise_0,noise_1;
    double del1,del2,deltheta1;
    gamma_i=2;
    int typebin=0;
    //Compute forces, evolve dynamics, return y factor
    //cout<<"Reached here\n";
    cout.flush();
    compute_forces_theta();
    for (int loopi=0;loopi<N_type;loopi++){
        for (int loopj=0;loopj<Ntype[loopi];loopj++){
        fd_term = sqrt( epsilon*2 * dt);
        noise_0 = fd_term*gsl_ran_gaussian(r,1);
        noise_1 = fd_term * gsl_ran_gaussian(r,1);
        deltheta1=dt *gamma_i*(ftheta[loopi][loopj])+ noise_0;
        del1=Pe*dt*cos(theta[loopi][loopj]);
        del2=Pe*dt*sin(theta[loopi][loopj]);    
	if (fabs(del1)>0.1*Lx|| fabs(del2)>0.1*Lx){
                cout<<del1<<"\t"<<del2<<"Kill program\n";
                //assert(0);//To check for wildly inappropriate moves.
                if (del1>1)
                    del1=1;
                if (del1<-1)
                    del1=-1;
                if (del2>1)
                    del2=1;
                if (del2<-1)
                    del2=-1;
           }
            theta[loopi][loopj]+=deltheta1;   
            pos[loopi][loopj][0] += del1;
            pos[loopi][loopj][1] += del2;
            if (pos[loopi][loopj][0]>Lx)
                pos[loopi][loopj][0]=pos[loopi][loopj][0]-Lx;
            if (pos[loopi][loopj][0]<0)
                pos[loopi][loopj][0]+=Lx;
            if (pos[loopi][loopj][1]>Ly)
                pos[loopi][loopj][1]=pos[loopi][loopj][1]-Ly;
            if (pos[loopi][loopj][1]<0)
                pos[loopi][loopj][1]+=Ly;
	    if (theta[loopi][loopj]>2*M_PI)
		theta[loopi][loopj]=theta[loopi][loopj]-2*M_PI;
            if (theta[loopi][loopj]<0)
                theta[loopi][loopj]=theta[loopi][loopj]+2*M_PI;
        }
	
    }
    //cout<<"Leaving loop\n";
    double y1=computey(Pe,t_c,tau);
    double y2=exp(-dt*S*y1);
    //cout<<y<<"\n";
    int yc;
    if (fabs(y1*S*dt)>3){
        cout<<"Panic in system"<<"\t"<<y1<<"\n";
        //cout.flush();
        if (y1*S*dt<0)
            return 20.0;
        if (y1*S*dt>0)
            return 0;
        //Absurd value of y are not returned.
    }
    else{
        return y2;//yc is used to compute cloning probabilities.
    }
    return y2;
}

void Langevin_dynamics::compute_forces_theta(){
    double ftheta_t=0;
    double x12=0;
    double y12=0;
    double r12=0;    
    fpos=zerovec;
    fpostype=zerovec;
    upostype=zerovec;
    ftheta=zerovectheta;
    double typebin=1.0;
    for (int loopi=0;loopi<N_type;loopi++){
        for(int loopk=0;loopk<Ntype[loopi];loopk++){
            for(int loopj=0;loopj<N_type;loopj++){
                for(int loopl=0;loopl<Ntype[loopj];loopl++){
                    //if (loopi!=loopj || loopi==loopj){
                    //if (loopi!=loopj || loopk!=loopl){
                        x12=pos[loopi][loopk][0]-pos[loopj][loopl][0];
                        y12=pos[loopi][loopk][1]-pos[loopj][loopl][1];
                        if (x12>0.5*Lx)
                            x12=x12-Lx;
                        if (x12<-0.5*Lx)
                            x12=x12+Lx;
                        if (y12>0.5*Ly)
                            y12=y12-Ly;
                        if (y12<-0.5*Ly)
                            y12=y12+Ly;
                        r12=pow(x12*x12+y12*y12,0.5);
			if (r12>RCUT){
		            ftheta_t=0.0;}
			else{
			    ftheta_t=-sin(theta[loopi][loopk]-theta[loopj][loopl])/(M_PI*RCUT*RCUT);
			}
                        ftheta[loopi][loopk]+=ftheta_t;
		    //}
		   //}
		  }
		}
	       }
	     }
}



void Langevin_dynamics::compute_forces(){
    double f12x=0;
    double f12y=0;
    double x12=0;
    double y12=0;
    double r12=0;
    double f122x=0;
    double f122y=0;
    fpos=zerovec;
    fpostype=zerovec;
    upostype=zerovec;
    double typebin=1.0;
    for (int loopi=0;loopi<N_type;loopi++){
        for(int loopk=0;loopk<Ntype[loopi];loopk++){
            for(int loopj=0;loopj<N_type;loopj++){
                for(int loopl=0;loopl<Ntype[loopj];loopl++){
                    if (loopi!=loopj || loopi==loopj){
                    f12x=0;
                    f12y=0;
                    f122x=0;
                    f122y=0;
                    if (loopi!=loopj || loopk!=loopl){
                        x12=pos[loopi][loopk][0]-pos[loopj][loopl][0];
                        y12=pos[loopi][loopk][1]-pos[loopj][loopl][1];
                        if (x12>0.5*Lx)
                            x12=x12-Lx;
                        if (x12<-0.5*Lx)
                            x12=x12+Lx;
                        if (y12>0.5*Ly)
                            y12=y12-Ly;
                        if (y12<-0.5*Ly)
                            y12=y12+Ly;
                        r12=pow(x12*x12+y12*y12,0.5);
			if (loopi!=loopj)
				typebin=k12*10;
			else 
				typebin=10.0;
            f12x=4*(12*pow(r12,-14)-6*pow(r12,-8))*x12*typebin;
            f12y=4*(12*pow(r12,-14)-6*pow(r12,-8))*y12*typebin;
            f122x=-4*(12*pow(r12,-13)-6*pow(r12,-7))*pow(r12,-1.0)+4*(168*pow(r12,-16)-48*pow(r12,-10.0))*(x12)*x12;
            f122y=-4*(12*pow(r12,-13)-6*pow(r12,-7))*pow(r12,-1.0)+4*(168*pow(r12,-16)-48*pow(r12,-10.0))*(y12)*y12;
			if (soft==1){
			  f12x=-typebin*2*ksoft*exp(-1*pow(r12-asoft,-2.0))*pow(r12-asoft,-3.0)*x12/r12;
			  f12y=-typebin*2*ksoft*exp(-1*pow(r12-asoft,-2.0))*pow(r12-asoft,-3.0)*y12/r12;
			  //f12x=ksoft*(exp(-r12/asoft))*x12/(r12*asoft);
			  //f12y=ksoft*(exp(-r12/asoft))*y12/(r12*asoft);
			  f122x=typebin*ksoft*exp(-1*pow(r12-asoft,-2.0))*(2*pow(r12,8.0)+pow(asoft,3.0)*(2*x12*x12*pow(r12,3.0)-2*pow(r12,5.0))+x12*x12*(4*pow(r12,4.0)-8*pow(r12,6.0))+pow(asoft,2.0)*(-12*x12*x12*pow(r12,4.0)+6*pow(r12,6.0))+asoft*(18*x12*x12*pow(r12,5.0)-6*pow(r12,7.0)))/(pow(r12,6.0)*pow(-asoft+r12,6.0));
              f122y=typebin*ksoft*exp(-1*pow(r12-asoft,-2.0))*(2*pow(r12,8.0)+pow(asoft,3.0)*(2*y12*y12*pow(r12,3.0)-2*pow(r12,5.0))+y12*y12*(4*pow(r12,4.0)-8*pow(r12,6.0))+pow(asoft,2.0)*(-12*y12*y12*pow(r12,4.0)+6*pow(r12,6.0))+asoft*(18*y12*y12*pow(r12,5.0)-6*pow(r12,7.0)))/(pow(r12,6.0)*pow(-asoft+r12,6.0));
			  //f122x=ksoft*(exp(-r12/asoft)*(asoft*x12*x12*pow(r12,3.0)+x12*x12*pow(r12,4.0)-asoft*pow(r12,5.0))*pow(r12,-6.0)*pow(asoft,-2.0));
			  //f122y=ksoft*(exp(-r12/asoft)*(asoft*y12*y12*pow(r12,3.0)+y12*y12*pow(r12,4.0)-asoft*pow(r12,5.0))*pow(r12,-6.0)*pow(asoft,-2.0));

			}
                        if (r12!=0 && r12<pow(2,1.0/6.0) && soft==0 ){
                         fpos[loopi][loopk][0]+=f12x;
                         fpos[loopi][loopk][1]+=f12y;
                            //if(r12<0.75){
                                //cout<<"Panic\t"<<pos[loopi][loopk][0]<<"\t"<<pos[loopj][loopl][0]<<"\t"<<pos[loopi][loopk][1]<<"\t"<<pos[loopj][loopl][1]<<"\t"<<f12x<<"\t"<<f12y<<"\n";
                                //Just to check for big overlaps. Might happen during equilibration.
                            //}
                        }
			if (r12!=0 && r12<asoft && soft==1 ){
                          fpos[loopi][loopk][0]+=f12x;
			  fpos[loopi][loopk][1]+=f12y;
			}
                        if (loopi!=loopj){
                            if (r12!=0 && r12<pow(2,1.0/6.0) && soft==0){
                             fpostype[loopi][loopk][0]+=f12x;
                             fpostype[loopi][loopk][1]+=f12y;
                             upostype[loopi][loopk][0]+=f122x;
                             upostype[loopi][loopk][1]+=f122y;
                            }
                            if (r12!=0 && r12<asoft && soft==1){
			     fpostype[loopi][loopk][0]+=f12x;
			     fpostype[loopi][loopk][1]+=f12y;
			     upostype[loopi][loopk][0]+=f122x;
		             upostype[loopi][loopk][1]+=f122y;
			    }
                        }
                    }
                    }
                }
            }
        }
    }
}

double Langevin_dynamics::computey(double Pe, double t_c, double tau){
    double x12=0;
    double y12=0;
    double y[2];
    y[0]=0;
    y[1]=0;
    double xy[2];
    xy[0]=0;
    xy[1]=0;
    int typebin=0;
    double ystore=0;
    double ftheta_t1;
    double ftheta_t2;
    double r12=0;
    double gamma_i=2;
    //cout<<"Entering compute loop\n";
    for (int loopi=0;loopi<N_type;loopi++){
        for(int loopk=0;loopk<Ntype[loopi];loopk++){
            for(int loopj=0;loopj<N_type;loopj++){
                for(int loopl=0;loopl<Ntype[loopj];loopl++){
                    //if (loopi!=loopj || loopi==loopj){
                    //if (loopi!=loopj || loopk!=loopl){
			x12=pos[loopi][loopk][0]-pos[loopj][loopl][0];
                        y12=pos[loopi][loopk][1]-pos[loopj][loopl][1];
                        if (x12>0.5*Lx)
                            x12=x12-Lx;
                        if (x12<-0.5*Lx)
                            x12=x12+Lx;
                        if (y12>0.5*Ly)
                            y12=y12-Ly;
                        if (y12<-0.5*Ly)
                            y12=y12+Ly;
                        r12=pow(x12*x12+y12*y12,0.5);
                        if (r12>RCUT){
                            ftheta_t1=0.0;
			    ftheta_t2=0.0;
			}
                        else{
                            ftheta_t1=-gamma_i*sin(theta[loopi][loopk]-theta[loopj][loopl])/(M_PI*RCUT*RCUT);
			    ftheta_t2=-gamma_i*cos(theta[loopi][loopk]-theta[loopj][loopl])/(M_PI*RCUT*RCUT);
                        }
	                ystore+=-0.5*epsilon*(ftheta_t2+ftheta_t1*ftheta_t1);
		//}
	    //}
	   }
	 }
	}	
    }
    return ystore;//yfactor computed for cloning.
}



