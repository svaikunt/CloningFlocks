
/*
CLONING simulation. */

#define PI 3.14159
//#define CLUSTER
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#ifdef CLUSTER
#include "/home/svaikunt/local/lib/include/fftw3.h"
#endif 
#ifndef CLUSTER
#include <fftw3.h>
#endif
#include <assert.h>
#include "Langevin_dynamics.h"

#define ARGS 8

using namespace std;
int main( int argc,char *argv[]){

    int loopi,loopj;
    double growthcgf=1;;
    if (argc==1){
	  cout<<"Cloning simulation\n";
	  cout<<"Usage N_max(number of particles of type1, there are 4 of type 0) t_analysis randomseed timestep  S(biasing function) Pe tau soft ksoft_solute_solvent N_Snapshots\n";
	  exit(1);
	}
	int N_max=atoi(argv[1]);
	double t_analysis=atof(argv[2]);
	long int randomseed=atoi(argv[3]);
	double dt=atof(argv[4]);
	double S=atof(argv[5]);
	double Pe=atof(argv[6]);
	double tau=atof(argv[7]);
	int soft=atoi(argv[8]);
	double ksoft=atof(argv[9]);
	int N_snapshots=atoi(argv[10]);//number of clones
        double Lx=4.0;
        double Ly=4.0;
    //Initiating gsl
    gsl_rng *r=gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(r,randomseed);
    

    
    //Iniatiate N_snapshots or clones of mylattice.
    Langevin_dynamics mylattice[N_snapshots];
    double y[N_snapshots];
    double yc[N_snapshots];
    double yc2[N_snapshots];
    double sumy=0;
	
    vector< vector < vector < vector<double> > > > tpos;
    
    
    int N_type=1;// number of particle types
    int Ntype[0];
    Ntype[0]=N_max;//Number of particles of type 1 input by user.
    for (int i=0;i<N_snapshots;i++){
        vector< vector< vector<double> > > pos;
        for (int loopi=0;loopi<N_type;loopi++){
            vector< vector<double> > pos2;
            for (int loopj=0;loopj<Ntype[loopi];loopj++){
                vector<double> pos1;
                for(int loopk=0;loopk<2;loopk++){
                    pos1.push_back(0);
                }
                pos2.push_back(pos1);
            }
            pos.push_back(pos2);
        }
        tpos.push_back(pos);
    }
    
    Langevin_dynamics masterlattice;
    masterlattice.initialize(0,Ntype[0],randomseed,S,soft,ksoft);
    cout<<"Clone>>>>>>>>>>>>>>>Equilibrating master\n";
    masterlattice.equilibrate(Pe,tau);

     //////
    char outputfile[100];
    sprintf(outputfile,"SnapshotsN%d.N2%d.S%.3f.Lx%.2f.Clone%d.Soft%d.E10.ktracer%.2f.lammpstrj",N_max,Ntype[0],S,Lx,N_snapshots,soft,ksoft);
    ofstream fileout;
    fileout.open(outputfile);
    sprintf(outputfile,"NCloneSnapshotsN%d.N2%d.S%.3f.Lx%.2f.Clone%d.Soft%d.E10.ktracer%.2f.lammpstrj",N_max,Ntype[0],S,Lx,N_snapshots,soft,ksoft);
    ofstream fileoutclone;
    fileoutclone.open(outputfile);
    //////

    for (int i=0;i<N_snapshots;i++){//loop to generate i clones.
        cout<<"Clone>>>>>>>>>>>>>>>\t"<<i<<"\n";
        for (int j=0;j<10.0/(dt*N_snapshots);j++){//the master lattice is run for 10/(dt*N_snapshots) steps.
          double dump=masterlattice.propogate_dynamics(dt,Pe,0,tau);
	  if (j%5==0){
           fileout<<"ITEM: TIMESTEP\n";
	   fileout<<j<<"\n";
	   fileout<<"ITEM: NUMBER OF ATOMS\n";
	   fileout<<N_max<<"\n";
	   fileout<<"ITEM: BOX BOUNDS\n";
	   fileout<<"0\t"<<Lx<<"\n";
	   fileout<<"0\t"<<Ly<<"\n";
	   fileout<<"0.1  0.1\n";
           fileout<<"ITEM: ATOMS id type x y z vx vy vz\n";
           for (int i=0;i<Ntype[0];i++){
            fileout<<i<<"\t"<<"0"<<"\t"<<masterlattice.pos[0][i][0]<<"\t"<<masterlattice.pos[0][i][1]<<"\t"<<0.1<<"\t"<<Pe*cos(masterlattice.theta[0][i])<<"\t"<<Pe*sin(masterlattice.theta[0][i])<<"\t"<<"0"<<"\n";
           }//Output is written out in Snapshot.....XYZ
	  }
	}
        cout<<"Clone>>>>>>>>>>>>>>>\t Initializing"<<i<<"\n";
	cout.flush();
	mylattice[i].initialize(0,Ntype[0],randomseed+i,S,soft,ksoft);//i^th close is created
        //mylattice[i].equilibrate(Pe,tau);
        mylattice[i].pos=masterlattice.pos;//Positions 
	mylattice[i].theta=masterlattice.theta;//and thetas of master lattice are stored in the i^th clone.  
        tpos[i]=mylattice[i].pos;
    }
    

    ///
    char outputfile1[100];
    sprintf(outputfile1,"g_rN%d.N2%d.S%.3f.Lx%.2f.Clone%d.Soft%d.ktracer%.2f.E10.XYZ",N_max,Ntype[0],S,Lx,N_snapshots,soft,ksoft);
    ofstream fileoutg;
    fileoutg.open(outputfile1);//File to potentially write out information related to correlations 
    ///
    
    
    double xdis,ydis,r2;
    int countdt=0;
    int dumpdt=0;//Main cloning loop starts here.
    for (double teetotaler=0;teetotaler<t_analysis;teetotaler+=dt){
        countdt+=1;
        sumy=0;
        double sumyc=0;
        if (countdt>dumpdt){
         for (int loopj=0;loopj<N_snapshots;loopj++){
            //j=pick randomly from the clones
	    y[loopj]=0;
	    for (int localtime=0;localtime<20;localtime++){
             y[loopj]+=mylattice[loopj].propogate_dynamics(dt,Pe,teetotaler,tau);//propogate clone.Store value of the trajectory dependent variable \cal E in y[loopi]
	    }
            tpos[loopj]=mylattice[loopj].pos;//storing clone positions for switching.
            y[loopj]=exp(y[loopj]);//The output from propogate_dynamcis is -S dt \dot{E}. This is added together and then exponentiated. 
	    sumy+=y[loopj];
         }
         for (int loopj=0;loopj<N_snapshots;loopj++){
            yc[loopj]=(int)(floor(y[loopj]/sumy*N_snapshots+gsl_rng_uniform(r)));//use procedure in PRE paper to generate random numbers according to \cal E 
            //cout<<"Testing yc\t"<<yc[loopj]<<"\n";
            sumyc+=yc[loopj];
            yc2[loopj]=sumyc;
         }
         if(sumyc==N_snapshots){
            int clonecount=0;
            for (int loopj=0;loopj<N_snapshots;loopj++){
                for (int loopk=0;loopk<yc[loopj];loopk++){
                    mylattice[loopk+clonecount].pos=tpos[loopj];//Cloning performed !!!
                }
                clonecount+=yc[loopj];
                //cout<<loopj<<"\t"<<clonecount<<"\n";
                //cout.flush();
            }
         }
         if (sumyc!=N_snapshots){
            for (int loopj=0;loopj<N_snapshots;loopj++){
                double randtemp=gsl_rng_uniform(r)*sumyc;
                int iterate=-1;
                do{
                    iterate+=1;
                }while(yc2[iterate]<=randtemp);
                if (iterate>N_snapshots-1){
                    cout<<"Iterate\t"<<iterate<<"\t"<<randtemp<<"\t"<<yc[N_snapshots-1]<<"\n";
                    iterate=N_snapshots-1;
                }
                mylattice[loopj].pos=tpos[iterate];//Cloning performed !! 
                if(teetotaler>0.0*t_analysis)
                    cout<<"Clone\t"<<loopj<<"\t is now clone \t"<<iterate<<"\t"<<yc2[iterate]<<"\t"<<randtemp<<"\n";
            }
         }
         double ratio=((double)(sumy)*pow(N_snapshots,-1.0));
         growthcgf=growthcgf*ratio;//This variable stores a running value for the cumulant generating function. Not crucial for our current project since we are looking for something qualitative.
        //Overall scheme: compute y for clone. Kill (copy some other clone into this) or select clones randomly to copy into
        //keep track of growth function.
       }
       if (countdt%10==0){
        fileoutclone<<"ITEM: TIMESTEP\n";
        fileoutclone<<countdt<<"\n";
        fileoutclone<<"ITEM: NUMBER OF ATOMS\n";
        fileoutclone<<N_max<<"\n";
        fileoutclone<<"ITEM: BOX BOUNDS\n";
        fileoutclone<<"0\t"<<Lx<<"\n";
        fileoutclone<<"0\t"<<Ly<<"\n";
        fileoutclone<<"0.1  0.1\n";
        fileoutclone<<"ITEM: ATOMS id type x y z vx vy vz\n";
        for (int i=0;i<Ntype[0];i++){
          fileoutclone<<i<<"\t"<<"0"<<"\t"<<mylattice[0].pos[0][i][0]<<"\t"<<mylattice[0].pos[0][i][1]<<"\t"<<0.1<<"\t"<<Pe*cos(mylattice[0].theta[0][i])<<"\t"<<Pe*sin(mylattice[0].theta[0][i])<<"\t"<<"0"<<"\n";
        }//Writeout the data from clone 0 into a movie CloneSnapshot.......XYZ
       }

     }
     cout<<"CGF:"<<log(growthcgf)/t_analysis<<"\n";
    
    //Compute averages over the cloned lattices.
    double avgy=0; //Average value of dU12/dt will be stored in avgy.
    sprintf(outputfile,"ystatsNmax%d.N2%d.S%.3f.Lx%.2f.Clone%d.Soft%d.E10.ktracer%.2f.XYZ",N_max,Ntype[0],S,Lx,N_snapshots,soft,ksoft);
    ofstream fileoutystats;
    fileoutystats.open(outputfile);
    fileoutystats<<"CGF\t"<<log(growthcgf)/t_analysis<<"\n";
    for (int loopi=0;loopi<N_snapshots;loopi++){
        mylattice[loopi].compute_forces();
        double tempy=mylattice[loopi].computey(Pe,t_analysis,tau);
        fileoutystats<<tempy<<"\n";
        cout<<"Compute y is "<<tempy<<"\t Clone is \t"<<loopi<<"\n";
        avgy+=tempy;
    }
    cout<<"avgenergy:"<<avgy*pow(N_snapshots,-1.0)<<"\n";
    
    
}

