/*
Library of functions for lattice KMC simulations and otherwise. C++ implemntations, making them classes. 
*/
#ifndef SOSGAURD
#define SOSGAURD
#include <cstdlib>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>
#include <fstream>
#include <vector>

using namespace std;


class Langevin_dynamics{
private:
	double t_analysis;
	double Lx;
	double Ly;
	double Ntype[2];
	double N_type;
	int N_max;
	double S;
	double gamma_i;
	double ksoft;
	double asoft;
	int soft;
	double k12;
	double epsilon;
	double RCUT;
public:
  
  ofstream fileout_traj,fileout_activity;
  
  
  gsl_rng *r;
  long int randomseed;
  
  vector < vector <double> > zerovectheta; 
  vector < vector <double> > ftheta;
  vector < vector <double> > theta;
  vector < vector < vector<double> > > pos;
  vector < vector < vector<double> > > fpos;
  vector < vector < vector<double> > > fpostype; 
  vector < vector < vector<double> > > upostype;
  vector < vector < vector<double> > > zerovec;

  
  Langevin_dynamics();
  ~Langevin_dynamics();
  void initialize(int, int, long int,double,int,double);
  double propogate_dynamics(double,double, double, double);
  void compute_forces();
  void compute_forces_theta();
  double computey(double, double, double);
  void equilibrate(double, double);
};

#endif



