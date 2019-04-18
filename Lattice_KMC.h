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
#include <stdio.h>
#include <fstream>
#include <vector>

using namespace std;


class Lattice_KMC{
private:
	double t_analysis;
	int N_max;
	int count;
	double time;
	double t_prev_analysis;
	int winding;
	double S;
public:
  
  ofstream fileout_traj,fileout_activity;
  
  
  gsl_rng *r;
  long int randomseed;
  char outputfile_traj[100];
  char outputfile_activity[100];
  double Loop_left_right;
  double Loop_left_right_0;
  double Loop_right_left;
  double Loop_right_left_0;
  double S_current;
  double S_old;
  double S_newchunk;
  double S_oldchunk;
  double time_test;
  vector < vector<double> > trajectory;
  vector < vector<double> > trajectory_analysis;
  vector < vector<double> > trajectory_TPSchunk;
  vector < vector<double> > trajectory_analysis_temp;
  vector < vector<double> > trajectory_oldTPSchunk;
  vector < vector<double> > W;
  vector < vector<double> > W0;
  vector <double> time_series;
  vector <double> timeseries_TPSchunk;
  vector <double> timeseries_oldTPSchunk;
  vector <double> time_series_temp;
  vector <double> system;
  vector <double> s_history;
  vector <double> t_history;
  vector <double> s2_history;
  vector <double> system_final;
  vector <double> system_initial;
  vector <double> system_temp;
  
  
  Lattice_KMC();
  ~Lattice_KMC();
  void initialize(int, double, long int,double);
  void initializelattice();
  double Sratematrix(vector <double>&, vector <double>& );
  //double adjacencymatrix(vector <double> , vector <double> );
  double kmc_transition();
  void transition_to(int );
  void KMC_evolve_extract(double);
  void analyze_traj(double );
  void writeout_traj();
  void writeout_activity();
  void propagate(int);
  void generate_trajectory_offshoot(double);
  void generate_trajectory_test(double);
  double traj_entropy(vector< vector<double> >&,double );
};

#endif



