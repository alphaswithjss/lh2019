//<<<<<<<<<< Output format: Observable Cumulative-Resummation Cummulative-LO Cummulative-delta-NLO
//<<<<<<<<<< Code can be run by: 
//<<<<<<<<<< ./SD-Thrust Observable alpha alphas gluons muR/(pT R) x_L fixedorder zcut beta ppar
//<<<<<<<<<< The default setup is:
//<<<<<<<<<< ./SD-Thrust 0.1 2. 0.118 0 1 1 2 0.1 0 1 
//<<<<<<<<<< Here the variables are given by:
//<<<<<<<<<< Observable: The value of the observable at which the cumulative distribution is computed
//<<<<<<<<<< alpha: power of the angle used in the observable
//<<<<<<<<<< alphas: the stong coupling constant
//<<<<<<<<<< gluons: 1 for a gluon 0 for a quark
//<<<<<<<<<< muR: the scale used in the strong coupling constant
//<<<<<<<<<< x_L: the factor included in the resummed logarithm
//<<<<<<<<<< fixedorder: 1 for LO, 2 for NLO
//<<<<<<<<<< zcut and beta: the Soft drop parameters choose zcut <= 0. for ungroomed
//<<<<<<<<<< ppar: End point correction parameter
//<<<<<<<<<< 

#include <cstring>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <algorithm>

#include "Definitions.h"
#include "Integrand.h"
#include "HypGeo.h"


const double nl = 5.0;
const double Tf = 1.0/2.0;
const double CA = 3.0;
const double CF = 4.0/3.0;
const double pi = M_PI;
const double GammaEuler = 0.5772156649015329;
const double Zeta3 = 1.202056903159594;
const double Zeta2 = pi*pi/6.0;

int f_gluon = 0;
int exp_coef = 0;;
double alpha_Obs = 2.;

double Grid[10000];


int LLorder = 2;
int fixedorder = 1; //(LLorder-1)

double ppar = 1.;

int Delta_PS_flag = 0.;

double  LO_incl = 2.;
double NLO_incl = 5.636921636;


int n_bins_NP = 50;

// double NLO_incl = 0.0;

// char* wd;

bool end_point = true;

bool F_HypGeo = true;

using namespace std;

const int points_grid = 9;

string DoubleToStr(double n)
{
    stringstream result;
    double tiny = 1E-11*n;
    result/* << fixed*/ << setprecision(8) << n+tiny;
//     cout << n << "\t" << setprecision(8) << n << "\t" << result.str() << endl;
    return result.str();
}

void readingrids(void){
        ifstream ifilegg("Grid-range.dat",ios::in);
	if (!ifilegg.is_open()) {
	  cout << "There was a problem opening the input file!" << endl;
	  exit(1);
	}
	int n = 0.;
	while(ifilegg >> Grid[n]) n++;	
        
	return;
}

void test_expansion(double bin_min, double bin_max, double n_bins, bool b_log, double zc, double beta, double alphas = 0.117, int order=-1){
    fixedorder = 2;
//     cout << "#taum taump exp1 exp2 res" << endl;
    for(double i=0;i<n_bins;i++){
        double result_exp1_cum_p,result_exp1_cum_m;
        double result_exp2_cum_p,result_exp2_cum_m;
        double result_exp3_cum_p,result_exp3_cum_m;
        double result_res_cum_p,result_res_cum_m;
        double tau_valp = 0.;
        double tau_valm = 0.;
        double Ltau_valp = 0.;
        double Ltau_valm = 0.;
        if(b_log){
            tau_valp = bin_min*pow(bin_max/bin_min, (i+1.)/n_bins);
            tau_valm = bin_min*pow(bin_max/bin_min, i/n_bins);
//             Ltau_valp = (bin_max-bin_min)*(i+1.)/n_bins+bin_min;
//             Ltau_valm = (bin_max-bin_min)*i/n_bins+bin_min;
//             tau_valp = exp(-Ltau_valm);
//             tau_valm = exp(-Ltau_valp);
        }
        else{
            tau_valp = (bin_max-bin_min)*(i+1.)/n_bins+bin_min;
            tau_valm = (bin_max-bin_min)*i/n_bins+bin_min;
        }
        para temp_p = {alphas,1.,1.,tau_valp};
        paraSD tempSD_p = {alphas,1.,1.,tau_valp,zc,beta};
        para temp_m = {alphas,1.,1.,tau_valm};
        paraSD tempSD_m = {alphas,1.,1.,tau_valm,zc,beta};
//         cout << zc << "\t" << beta << endl;
        if(zc <= 0.){
//             cout << "zc = 0" << endl;
            exp_coef = 1;
            result_exp1_cum_p = Expansion_cumu(&temp_p);
            result_exp1_cum_m = Expansion_cumu(&temp_m);
            exp_coef = 2;
            result_exp2_cum_p = Expansion_cumu(&temp_p);
            result_exp2_cum_m = Expansion_cumu(&temp_m);
            
            exp_coef = 0;
            result_res_cum_p = dQCD_cumu(&tempSD_p);
            result_res_cum_m = dQCD_cumu(&tempSD_m);
        }
        else{
            exp_coef = 1;
            result_exp1_cum_p = ExpansionSD_cumu(&tempSD_p);
            result_exp1_cum_m = ExpansionSD_cumu(&tempSD_m);
            exp_coef = 2;
            result_exp2_cum_p = ExpansionSD_cumu(&tempSD_p);
            result_exp2_cum_m = ExpansionSD_cumu(&tempSD_m);
            exp_coef = 0;
            result_res_cum_p = dQCDSD_cumu(&tempSD_p);
            result_res_cum_m = dQCDSD_cumu(&tempSD_m);
        }
        if(tau_valm==0.){
            result_exp1_cum_m = 0.;
            result_exp2_cum_m = 0.;
            result_res_cum_m = 0.;
        }
        if(result_res_cum_m!=result_res_cum_m) result_res_cum_m = 0.;
        if(result_res_cum_p!=result_res_cum_p) result_res_cum_p = 0.;
        
        if(i==0 && tau_valm!=0.) cout << 0. << "\t" << tau_valm << "\t" << result_exp1_cum_m << "\t" << result_exp2_cum_m << "\t" << result_res_cum_m << endl;
//                 cout << tau_valp << "\t" << result_exp1_cum_p  << "\t" << result_exp2_cum_p << "\t" << result_res_cum_p << endl;
        cout << tau_valm << "\t" << tau_valp << "\t" << (result_exp1_cum_p-result_exp1_cum_m)  << "\t" << (result_exp2_cum_p-result_exp2_cum_m) << "\t" << (result_res_cum_p-result_res_cum_m) << endl;
        
    }
//     exit(0);
    return;
}
       
void test_transition(double zc, double beta, bool f_end = true, double alphas = 0.118, double epsi=1E-4){
    if(F_HypGeo) readingrids();
    
    double result_exp1_cum[5], result_exp2_cum[5], result_exp3_cum[5], result_res_cum[5];
    double dresult_exp1_cum[5], dresult_exp2_cum[5], dresult_exp3_cum[5], dresult_res_cum[5];
    double taum = 1./4.;
    if(fixedorder == 1) taum = 1./4.;
    if(fixedorder == 2) taum = 1./3.;
    
    double tau_trans = taum*zc/(zc+taum-taum*zc);
    if(!end_point) tau_trans = zc;
    
    if(f_end){
        tau_trans = taum;
    }
    
    for(double i=0;i<5;i++){
        
        double tau_val = tau_trans-2.*epsi+epsi*i;
        cout << tau_val << endl;
        
        
        if(zc <= 0.){
            paraSD temp = {alphas,1.,1.,tau_val};
            
            exp_coef = 1;
            result_exp1_cum[(int) i] = Expansion_cumu(&temp);
            exp_coef = 2;
            result_exp2_cum[(int) i] = Expansion_cumu(&temp);
            exp_coef = 0;
            result_exp3_cum[(int) i] = Expansion_cumu(&temp)-1.-result_exp1_cum[(int) i]-result_exp2_cum[(int) i];
            
            result_res_cum[(int) i] = dQCD_cumu(&temp);
        }
        else{
            paraSD tempSD = {alphas,1.,1.,tau_val,zc,beta};
        
            exp_coef = 1;
            result_exp1_cum[(int) i] = ExpansionSD_cumu(&tempSD);
            exp_coef = 2;
            result_exp2_cum[(int) i] = ExpansionSD_cumu(&tempSD);
            exp_coef = 0;
            result_exp3_cum[(int) i] = ExpansionSD_cumu(&tempSD)-1.-result_exp1_cum[(int) i]-result_exp2_cum[(int) i];
            
            result_res_cum[(int) i] = dQCDSD_cumu(&tempSD);
        }
        if(i>0){
            dresult_exp1_cum[(int) i] = (result_exp1_cum[(int) i]-result_exp1_cum[(int) i-1])/epsi;
            dresult_exp2_cum[(int) i] = (result_exp2_cum[(int) i]-result_exp2_cum[(int) i-1])/epsi;
            dresult_exp3_cum[(int) i] = (result_exp3_cum[(int) i]-result_exp3_cum[(int) i-1])/epsi;
            dresult_res_cum[(int) i] = (result_res_cum[(int) i]-result_res_cum[(int) i-1])/epsi;
        }
    }
    cout << endl;
    for(int i=0;i<5;i++){
        if(i>0){
            cout << result_exp1_cum[i] << "\t" << dresult_exp1_cum[i] << endl;
        }
        else{
            cout << result_exp1_cum[i] << "\t" << 0. << endl;
        }
    }
    cout << endl;
    for(int i=0;i<5;i++){
        if(i>0){
            cout << result_exp2_cum[i] << "\t" << dresult_exp2_cum[i] << endl;
        }
        else{
            cout << result_exp2_cum[i] << "\t" << 0. << endl;
        }
    }
    cout << endl;
    for(int i=0;i<5;i++){
        if(i>0){
            cout << result_exp3_cum[i] << "\t" << dresult_exp3_cum[i] << endl;
        }
        else{
            cout << result_exp3_cum[i] << "\t" << 0. << endl;
        }
    }
    cout << endl;
    for(int i=0;i<5;i++){
        if(i>0){
            cout << result_res_cum[i] << "\t" << dresult_res_cum[i] << endl;
        }
        else{
            cout << result_res_cum[i] << "\t" << 0. << endl;
        }
    }
        
//     exit(0);
    return;
}

int main(int argc,char** argv){
  double chiR,chiQ;
  double Observable;
  double zcut;
  double beta;
  double alphas;
  
  
  if(argc > 1) sscanf(argv[1],"%lf",&Observable);
  else Observable = 0.1;
  if(argc > 2) sscanf(argv[2],"%lf",&alpha_Obs);
  else alpha_Obs = 2.;
  if(argc > 3) sscanf(argv[3],"%lf",&alphas);
  else alphas = 0.118;
  if(argc > 4) sscanf(argv[4],"%d",&f_gluon);
  else f_gluon = 0;
  if(argc > 5) sscanf(argv[5],"%lf",&chiR);
  else chiR = 1.;
  if(argc > 6) sscanf(argv[6],"%lf",&chiQ);
  else chiQ = 1.;
  if(argc > 7) sscanf(argv[7],"%d",&fixedorder);
  else fixedorder = 2;
  if(argc > 8) sscanf(argv[8],"%lf",&zcut);
  else zcut = 0.1;
  if(argc > 9) sscanf(argv[9],"%lf",&beta);
  else beta = 0.0;
  if(argc > 10) sscanf(argv[10],"%lf",&ppar);
  else ppar = 1.;
          
  double rscaleR, rscaleQ;
  
  rscaleQ = 1.0/chiQ/chiQ;
  rscaleR = 1.0/chiR/chiR;
  
  double b0 = (11.*CA-2.*nl)/12./pi;
  double b1 = (34./3.*CA*CA-4.*Tf*nl*(CF+5./3.*CA))/16./pi/pi;
  
  double Resummation;
  double Expansion1;
  double Expansion2;
  
//   test_transition(zcut, beta, true, alphas, 1E-8);
  
  if(zcut <= 0.){
      
      para temp = {alphas,rscaleR,rscaleQ,Observable};
      
      exp_coef = 1;
      Expansion1 = Expansion_cumu(&temp);
      exp_coef = 2;
      Expansion2 = Expansion_cumu(&temp);
      
      Resummation = dQCDSD_cumu(&temp);
  }
  else{
      
      readingrids();
      paraSD tempSD = {alphas,rscaleR,rscaleQ,Observable,zcut,beta};
      
      exp_coef = 1;
      Expansion1 = Expansion_cumu(&tempSD);
      exp_coef = 2;
      Expansion2 = Expansion_cumu(&tempSD);
      
      Resummation = dQCDSD_cumu(&tempSD);
  }
  
  cout << Observable << "  " <<  Resummation << "  " << Expansion1 << "  " << Expansion2 << endl;
        
  return 0;
}
