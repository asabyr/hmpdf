#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "utils.h"
#include "configs.h"
#include "object.h"
#include "halo_model.h"

int
null_halo_model(hmpdf_obj *d)
{//{{{
    STARTFCT

    d->h->inited_halo = 0;
    d->h->hmf = NULL;
    d->h->bias = NULL;
    d->h->c_interp = NULL;
    d->h->c_accel = NULL;

    ENDFCT
}//}}}

int
reset_halo_model(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\treset_halo_model\n");

    if (d->h->hmf != NULL)
    {
        for (int z_index=0; z_index<d->n->Nz; z_index++)
        {
            if (d->h->hmf[z_index] != NULL) { free(d->h->hmf[z_index]); }
        }
        free(d->h->hmf);
    }
    if (d->h->bias != NULL)
    {
        for (int z_index=0; z_index<d->n->Nz; z_index++)
        {
            if (d->h->bias[z_index] != NULL) { free(d->h->bias[z_index]); }
        }
        free(d->h->bias);
    }
    if (d->h->c_interp != NULL) { gsl_spline_free(d->h->c_interp); }
    if (d->h->c_accel != NULL)
    {
        for (int ii=0; ii<d->Ncores; ii++)
        {
            if (d->h->c_accel[ii] != NULL)
            {
                gsl_interp_accel_free(d->h->c_accel[ii]);
            }
        }
        free(d->h->c_accel);
    }

    ENDFCT
}//}}}



static int
splint(double x_arr[],
       double y_arr[],
       double y2_arr[],
       int n_points,
       double x_value,
       double *y_value)
{
  // spline interpolation function from class_sz: splint in class_sz_tools.c https://github.com/CLASS-SZ/class_sz/blob/master/class-sz/tools/class_sz_tools.c
  // it's from Sec 3.3 of Press+1992 as cited in Tinker+2008: https://arxiv.org/pdf/0803.2706
  int k_low, k_high, k;
  float h,b,a;

  k_low=0;
  k_high = n_points-1;
  
  while (k_high-k_low > 1)
  {
    k = (k_high+k_low) >> 1;//shift operator divides by 2 so k=(k_high+k_low)/2
    if (x_arr[k] > x_value)
      k_high = k;
    else
      k_low = k;
  }

  h = x_arr[k_high] - x_arr[k_low];
  if (h == 0.0) 
    return 0;

  a = (x_arr[k_high] - x_value)/h;
  b = (x_value-x_arr[k_low])/h;

  *y_value =
  a*y_arr[k_low]
  +b*y_arr[k_high]
  +((a*a*a-a)*y2_arr[k_low]
    +(b*b*b-b)*y2_arr[k_high])
  *(h*h)
  /6.0;

  return 1;
}

static int 
DeltaVir_BryanNorman98(hmpdf_obj *d, int z_index, double *out)
{//{{{
    STARTFCT

    double x = d->c->Om[z_index] - 1.0;
    *out = 18.0*M_PI*M_PI + 82.0*x - 39.0*x*x;

    ENDFCT
}//}}}

static int
density_threshold(hmpdf_obj *d, int z_index, hmpdf_mdef_e mdef, double *out)
{//{{{{
    STARTFCT

    double dvir = 0.0; // to avoid maybe-uninitialized
    switch (mdef)
    {
        case hmpdf_mdef_c : *out = 200.0 * d->c->rho_c[z_index];
                            break;
        case hmpdf_mdef_v : SAFEHMPDF(DeltaVir_BryanNorman98(d, z_index, &dvir));
                            *out = dvir * d->c->rho_c[z_index];
                            break;
        case hmpdf_mdef_m : *out = 200.0 * d->c->rho_m[z_index];
                            break;
        default           : *out = 0.0; // to avoid maybe-uninitialized
                            HMPDFERR("Unknown mass definition.");
    }

    ENDFCT
}//}}}

int 
RofM(hmpdf_obj *d, int z_index, int M_index, double *out,
     double mass_resc)
{//{{{
    STARTFCT

    double dt;
    SAFEHMPDF(density_threshold(d, z_index, d->h->HMF_mdef, &dt));
    *out = cbrt(3.0*mass_resc*d->n->Mgrid[M_index] / 4.0 / M_PI / dt);

    ENDFCT
}//}}}

static int
MofR(hmpdf_obj *d, int z_index, double R, hmpdf_mdef_e mdef, double *out)
{//{{{
    STARTFCT

    double dt;
    SAFEHMPDF(density_threshold(d, z_index, mdef, &dt));
    *out = 4.0 * M_PI * dt * gsl_pow_3(R) / 3.0;

    ENDFCT
}//}}}

static inline double
c_Duffy08_1(hmpdf_obj *d, double z, double M, hmpdf_mdef_e mdef, double *conc_params)
// for conc_params, can pass NULL to use the default d->h->Duffy08_params
{//{{{
    // convert to Msun/h
    M *= d->c->h / 2e12;

    if (conc_params == NULL)
        conc_params = d->h->Duffy08_params;

    // take the correct mass definition
    conc_params += (int)(mdef) * 3;

    switch (mdef)
    {
        // this case is special, since it has more parameters
        //     (the other definitions are not used for important things)
        case (hmpdf_mdef_m) :
        {
            // we harden this a bit against weird outputs during MCMC sampling
            static const double min_out=CINTERP_CMIN*1.01, max_out=CINTERP_CMAX*0.99;

            double pref = conc_params[0]
                          * pow(M,     conc_params[1])
                          * pow(1.0+z, conc_params[2]);
            double arg = conc_params[3]
                         * gsl_pow_2(z-conc_params[4])
                         * gsl_pow_2(log(M) - conc_params[5]);
            // avoid underflow and overflow
            if (arg < log(fmin(1.0, pref)) + log((double)FLT_RADIX) * (double)DBL_MIN_EXP + 2.0)
                return min_out;
            else if (arg > -log(fmax(1.0, pref)) + log((double)FLT_RADIX) * (double)DBL_MAX_EXP - 2.0)
                return max_out;
            double out = pref * exp( arg );
            return (out<min_out) ? min_out : (out>max_out) ? max_out : out;
        }
        // use the Duffy08 parameterization
        default :
            return conc_params[0]
                   * pow(M,      conc_params[1])
                   * pow(1.0+z,  conc_params[2]);
    }
}//}}}
static inline double
c_Duffy08(hmpdf_obj *d, int z_index, int M_index,
          double mass_resc, double *conc_params)
{//{{{
    double out = c_Duffy08_1(d, d->n->zgrid[z_index],
                             mass_resc * d->n->Mgrid[M_index],
                             d->h->HMF_mdef, conc_params);

    if (d->h->conc_resc != NULL)
    {
        out *= d->h->conc_resc(d->n->zgrid[z_index],
                               mass_resc * d->n->Mgrid[M_index] * d->c->h,
                               d->h->conc_resc_params);
    }

    return out;
}//}}}

int
NFW_fundamental(hmpdf_obj *d, int z_index, int M_index,
                double mass_resc, double *conc_params,
                double *rhos, double *rs)
// returns rhos via function call and rs via return value
// this function is tested against Colossus --> everything here works
{//{{{
    STARTFCT

    double c = c_Duffy08(d, z_index, M_index, mass_resc, conc_params);
    SAFEHMPDF(RofM(d, z_index, M_index, rs, mass_resc));
    *rs /= c;
    *rhos = mass_resc * d->n->Mgrid[M_index]/4.0/M_PI/gsl_pow_3(*rs)
            / (log1p(c)-c/(1.0+c));

    ENDFCT
}//}}}

static int
create_c_of_y(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\tcreate_c_of_y\n");

    double *logc_grid;
    double *logy_grid;
    SAFEALLOC(logc_grid, malloc(CINTERP_NC * sizeof(double)));
    SAFEALLOC(logy_grid, malloc(CINTERP_NC * sizeof(double)));
    for (int ii=0; ii<CINTERP_NC; ii++)
    {
        logc_grid[ii] = log(CINTERP_CMAX)
                        - (double)(ii)*log(CINTERP_CMAX/CINTERP_CMIN)/(double)(CINTERP_NC-1);
        double _c = exp(logc_grid[ii]);
        logy_grid[ii] = log(3.0) - 3.0*logc_grid[ii] + log(log1p(_c) - _c/(1.0+_c));
    }
    SAFEALLOC(d->h->c_interp, gsl_spline_alloc(gsl_interp_cspline, CINTERP_NC));
    SAFEALLOC(d->h->c_accel,  malloc(d->Ncores * sizeof(gsl_interp_accel *)));
    SETARRNULL(d->h->c_accel, d->Ncores);
    for (int ii=0; ii<d->Ncores; ii++)
    {
        SAFEALLOC(d->h->c_accel[ii], gsl_interp_accel_alloc());
    }
    SAFEGSL(gsl_spline_init(d->h->c_interp, logy_grid, logc_grid, CINTERP_NC));
    free(logc_grid);
    free(logy_grid);

    ENDFCT
}//}}}

static int
c_of_y(hmpdf_obj *d, double y, double *out)
// inverts the function y(c) = 3/c^3 * (log(1+c)-c/(1+c))
// this function is much smoother on loglog scale, so do the interpolation this way
{//{{{
    STARTFCT

    SAFEGSL(gsl_spline_eval_e(d->h->c_interp, log(y),
                              d->h->c_accel[THIS_THREAD], out));
    *out = exp(*out);

    ENDFCT
}//}}}

int
Mconv(hmpdf_obj *d, int z_index, int M_index, hmpdf_mdef_e mdef_out,
      double mass_resc,
      double *M, double *R, double *c)
// returns the converted mass via function call and the new radius and concentration via return value
// this function is tested against Colossus --> everything here works
{//{{{
    STARTFCT

    double rhos, rs;

    // use default concentration model here
    SAFEHMPDF(NFW_fundamental(d, z_index, M_index, mass_resc, NULL, &rhos, &rs));

    double dt;
    SAFEHMPDF(density_threshold(d, z_index, mdef_out, &dt));
    SAFEHMPDF(c_of_y(d, dt/rhos, c));
    *R = rs * *c;
    SAFEHMPDF(MofR(d, z_index, *R, mdef_out, M));

    ENDFCT
}//}}}

static inline double
fnu_Tinker10_primitive(hmpdf_obj *d, int n, double z)
{//{{{
    return d->h->Tinker10_params[n*2]
           *pow(1.0+z, d->h->Tinker10_params[n*2 + 1]);
}//}}}

static double
fnu_Tinker10(hmpdf_obj *d, double nu, double z)
{//{{{
    z = (z<3.0) ? z : 3.0;
    double beta  = fnu_Tinker10_primitive(d, 0, z);
    double phi   = fnu_Tinker10_primitive(d, 1, z);
    double eta   = fnu_Tinker10_primitive(d, 2, z);
    double gamma = fnu_Tinker10_primitive(d, 3, z);
    double alpha = fnu_Tinker10_primitive(d, 4, z);
    return nu * alpha*(1.0+pow(beta*nu, -2.0*phi))
           * pow(nu, 2.0*eta) * exp(-0.5*gamma*gsl_pow_2(nu));
}//}}}

static double
bnu_Tinker10(double nu)
{//{{{
    double y = 2.0 + M_LN2/M_LN10; // y = log_10(200)
    double A = 1.0 + 0.24 * y * exp(-gsl_pow_4(4.0/y));
    double a = 0.44 * y - 0.88;
    double B = 0.183;
    double b = 1.5;
    double C = 0.019 + 0.107 * y + 0.19 * exp(-gsl_pow_4(4.0/y));
    double c = 2.4;
    return 1.0 - A*pow(nu, a)/(pow(nu, a) + pow(1.686, a)) + B*pow(nu, b) + C*pow(nu, c);
}//}}}

static double
get_Omega_m_nonu_at_z(hmpdf_obj *d, double z){
    
    //get_Omega_m_nonu_at_z function in class_sz_tools.c https://github.com/CLASS-SZ/class_sz/blob/master/class-sz/tools/class_sz_tools.c
    
    double Om_0 = d->c->Om_0;
    double Om_0_nonu = d->c->Oc_0 + d->c->Ob_0;
    double Or_0 = d->c->Or_0;
    double Ol_0 = 1. - Om_0 - Or_0;
    double Om_z = Om_0_nonu * pow(1. + z, 3.) / (Om_0 * pow(1. + z, 3.) + Ol_0 + Or_0 * pow(1. + z, 4.));
    return Om_z;}

static double
get_delta_mean_from_delta_crit_at_z(hmpdf_obj *d, double delta_crit,
                                           double z){
    
    //get_delta_mean_from_delta_crit_at_z function in class_sz_tools.c https://github.com/CLASS-SZ/class_sz/blob/master/class-sz/tools/class_sz_tools.c
    
    double Omega_m_z = get_Omega_m_nonu_at_z(d,z);//get matter density
    double delta_mean = delta_crit / Omega_m_z; //compute \Delta_m from \Delta_c 

    return delta_mean;
}


static double
fnu_Tinker08(hmpdf_obj *d, double sigma, double z){

      //adapted from class_sz: MF_T08_m500 function in class_sz_tools.c https://github.com/CLASS-SZ/class_sz/blob/master/class-sz/tools/class_sz_tools.c

      z = (z<2.5) ? z : 2.5;
      
      double delta_mean_not_log;

      //!NOTE HARDCODED \Delta=200!
      if (d->h->HMF_mdef==hmpdf_mdef_c){
      
      //get delta_mean from delta_crit
      delta_mean_not_log=get_delta_mean_from_delta_crit_at_z(d, 200.,z);}

      else if (d->h->HMF_mdef==hmpdf_mdef_m){
    
      delta_mean_not_log=200.0;}
    
      double *A_z0=malloc(sizeof(double));
      double *a_z0=malloc(sizeof(double));
      double *b_z0=malloc(sizeof(double));
      double *c_z0=malloc(sizeof(double));      
  
      if (d->h->interp_T08==0){
      //params from the paper, Table 2, row 1 for \Delta = 200
      *A_z0=0.186;
      *a_z0=1.47;
      *b_z0=2.57;
      *c_z0=1.19;
      }
      else if (d->h->interp_T08>0){
      
      //interpolate
      
      double delta_mean = log10(delta_mean_not_log); //interp in log
         
      double delta_mean_arr[9]={200., 300., 400., 600., 800., 1200., 1600., 2400., 3200.};
      
      int i;
      for (i=0;i<9;i++)//interp in log
      delta_mean_arr[i] =
      log10(delta_mean_arr[i]);

      //Table 2 in Tinker+08
      double A_arr[9]={0.186,0.200, 0.212, 0.218, 0.248, 0.255, 0.260, 0.260, 0.260};
      double a_arr[9]={1.47, 1.52, 1.56, 1.61, 1.87, 2.13, 2.30, 2.53, 2.66};
      double b_arr[9]={2.57, 2.25, 2.05, 1.87, 1.59, 1.51, 1.46, 1.44, 1.41};
      double c_arr[9]={1.19, 1.27, 1.34, 1.45, 1.58, 1.80, 1.97, 2.24, 2.44};
      double d2_A_arr[9]={0.00, 0.50, -1.56, 3.05, -2.95, 1.07, -0.71, 0.21, 0.00};
      double d2_a_arr[9]={0.00,1.19,-6.34,21.36,-10.95,2.59,-0.85,-2.07,0.00};
      double d2_b_arr[9]={0.00, -1.08, 12.61,-20.96,24.08, -6.64, 3.84, -2.09,0.00};
      double d2_c_arr[9]={0.00, 0.94, -0.43, 4.61, 0.01, 1.21, 1.43, 0.33, 0.00};
        
      //interp values
      splint(delta_mean_arr, A_arr,d2_A_arr,9,delta_mean,A_z0);
      splint(delta_mean_arr, a_arr,d2_a_arr,9,delta_mean,a_z0);
      splint(delta_mean_arr, b_arr,d2_b_arr,9,delta_mean,b_z0);
      splint(delta_mean_arr,c_arr,d2_c_arr,9,delta_mean,c_z0);
      }
        
      //compute at z
      double alphaT08 =pow(10.,-pow(0.75/log10(delta_mean_not_log/75.),1.2));
      double A=*A_z0*pow(1.+z,-0.14);
      double a=*a_z0*pow(1.+z,-0.06);
      double b=*b_z0*pow(1.+z,-alphaT08);
      double c=*c_z0;
      
      free(A_z0);
      free(a_z0);
      free(b_z0);
      free(c_z0);

      return (A*(pow(sigma/b,-a)+1.)*exp(-c/pow(sigma,2.)));
}

static int
dndlogM(hmpdf_obj *d, int z_index, int M_index, double *hmf, double *bias)
{//{{{
    STARTFCT

    // case where we are above the possible mass cut
    if (d->n->mass_cuts != NULL
        && d->n->Mgrid[M_index]
           > d->n->mass_cuts(d->n->zgrid[z_index], d->n->mass_cuts_params)
             / d->c->h)
    {
        *hmf = 0.0;
        *bias = 0.0;
    }
    // we are below the mass cut (or none is given)
    else
    {  //original, using T10 mass function
        if (strcmp(d->h->HMF_func, "T10")==0.0){
        double sigma_squared = d->pwr->ssq[M_index][0];
        double sigma_squared_prime = d->pwr->ssq[M_index][1];
        //double nu = 1.686/sqrt(d->c->Dsq[z_index] * sigma_squared);
    	double dc=3.0/20.0*pow(12.0*M_PI,2.0/3.0);
        double nu=dc/sqrt(d->c->Dsq[z_index]*sigma_squared);
        double fnu = fnu_Tinker10(d, nu, d->n->zgrid[z_index]);

        //#ifdef SAVE_SIGMA_NU
	    //FILE *fp = fopen("/scratch/07833/tg871330/tSZ_maps/hmpdf_maps/sigma_nu/hmf_sigma_nu.txt", "a");
	    //fprintf(fp, "%.8f %.18e %.18e %.18e\n",d->n->zgrid[z_index], d->n->Mgrid[M_index],sigma_squared*d->c->Dsq[z_index], nu);
	    //fclose(fp);
	    //#endif

        *hmf = -fnu * d->c->rho_m_0 * sigma_squared_prime
               / (2.0 * sigma_squared * d->n->Mgrid[M_index]);

        *bias = bnu_Tinker10(nu);}
        
        else if (strcmp(d->h->HMF_func, "T08")==0.0){
        //using T08 mass function at M200c
        double sigma_squared = d->pwr->ssq[M_index][0];
        double sigma=sqrt(d->c->Dsq[z_index]*d->pwr->ssq[M_index][0]);
        double sigma_squared_prime = d->pwr->ssq[M_index][1];
        double dc=3.0/20.0*pow(12.0*M_PI,2.0/3.0);
        double nu=dc/sqrt(d->c->Dsq[z_index]*sigma_squared);
        double fnu = fnu_Tinker08(d, sigma, d->n->zgrid[z_index]);

        *hmf = -fnu * d->c->rho_m_0 * sigma_squared_prime
               / (2.0 * sigma_squared * d->n->Mgrid[M_index]);

        *bias = bnu_Tinker10(nu);
        
        }
        
        else {
        HMPDFERR("Pick mass function between T10 or T08");
 
        }
        if (d->h->massfunc_corr != NULL)
        {
            *hmf *= d->h->massfunc_corr(d->n->zgrid[z_index],
                                        d->n->Mgrid[M_index] * d->c->h,
                                        d->h->massfunc_corr_params);
        }

        if (d->h->bias_resc != NULL)
        {
            *bias *= d->h->bias_resc(d->n->zgrid[z_index],
                                     d->n->Mgrid[M_index] * d->c->h,
                                     d->h->bias_resc_params);
        }
    }

    ENDFCT
}//}}}

static int
create_dndlogM(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\tcreate_dndlogM\n");

    SAFEALLOC(d->h->hmf,  malloc(d->n->Nz * sizeof(double *)));
    SETARRNULL(d->h->hmf, d->n->Nz);
    SAFEALLOC(d->h->bias, malloc(d->n->Nz * sizeof(double *)));
    SETARRNULL(d->h->hmf, d->n->Nz);

//    #ifdef SAVE_SIGMA_NU
//    FILE *fp = fopen("/scratch/07833/tg871330/tSZ_maps/hmpdf_maps/sigma_nu/hmf_sigma_nu.txt", "w");
//    fclose(fp);
//    #endif

    for (int z_index=0; z_index<d->n->Nz; z_index++)
    {
        SAFEALLOC(d->h->hmf[z_index],  malloc(d->n->NM * sizeof(double)));
        SAFEALLOC(d->h->bias[z_index], malloc(d->n->NM * sizeof(double)));
        for (int M_index=0; M_index<d->n->NM; M_index++)
        {
            SAFEHMPDF(dndlogM(d, z_index, M_index,
                              d->h->hmf[z_index]+M_index,
                              d->h->bias[z_index]+M_index));
        }
    
    if (strcmp(d->h->hmf_file, "none")!=0){ 
    //printf("%s hmf file\n",d->h->hmf_file);
    double hmf_input[d->n->Nz * d->n->NM];
    FILE *fp_hmf=fopen(d->h->hmf_file, "rb");
    fread(hmf_input, sizeof(double), d->n->Nz * d->n->NM, fp_hmf);
    }
    

	#ifdef SAVE_HMF
	char buffer[512];
    sprintf(buffer, "/scratch/07833/tg871330/software_scratch/hmpdf/hmf/hmf_%.8f.bin", d->n->zgrid[z_index]);
    FILE *fp = fopen(buffer, "w");
	fwrite(d->n->Mgrid, sizeof(double), d->n->NM, fp);
	fwrite(d->h->hmf[z_index], sizeof(double), d->n->NM, fp);
	fclose(fp);
	#endif
    }

    ENDFCT
}//}}}

int
init_halo_model(hmpdf_obj *d)
{//{{{
    STARTFCT

    if (d->h->inited_halo) { return 0; }

    HMPDFPRINT(1, "init_halo_model\n");

    SAFEHMPDF(create_c_of_y(d));
    SAFEHMPDF(create_dndlogM(d));
    d->h->inited_halo = 1;

    ENDFCT
}//}}}

