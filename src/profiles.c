#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_dht.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_fit.h>

#include "configs.h"
#include "utils.h"
#include "data.h"
#include "cosmology.h"
#include "halo_model.h"
#include "filter.h"
#include "profiles.h"

static
void create_angle_grids(struct all_data_s *d)
{//{{{
    printf("\tcreate_angle_grids\n");
    d->p->prtilde_Ntheta = PRTILDE_INTEGR_NTHETA;

    d->p->decr_tgrid = (double *)malloc(d->p->Ntheta * sizeof(double));
    d->p->incr_tgrid = (double *)malloc(d->p->Ntheta * sizeof(double));
    d->p->decr_tsqgrid = (double *)malloc(d->p->Ntheta * sizeof(double));
    d->p->reci_tgrid = (double *)malloc(d->p->Ntheta * sizeof(double));
    d->p->prtilde_thetagrid = (double *)malloc(d->p->prtilde_Ntheta * sizeof(double));
    for (int ii=0; ii<d->p->Ntheta; ii++)
    {
        // reverse order, maximum angle is the first one
        d->p->decr_tgrid[d->p->Ntheta-1-ii] = gsl_sf_bessel_zero_J0(ii+1)/gsl_sf_bessel_zero_J0(d->p->Ntheta);
        d->p->decr_tsqgrid[d->p->Ntheta-1-ii] = gsl_pow_2(d->p->decr_tgrid[d->p->Ntheta-1-ii]);

        d->p->reci_tgrid[ii] = gsl_sf_bessel_zero_J0(ii+1);
    }
    reverse(d->p->Ntheta, d->p->decr_tgrid, d->p->incr_tgrid);
    linspace(d->p->prtilde_Ntheta, 0.0, 1.0, d->p->prtilde_thetagrid);
    d->p->incr_tgrid_accel = gsl_interp_accel_alloc();
    d->p->reci_tgrid_accel = gsl_interp_accel_alloc();
}//}}}

static
void kappa_profile(all_data *d, int z_index, int M_index,
                   double theta_out, double Rout, double *p)
{//{{{
    // find the NFW parameters
    double rhos, rs;
    rhos = NFW_fundamental(d, z_index, M_index, &rs);

    // fill the profile
    for (int ii=0; ii<d->p->Ntheta; ii++)
    {
        double t = d->p->decr_tgrid[ii] * theta_out;
        double Rproj = tan(t) * d->c->angular_diameter[z_index];
        double lin = 0.0;
        double lout = sqrt(Rout*Rout - Rproj*Rproj);

        if (Rproj > rs)
        {
            p[ii] = 2.0*rhos*gsl_pow_3(rs)/(2.0*(Rout+rs)*pow((Rproj-rs)*(Rproj+rs),1.5))
                            *(+M_PI*rs*(Rout+rs)
                              +2.0*sqrt((Rout-Rproj)*(Rout+Rproj)*(Rproj-rs)*(Rproj+rs))
                              -2.0*rs*(Rout+rs)*atan2(Rproj*Rproj+Rout*rs,                
                                                      -sqrt((Rout-Rproj)*(Rout+Rproj)          
                                                            *(Rproj-rs)*(Rproj+rs))));
        }
        else if (Rproj < rs)
        {
            p[ii] = 2.0*rhos*gsl_pow_3(rs)/((Rout+rs)*pow(rs*rs-Rproj*Rproj,1.5))
                            *(-sqrt((Rout-Rproj)*(Rout+Rproj)*(rs*rs-Rproj*Rproj))
                              +rs*(Rout+rs)*log(Rproj*(Rout+rs)/(Rproj*Rproj+Rout*rs
                                                                 -sqrt((Rout-Rproj)*(Rout+Rproj)
                                                                       *(rs*rs-Rproj*Rproj)))));
        }
        else
        {
            p[ii] = 2.0*rhos*sqrt(Rout-rs)*rs*(Rout+2.0*rs)/(3.0*pow(Rout+rs,1.5));
        }
        // TODO check if rho_m here physical or comoving
        p[ii] -= 2.0*(lout-lin)*d->c->rho_m[z_index];
        p[ii] /= d->c->Scrit[z_index];
    }
    p[0] = 0.0; // FIXME without this there is a bug
}//}}}

// Battaglia profiles{{{
typedef struct
{
    double alpha;
    double beta;
    double gamma;
    double rproj;
}
Battmodel_params;
                                  //   A0    alpham  alphaz
double Battmodel_fitparams[][3] = {{ 18.1,    0.154, -0.758}, // P0
                                   {0.497, -0.00865,  0.731}, // xc
                                   {  1.0,      0.0,    0.0}, // alpha
                                   { 4.35,   0.0393,  0.415}, // beta
                                   { -0.3,      0.0,    0.0}};// gamma
static
double Battmodel_primitive(double M200c, double z, int numparam)
{
    return Battmodel_fitparams[numparam][0]
           * pow(M200c/1e14, Battmodel_fitparams[numparam][1])
           * pow(1.0+z, Battmodel_fitparams[numparam][2]);
}
static
double Battmodel_integrand(double z, void *params)
{
    Battmodel_params *p = (Battmodel_params *)params;
    double r = hypot(z, p->rproj);
    return pow(r, p->gamma) / pow(1.0 + pow(r, p->alpha), p->beta);
}
static
void tsz_profile(all_data *d, int z_index, int M_index,
                 double theta_out, double Rout, double *p)
{
    // convert to 200c
    double M200c, R200c, c200c;
    M200c = Mconv(d, z_index, M_index, mdef_c, &R200c, &c200c);
    double P0 = Battmodel_primitive(M200c, d->n->gr->zgrid[z_index], 0);
    double xc = Battmodel_primitive(M200c, d->n->gr->zgrid[z_index], 1);
    Rout /= R200c * xc;
    // prepare the integration
    Battmodel_params par;
    par.alpha = Battmodel_primitive(M200c, d->n->gr->zgrid[z_index], 2);
    par.beta  = Battmodel_primitive(M200c, d->n->gr->zgrid[z_index], 3);
    par.gamma = Battmodel_primitive(M200c, d->n->gr->zgrid[z_index], 4);
    gsl_function integrand;
    integrand.function = &Battmodel_integrand;
    integrand.params = &par;
    // loop over angles
    for (int ii=0; ii<d->p->Ntheta; ii++)
    {
        double t = d->p->decr_tgrid[ii] * theta_out;
        par.rproj = tan(t) * d->c->angular_diameter[z_index] / R200c / xc;
        double lin  = 0.0;
        double lout = sqrt(Rout*Rout - par.rproj*par.rproj);
        double err;
        size_t neval;
        gsl_integration_qng(&integrand, lin, lout,
                            BATTINTEGR_EPSABS, BATTINTEGR_EPSREL,
                            p+ii, &err, &neval);
        // normalize
        p[ii] *= P0 * xc * M200c * 200.0
                 // TODO check if rho_c here comoving or physical
                 * d->c->rho_c[z_index] * d->c->Ob_0/d->c->Om_0
                 * GNEWTON * SIGMATHOMSON / MELECTRON / gsl_pow_2(SPEEDOFLIGHT)
                 / 1.932; // convert from thermal to electron pressure
    }
}
//}}}

static
double profile(all_data *d, int z_index, int M_index, double *p)
// returns theta_out and writes the profile into return value
{//{{{
    // find the outer radius on the sky
    double Rout, c;
    Mconv(d, z_index, M_index, d->p->rout_def, &Rout, &c);
    Rout *= d->p->rout_scale;
    double theta_out = atan(Rout/d->c->angular_diameter[z_index]);

    if (d->p->stype == kappa)
    {
        kappa_profile(d, z_index, M_index, theta_out, Rout, p);
    }
    else if (d->p->stype == tsz)
    {
        tsz_profile(d, z_index, M_index, theta_out, Rout, p);
    }

    return theta_out;
}//}}}

static
void create_profiles(all_data *d)
{//{{{
    printf("\tcreate_profiles\n");
    d->p->profiles = (double ***)malloc(d->n->gr->Nz * sizeof(double **));
    for (int z_index=0; z_index<d->n->gr->Nz; z_index++)
    {
        d->p->profiles[z_index] = (double **)malloc(d->n->gr->NM * sizeof(double *));
        for (int M_index=0; M_index<d->n->gr->NM; M_index++)
        {
            d->p->profiles[z_index][M_index] = (double *)malloc((d->p->Ntheta+1) * sizeof(double));
            d->p->profiles[z_index][M_index][0] = profile(d, z_index, M_index,
                                                          d->p->profiles[z_index][M_index]+1);
        }
    }
}//}}}

void create_conj_profiles(all_data *d)
// computes the conjugate space profiles
{//{{{
    if (d->p->created_conj_profiles) { return; }
    printf("\tcreate_conj_profiles\n");
    // prepare the Hankel transform work space
    d->p->dht_ws = gsl_dht_new(d->p->Ntheta, 0, 1.0);
    // need buffer to store the profiles with theta increasing
    double *temp = (double *)malloc(d->p->Ntheta * sizeof(double));
    d->p->conj_profiles = (double ***)malloc(d->n->gr->Nz * sizeof(double **));
    for (int z_index=0; z_index<d->n->gr->Nz; z_index++)
    {
        d->p->conj_profiles[z_index] = (double **)malloc(d->n->gr->NM * sizeof(double *));
        for (int M_index=0; M_index<d->n->gr->NM; M_index++)
        {
            d->p->conj_profiles[z_index][M_index] = (double *)malloc((d->p->Ntheta+1) * sizeof(double));
            reverse(d->p->Ntheta, d->p->profiles[z_index][M_index]+1, temp);
            gsl_dht_apply(d->p->dht_ws, temp, d->p->conj_profiles[z_index][M_index]+1);
            d->p->conj_profiles[z_index][M_index][0] = 1.0/d->p->profiles[z_index][M_index][0];
        }
    }
    free(temp);

    d->p->created_conj_profiles = 1;
}//}}}

static
void create_filtered_profiles(all_data *d)
{//{{{
    if (d->f->Nfilters == 0 ) { return; }
    printf("\tcreate_filtered_profiles\n");

    double *ell = (double *)malloc(d->p->Ntheta * sizeof(double));
    double *temp = (double *)malloc(d->p->Ntheta * sizeof(double)); // buffer
    for (int z_index=0; z_index<d->n->gr->Nz; z_index++)
    {
        for (int M_index=0; M_index<d->n->gr->NM; M_index++)
        {
            for (int ii=0; ii<d->p->Ntheta; ii++)
            {
                ell[ii] = d->p->reci_tgrid[ii] * d->p->conj_profiles[z_index][M_index][0];
            }
            apply_filters(d, d->p->Ntheta, ell, d->p->conj_profiles[z_index][M_index]+1, temp, filter_pdf);
            gsl_dht_apply(d->p->dht_ws, temp, d->p->profiles[z_index][M_index]+1);
            // reverse the profile
            reverse(d->p->Ntheta, d->p->profiles[z_index][M_index]+1, d->p->profiles[z_index][M_index]+1);
            // normalize properly
            for (int ii=0; ii<d->p->Ntheta; ii++)
            {
                d->p->profiles[z_index][M_index][ii+1] *= gsl_pow_2(d->p->reci_tgrid[d->p->Ntheta-1]);
            }
            // set the value at maximum theta to zero
            d->p->profiles[z_index][M_index][1] = 0.0;
        }
    }
    free(temp);
    free(ell);
}//}}}

static
int not_monotonic(int N, double *x, int *problems)
// checks if x is monotonically increasing
// (this is what we require for the FFTs to work)
{//{{{
    int out = 0;
    for (int ii=1; ii<N; ii++)
    {
        if (x[ii] < x[ii-1])
        {
            problems[out] = ii;
            out += 1;
        }
    }
    return out;
}//}}}

static
int find_breakpoint(int N, int *x)
// computes the first index from which we can apply the FFTs up to Mmax
{//{{{
    int out = 0;
    for (int ii=0; ii<N; ii++)
    {
        if (x[ii]) // not monotonic
        {
            out = ii + 1;
        }
    }
    return out;
}//}}}

static
void adjust_breakpoint(all_data *d, int *bp)
{//{{{
    if (GSL_IS_ODD(*bp))
    // we want to do simpson integrations with odd number of sample points
    {
        *bp += 1;
    }
    if (*bp < INTEGR_MINSAMPLES-1)
    // do more slow integrals
    {
        *bp = INTEGR_MINSAMPLES-1;
    }
    if (d->n->gr->NM - *bp < INTEGR_MINSAMPLES-1)
    // remove the small FFT integral
    {
        *bp = d->n->gr->NM-1;
    }
}//}}}

static
void monotonize(int Nx, int Nproblems, int *problems, double *x, double *y)
// problems is an int[Nproblems], holding the indices where y was decreasing,
// in increasing order
{//{{{
    int ledge = 0;
    int redge = 0;
    for (int ii=0; ii<Nproblems; ii++)
    {
        // check if we have already covered this problem
        if (redge > problems[ii])
        {
            continue;
        }
        ledge = problems[ii] - 1; // >= 0
        
        // find the right edge
        redge = Nx;
        for (int jj=ledge+1; jj<Nx; jj++)
        {
            if (y[jj] > y[ledge])
            {
                redge = jj;
                break;
            }
        }

        double a, b;
        // treat the case when no redge was found (only decreasing from here onwards)
        // in this case, we build a simple linear ramp, and hope that this profile is not important
        if (redge == Nx)
        {
            if (ledge == 0)
            // case when we have no valid points
            // use stdev as relevant scale
            // profiles using this are completely irrelevant
            {
                double scale = 0.0;
                for (int jj=1; jj<Nx; jj++)
                {
                    scale += y[jj] * y[jj];
                }
                scale = sqrt(scale/(double)(Nx-1));
                a = scale / (x[Nx-1] - x[ledge]);
                b = - a * x[ledge];
            }
            else
            // perform linear fit through the last few datapoints
            // x       x      x      x    x
            // _start ...   ....   ....  ledge
            {
                int _start = GSL_MAX(0, ledge-3);
                double cov00, cov01, cov11, sumsq;
                gsl_fit_linear(x+_start, 1, y+_start, 1, ledge-_start+1,
                               &b, &a, &cov00, &cov01, &cov11, &sumsq);
                // as is by construction positive, but we still need to ensure
                // b doesn't spoil the party
                b = GSL_MAX(y[ledge] - a * x[ledge], b);
            }
        }
        else
        {
            a = (y[redge]-y[ledge]) / (x[redge] - x[ledge]);
            b = y[ledge] - a * x[ledge];
        }

        for (int jj=ledge+1; jj<redge; jj++)
        {
            y[jj] = a * x[jj] + b;
        }
    }
}//}}}

static
void create_breakpoints_or_monotonize(all_data *d)
{//{{{
    printf("\tcreate_breakpoints\n");
    d->p->breakpoints = (int *)malloc(d->n->gr->Nz * sizeof(int));
    for (int z_index=0; z_index<d->n->gr->Nz; z_index++)
    {
        int _not_monotonic[d->n->gr->NM];
        for (int M_index=0; M_index<d->n->gr->NM; M_index++)
        {
            int _problems[d->p->Ntheta];
            _not_monotonic[M_index] = not_monotonic(d->p->Ntheta,
                                                    d->p->profiles[z_index][M_index]+1,
                                                    _problems);
            if (d->n->monotonize)
            {
                monotonize(d->p->Ntheta, _not_monotonic[M_index], _problems,
                           d->p->decr_tgrid, d->p->profiles[z_index][M_index]+1);
                _not_monotonic[M_index] = 0;
            }
        }

        d->p->breakpoints[z_index] = find_breakpoint(d->n->gr->NM, _not_monotonic);
        if (d->p->breakpoints[z_index] > 0)
        // we potentially need to make adjustments to breakpoint
        {
            adjust_breakpoint(d, d->p->breakpoints+z_index);
        }
    }
}//}}}

int breakpoints(all_data *d, int z_index)
// TODO remove this function
{//{{{
    if (d->n->monotonize)
    {
        return 0;
    }
    else
    {
        return d->p->breakpoints[z_index];
    }
}//}}}

void s_of_t(all_data *d, int z_index, int M_index, int Nt, double *t, double *s)
// returns signal(t) at z_index, M_index
{//{{{
    double *temp = (double *)malloc(d->p->Ntheta * sizeof(double));
    reverse(d->p->Ntheta, d->p->profiles[z_index][M_index]+1, temp);
    interp1d *interp = new_interp1d(d->p->Ntheta, d->p->incr_tgrid, temp, temp[0], 0.0,
                                    PRINTERP_TYPE, d->p->incr_tgrid_accel);

    for (int ii=0; ii<Nt; ii++)
    {
        s[ii] = interp1d_eval(interp, t[ii]);
    }
    free(temp);
    delete_interp1d(interp);
}//}}}

void s_of_l(all_data *d, int z_index, int M_index, int Nl, double *l, double *s)
// returns the conjugate space profile at z_index, M_index
{//{{{
    interp1d *interp = new_interp1d(d->p->Ntheta, d->p->reci_tgrid,
                                    d->p->conj_profiles[z_index][M_index]+1,
                                    d->p->conj_profiles[z_index][M_index][1]/*high l*/, 0.0/*low l*/,
                                    SELL_INTERP_TYPE, d->p->reci_tgrid_accel);
    double hankel_norm = 0.0; // FIXME
    for (int ii=0; ii<Nl; ii++)
    {
        double l_normalized = l[ii] / d->p->conj_profiles[z_index][M_index][0];
        s[ii] = hankel_norm * interp1d_eval(interp, l_normalized);
    }
    delete_interp1d(interp);
}//}}}

void dtsq_of_s(all_data *d, int z_index, int M_index, double *dtsq)
// write dtheta(signal)/dsignal*dsignal into return values
{//{{{
    interp1d *interp = new_interp1d(d->p->Ntheta, d->p->profiles[z_index][M_index]+1,
                                    d->p->decr_tsqgrid, 0.0, 0.0, PRINTERP_TYPE, NULL);
    for (int ii=0; ii<d->n->gr->Nsignal; ii++)
    {
        dtsq[ii] = gsl_pow_2(d->p->profiles[z_index][M_index][0])
                   * (d->n->gr->signalgrid[1] - d->n->gr->signalgrid[0])
                   * interp1d_eval_deriv(interp, d->n->gr->signalgrid[ii]);
    }
    delete_interp1d(interp);
}//}}}

void t_of_s(all_data *d, int z_index, int M_index, double *t)
// only valid results if monotonize is True
{//{{{
    interp1d *interp = new_interp1d(d->p->Ntheta, d->p->profiles[z_index][M_index]+1,
                                    d->p->decr_tgrid, 0.0, 0.0, PRINTERP_TYPE, NULL);
    for (int ii=0; ii<d->n->gr->Nsignal; ii++)
    {
        t[ii] = d->p->profiles[z_index][M_index][0]
                * interp1d_eval(interp, d->n->gr->signalgrid[ii]);
    }
    delete_interp1d(interp);
}//}}}

void init_profiles(all_data *d)
{//{{{
    if (d->p->inited_profiles) { return; }
    printf("In profiles.h -> init_profiles :\n");
    create_angle_grids(d);
    create_profiles(d);

    if (d->f->Nfilters > 0)
    {
        create_conj_profiles(d);
        create_filtered_profiles(d);
    }
    create_breakpoints_or_monotonize(d); // have to keep this here!

    d->p->inited_profiles = 1;
}//}}}
