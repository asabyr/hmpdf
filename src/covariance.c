#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>

#include "utils.h"
#include "configs.h"
#include "object.h"
#include "power.h"
#include "profiles.h"
#include "noise.h"
#include "onepoint.h"
#include "twopoint.h"
#include "covariance.h"

#include "hmpdf.h"

int
null_covariance(hmpdf_obj *d)
{//{{{{
    STARTFCT

    d->cov->ws = NULL;
    d->cov->Cov = NULL;
    d->cov->Cov_noisy = NULL;
    d->cov->corr_diagn = NULL;
    d->cov->created_tp_ws = 0;
    d->cov->created_phigrid = 0;
    d->cov->created_cov = 0;
    d->cov->created_noisy_cov = 0;

    ENDFCT
}//}}}

int
reset_covariance(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\treset_covariance\n");

    if (d->cov->Cov != NULL) { free(d->cov->Cov); }
    if (d->cov->Cov_noisy != NULL) { free(d->cov->Cov_noisy); }
    if (d->cov->corr_diagn != NULL) { free(d->cov->corr_diagn); }
    if (d->cov->ws != NULL)
    {
        for (int ii=0; ii<d->cov->Nws; ii++)
        {
            if (d->cov->ws[ii] != NULL)
            {
                fftw_free(d->cov->ws[ii]->pdf_real);
                free(d->cov->ws[ii]->bc);
                fftw_free(d->cov->ws[ii]->tempc_real);
                fftw_destroy_plan(d->cov->ws[ii]->pu_r2c);
                fftw_destroy_plan(d->cov->ws[ii]->pc_r2c);
                fftw_destroy_plan(d->cov->ws[ii]->ppdf_c2r);
                free(d->cov->ws[ii]);
            }
        }
        free(d->cov->ws);
    }

    ENDFCT
}//}}}

static int
comp_int(const void *a, const void *b)
// to qsort an array of ints
{//{{{
    int int_a = *((int *)(a));
    int int_b = *((int *)(b));
    if (int_a == int_b)
    {
        return 0;
    }
    else if (int_a > int_b)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}//}}}

static int
phigrid_exact_part_centers(hmpdf_obj *d, double **grid, double **weights, int *buflen, int *Nexact)
{//{{{
    STARTFCT

    int *_rsq;
    SAFEALLOC(_rsq, malloc(d->n->pixelexactmax * d->n->pixelexactmax
                           * sizeof(int)));
    int N = 0; // counts possibly repeated sums of two squares
    // loop over one quadrant
    for (int ii=0; ii<=d->n->pixelexactmax; ii++)
    {
        for (int jj=1;
             jj*jj <= d->n->pixelexactmax * d->n->pixelexactmax - ii*ii;
             jj++)
        {
            _rsq[N++] = ii*ii + jj*jj;
        }
    }
    // sort and identify unique values
    qsort(_rsq, N, sizeof(int), comp_int);
    *Nexact = 0; // counts number of phi-values for which exact treatment is needed
    int realloced = 0;
    for (int ii=0; ii<N;)
    {
        int ctr = 0;
        int jj;
        for (jj=ii; _rsq[jj] == _rsq[ii]; jj++)
        {
            ctr += 1;
        }

        while (UNLIKELY(*Nexact >= *buflen))
        {
            *buflen *= 2;
            SAFEALLOC(*grid,    realloc(*grid,    *buflen * sizeof(double)));
            SAFEALLOC(*weights, realloc(*weights, *buflen * sizeof(double)));
            ++realloced;
            if (realloced > 2)
            {
                HMPDFERR("Failed to create phigrid, "
                         "expanded buffer too often");
            }
        }

        *grid[*Nexact] = sqrt((double)(_rsq[ii]))
                              * d->f->pixelside;
        *weights[(*Nexact)++] = (double)(4 * ctr);
        ii = jj;
    }
    free(_rsq);

    ENDFCT
}//}}}

static int
phigrid_exact_part_jitter(hmpdf_obj *d, double **grid, double **weights, int *buflen, int *Nexact)
{//{{{
    STARTFCT
    //
    // use a falling density of sample points
    // we transform phi = x^phipwr and sample uniformly in x
    // then the "density of states" is
    //      rho(phi) = rho0 * phi^(1/a - 1)
    // We want to fix the constant rho0 such that \int_phimin^phimax rho(phi) dphi = Nphi
    double integral = d->n->phipwr
                      * (pow(d->n->phimax, 1.0/d->n->phipwr)
                         - pow(d->f->pixelside, 1.0/d->n->phipwr));
    double rho0 = (double)(d->n->Nphi) / integral;
    int end = *Nexact; // index of the first free element in phigrid
    int realloced = 0;
    for (int ii=0; ii<*Nexact; ii++)
    {
        double deltaphi;
        if (ii == 0)
        {
            deltaphi = grid[1] - grid[0];
        }
        else if (ii == *Nexact-1)
        {
            deltaphi = grid[*Nexact-1] - grid[*Nexact-2];
        }
        else
        {
            deltaphi = 0.5 * (grid[ii+1] - grid[ii-1]);
        }
        double rho_here = pow(*grid[ii], 1.0/d->n->phipwr-1.0);
        
        int Nphi_here = (int)(ceil(deltaphi * rho0 * rho_here)) - 1;

        if (Nphi_here%2)
        // want Nphi_here to be even for convenience
        // also this averages out the ceil
        {
            --Nphi_here;
        }

        *weights[ii] /= (double)(Nphi_here+1);

        while (UNLIKELY(end+Nphi_here-1 >= *buflen))
        {
            *buflen *= 2;
            SAFEALLOC(*grid,    realloc(*grid,    *buflen * sizeof(double)));
            SAFEALLOC(*weights, realloc(*weights, *buflen * sizeof(double)));
            ++realloced;
            if (realloced > 2)
            {
                HMPDFERR("Failed to create phigrid, "
                         "expanded buffer too often.")
            }
        }

        for (int jj=1; jj<=Nphi_here/2; jj++)
        {
            *grid[end+2*jj-2] = *grid[ii]
                                + deltaphi * d->n->phijitter
                                  * (double)(jj) / (double)(Nphi_here/2);
            *grid[end+2*jj-1] = *grid[ii]
                                - deltaphi * d->n->phijitter
                                  * (double)(jj) / (double)(Nphi_here/2);
            *weights[end+2*jj-2] = *weights[end+2*jj-1] = *weights[ii];
        }
        end += Nphi_here;
    }
    *Nexact = end; // Nphi - Nexact is the number of continuum phi values

    ENDFCT
}//}}}

static int
phigrid_approx_part(hmpdf_obj *d, int Nexact, double **grid, double **weights, int *buflen, int *Napprox)
{//{{{
    STARTFCT

    // fill the rest of the phi-grid
    // the continuum version is
    //      \int_lo^\phi_max d\phi 2\pi\phi/(pixel sidelength)^2 * integrand
    double lo = pow((double)(d->n->pixelexactmax) * d->f->pixelside,
                    1.0/d->n->phipwr);
    double hi = pow(d->n->phimax, 1.0/d->n->phipwr);

    gsl_integration_fixed_workspace *t;
    SAFEALLOC(t, gsl_integration_fixed_alloc(gsl_integration_fixed_legendre,
                                             d->n->Nphi-Nexact,
                                             lo, hi, 0.0, 0.0));
    
    double *x = gsl_integration_fixed_nodes(t);
    double *w = gsl_integration_fixed_weights(t);
    
    *Napprox = 0;
    int realloced = 0;
    for (int ii=0; ii<(int)gsl_integration_fixed_n(t); ii++)
    {
        if (x[ii] < lo)
        {
            continue;
        }
        else if (x[ii] > hi)
        // nodes are sorted
        {
            break;
        }
        else
        {
            while (UNLIKELY(Nexact + *Napprox >= *buflen))
            {
                *buflen *= 2;
                SAFEALLOC(*grid,    realloc(*grid,    *buflen * sizeof(double)));
                SAFEALLOC(*weights, realloc(*weights, *buflen * sizeof(double)));
                if (realloced > 2)
                {
                    HMPDFERR("Failed to create phigrid, "
                             "expanded buffer too often.");
                }
            }

            *grid[Nexact + *Napprox] = pow(x[ii], d->n->phipwr);
            // fix the normalization
            *weights[Nexact + *Napprox] = w[ii]
                                          * 2.0 * M_PI * *grid[Nexact + *Napprox]          
                                          / gsl_pow_2(d->f->pixelside)                 
                                          // Jacobian of the transformation            
                                          * d->n->phipwr                               
                                          * pow(x[ii], d->n->phipwr-1.0);          
            ++(*Napprox);
        }
    }

    gsl_integration_fixed_free(t);

    ENDFCT
}//}}}

static int
phigrid_shuffle_and_copy(hmpdf_obj *d, double *grid, double *weights)
{//{{{
    STARTFCT

    int *indices;
    SAFEALLOC(indices, malloc(d->n->Nphi * sizeof(int)));

    for (int ii=0; ii<d->n->Nphi; ii++)
    {
        indices[ii] = ii;
    }

    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r, time(NULL));
    gsl_ran_shuffle(r, indices, d->n->Nphi, sizeof(int));
    gsl_rng_free(r);
    
    for (int ii=0; ii<d->n->Nphi; ii++)
    {
        d->n->phigrid[ii]    = grid[indices[ii]];
        d->n->phiweights[ii] = weights[indices[ii]];
    }

    free(indices);

    ENDFCT
}//}}}

static int
create_phigrid(hmpdf_obj *d)
{//{{{
    STARTFCT

    if (d->cov->created_phigrid) { return hmpdf_status; }

    HMPDFPRINT(2, "\tcreate_phigrid\n");

    // sanity checks
    if (d->n->pixelexactmax <= 0)
    {
        HMPDFERR("You set hmpdf_pixelexact_max=%d "
                 "(or did not set at all, "
                 "must be strictly positive.", d->n->pixelexactmax);
    }
    if (d->f->pixelside <= 0.0)
    {
        HMPDFERR("You set hmpdf_pixel_side=%.4e "
                 "(or did not set at all), "
                 "must be strictly positive.", d->f->pixelside);
    }

    // allocate twice as much as we probably need,
    // the continuum integral adds a bit of noise
    double *_phigrid;
    double *_phiweights;
    int buflen = 2 * d->n->Nphi;
    SAFEALLOC(_phigrid,    malloc(buflen * sizeof(double)));
    SAFEALLOC(_phiweights, malloc(buflen * sizeof(double)));

    // first treat the exact pixelization part
    int Nexact = 0; // to avoid maybe-uninitialized
    
    // find the exact pixel separations
    SAFEHMPDF(phigrid_exact_part_centers(d, &_phigrid, &_phiweights, &buflen, &Nexact));

    // add the jittered part for numerical stability
    SAFEHMPDF(phigrid_exact_part_jitter(d, &_phigrid, &_phiweights, &buflen, &Nexact));

    if (Nexact > d->n->Nphi)
    {
        HMPDFWARN("N_phi = %d is quite small "
                  "(suggested increase at least to %d)",
                  d->n->Nphi, 2*Nexact);
    }

    // now do the integral version for high pixel separations
    int Napprox = 0; // to avoid maybe-uninitialized
    SAFEHMPDF(phigrid_approx_part(d, Nexact, &_phigrid, &_phiweights, &buflen, &Napprox));

    // set Nphi to correct value
    d->n->Nphi = Nexact + Napprox;
    
    HMPDFPRINT(4, "\t\t\tNphi = %d, Nexact = %d\n", d->n->Nphi, Nexact);

    // now copy into the main grids
    // including shuffling to make parallel execution more efficient
    SAFEALLOC(d->n->phigrid,    malloc(d->n->Nphi * sizeof(double)));
    SAFEALLOC(d->n->phiweights, malloc(d->n->Nphi * sizeof(double)));

    SAFEHMPDF(phigrid_shuffle_and_copy(d, _phigrid, _phiweights));

    free(_phigrid);
    free(_phiweights);

    d->cov->created_phigrid = 1;

    ENDFCT
}//}}}

static int
create_tp_ws(hmpdf_obj *d)
{//{{{
    STARTFCT

    if (d->cov->created_tp_ws) { return hmpdf_status; }

    HMPDFPRINT(2, "\tcreate_tp_ws\n");
    HMPDFPRINT(3, "\t\ttrying to allocate workspaces for %d threads.\n", d->Ncores);
    
    SAFEALLOC(d->cov->ws, malloc(d->Ncores * sizeof(twopoint_workspace *)));
    d->cov->Nws = 0;
    // allocate workspaces until we run out of memory
    for (int ii=0; ii<d->Ncores; ii++)
    {
        SAFEHMPDF_NORETURN(new_tp_ws(d->n->Nsignal, d->cov->ws+ii));
        if (hmpdf_status) // failure to allocate a work space is not considered
                          // a critical error
        {
            hmpdf_status = 0;
            break;
        }
        else
        {
            ++d->cov->Nws;
        }
    }
    // NULL the remaining workspaces
    for (int ii=d->cov->Nws; ii<d->Ncores; ii++)
    {
        d->cov->ws[ii] = NULL;
    }

    if (d->cov->Nws < d->Ncores)
    {
        HMPDFPRINT(1, "Allocated only %d workspaces, "
                      "because memory ran out.\n", d->cov->Nws);
    }

    if (d->cov->Nws < 1)
    {
        HMPDFERR("Failed to allocate any workspaces.");
    }

    d->cov->created_tp_ws = 1;

    ENDFCT
}//}}}

static int
corr_diagn(hmpdf_obj *d, twopoint_workspace *ws, double *out)
{//{{{
    STARTFCT

    *out = 0.0;
    for (int ii=0; ii<d->n->Nsignal; ii++)
    {
        for (int jj=0; jj<d->n->Nsignal; jj++)
        {
            // assumes pdf_real to be properly normalized!
            *out += ws->pdf_real[ii*(d->n->Nsignal+2)+jj]
                    * (d->n->signalgrid[ii] - d->op->signalmeanc)
                    * (d->n->signalgrid[jj] - d->op->signalmeanc);
        }
    }

    ENDFCT
}//}}}

static int
status_update(hmpdf_obj *d, time_t t0, int done, int tot)
{//{{{
    STARTFCT

    time_t t1 = time(NULL);
    double delta_time = difftime(t1, t0);
    double remains = delta_time/(double)done
                     *(double)(tot - done);
    int hrs = (int)floor(remains/60.0/60.0);
    int min = (int)round(remains/60.0 - 60.0*(double)(hrs));
    int done_perc = (int)round(100.0*(double)(done)/(double)(tot));
    HMPDFPRINT(1, "\t\t%3d %% done, %.2d hrs %.2d min remaining "
                  "in create_cov.\n", done_perc, hrs, min);

    ENDFCT
}//}}}

static int
add_shotnoise_diag(int N, double *cov, double *p)
{//{{{
    STARTFCT

    for (int ii=0; ii<N; ii++)
    {
        cov[ii*(N+1)] += p[ii];
    }

    ENDFCT
}//}}}

static int
rescale_to_fsky1(hmpdf_obj *d, int N, double *cov)
// divides covariance matrix by number of pixels in the sky
{//{{{
    STARTFCT

    double Npixels = 4.0*M_PI/gsl_pow_2(d->f->pixelside);
    for (int ii=0; ii<N*N; ii++)
    {
        cov[ii] /= Npixels;
    }

    ENDFCT
}//}}}

static int
create_noisy_cov(hmpdf_obj *d)
{//{{{
    STARTFCT

    if (d->cov->created_noisy_cov) { return hmpdf_status; }

    HMPDFPRINT(2, "\tcreate_noisy_cov\n");

    SAFEALLOC(d->cov->Cov_noisy, malloc(d->n->Nsignal_noisy
                                        * d->n->Nsignal_noisy
                                        * sizeof(double)));
    SAFEHMPDF(noise_matr(d, d->cov->Cov, d->cov->Cov_noisy));

    d->cov->created_noisy_cov = 1;

    ENDFCT
}//}}}

static int
create_cov(hmpdf_obj *d)
{//{{{
    STARTFCT

    if (d->cov->created_cov) { return hmpdf_status; }

    HMPDFPRINT(2, "\tcreate_cov\n");

    // allocate storage
    SAFEALLOC(d->cov->Cov, malloc(d->n->Nsignal * d->n->Nsignal
                                  * sizeof(double)));
    SAFEALLOC(d->cov->corr_diagn, malloc(d->n->Nphi * sizeof(double)));

    // zero covariance
    zero_real(d->n->Nsignal * d->n->Nsignal, d->cov->Cov);

    // status
    int Nstatus = 0;
    time_t start_time = time(NULL);
    
    // loop over phi values
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(d->cov->Nws)
    #endif
    for (int pp=0; pp<d->n->Nphi; pp++)
    {
        if (hmpdf_status) { continue; }

        // create twopoint at this phi
        SAFEHMPDF_NORETURN(create_tp(d, d->n->phigrid[pp],
                                     d->cov->ws[this_core()]));

        if (hmpdf_status) { continue; }
        // compute the correlation function
        SAFEHMPDF_NORETURN(corr_diagn(d, d->cov->ws[this_core()],
                                      d->cov->corr_diagn+pp));

        if (hmpdf_status) { continue; }
        
        // add to covariance
        for (int ii=0; ii<d->n->Nsignal; ii++)
        {
            for (int jj=0; jj<d->n->Nsignal; jj++)
            {
                #ifdef _OPENMP
                // make sure no two threads add to the same element simultaneously
                #pragma omp atomic
                #endif
                d->cov->Cov[ii*d->n->Nsignal+jj]
                    += d->n->phiweights[pp]
                       * d->cov->ws[this_core()]->pdf_real[ii*(d->n->Nsignal+2)+jj];
            }
        }

        // status update
        #ifdef _OPENMP
        #pragma omp critical(Status)
        #endif
        {
            ++Nstatus;
            if ((Nstatus%COV_STATUS_PERIOD == 0) && (d->verbosity > 0))
            {
                SAFEHMPDF_NORETURN(status_update(d, start_time, Nstatus, d->n->Nphi));
            }
        }
        if (hmpdf_status) { continue; }
    }

    // subtract the one-point outer product
    double weight_sum = 1.0; // includes the zero-separation contribution here
    for (int pp=0; pp<d->n->Nphi; pp++)
    {
        weight_sum += d->n->phiweights[pp];
    }
    for (int ii=0; ii<d->n->Nsignal; ii++)
    {
        for (int jj=0; jj<d->n->Nsignal; jj++)
        {
            d->cov->Cov[ii*d->n->Nsignal+jj] -= weight_sum
                                                * d->op->PDFc[ii]
                                                * d->op->PDFc[jj];
        }
    }

    d->cov->created_cov = 1;

    ENDFCT
}//}}}

static int
prepare_cov(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(1, "prepare_cov\n");

    // run necessary code from other modules
    if (d->f->Nfilters > 0)
    {
        SAFEHMPDF(create_conj_profiles(d));
        SAFEHMPDF(create_filtered_profiles(d));
    }
    SAFEHMPDF(create_segments(d));
    SAFEHMPDF(create_op(d));

    SAFEHMPDF(create_corr(d));

    SAFEHMPDF(create_phi_indep(d));
    
    // create phi grid
    SAFEHMPDF(create_phigrid(d));
    
    // allocate the workspaces
    SAFEHMPDF(create_tp_ws(d));
    
    // create covariance matrix
    SAFEHMPDF(create_cov(d));

    // create noisy covariance matrix
    if (d->ns->noise > 0.0)
    {
        SAFEHMPDF(create_noisy_cov(d));
    }

    ENDFCT
}//}}}

int
hmpdf_get_cov(hmpdf_obj *d, int Nbins, double binedges[Nbins+1], double cov[Nbins*Nbins], int noisy)
{//{{{
    STARTFCT

    SAFEHMPDF(pdf_check_user_input(d, Nbins, binedges, noisy));

    if (d->f->pixelside < 0.0)
    {
        HMPDFERR("pixel sidelength must be set for covariance matrix.");
    }

    // perform the computation
    SAFEHMPDF(prepare_cov(d));

    // if kappa, adjust the bins
    double _binedges[Nbins+1];
    memcpy(_binedges, binedges, (Nbins+1) * sizeof(double));
    if (d->p->stype == hmpdf_kappa)
    {
        for (int ii=0; ii<=Nbins; ii++)
        {
            _binedges[ii] += d->op->signalmeanc;
        }
    }

    // perform the binning
    HMPDFPRINT(3, "\t\tbinning the covariance matrix\n");
    SAFEHMPDF(bin_2d((noisy) ? d->n->Nsignal_noisy : d->n->Nsignal,
                     (noisy) ? d->n->signalgrid_noisy : d->n->signalgrid,
                     (noisy) ? d->cov->Cov_noisy : d->cov->Cov,
                     COVINTEGR_N, Nbins, _binedges, cov, TPINTERP_TYPE));

    // compute the shot noise term
    double *temp;
    SAFEALLOC(temp, malloc(Nbins * sizeof(double)));
    SAFEHMPDF(hmpdf_get_op(d, Nbins, binedges, temp, 1, noisy));
    SAFEHMPDF(add_shotnoise_diag(Nbins, cov, temp));
    free(temp);

    // normalize properly
    SAFEHMPDF(rescale_to_fsky1(d, Nbins, cov));

    ENDFCT
}//}}}

int
_get_Nphi(hmpdf_obj *d, int *Nphi)
{//{{{
    STARTFCT
    *Nphi = d->n->Nphi;
    ENDFCT
}//}}}

int
_get_phi(hmpdf_obj *d, double *phi)
{//{{{
    STARTFCT
    memcpy(phi, d->n->phigrid, d->n->Nphi * sizeof(double));
    ENDFCT
}//}}}

int
_get_phiweights(hmpdf_obj *d, double *phiweights)
{//{{{
    STARTFCT
    memcpy(phiweights, d->n->phiweights, d->n->Nphi * sizeof(double));
    ENDFCT
}//}}}

int
_get_corr_diagn(hmpdf_obj *d, double *corr_diagn)
{//{{{
    STARTFCT
    memcpy(corr_diagn, d->cov->corr_diagn, d->n->Nphi * sizeof(double));
    ENDFCT
}//}}}

int
hmpdf_get_cov_diagnostics(hmpdf_obj *d, int *Nphi, double **phi, double **phiweights, double **corr_diagn)
{//{{{
    STARTFCT

    // perform the computation
    SAFEHMPDF(prepare_cov(d));

    if (Nphi != NULL)
    {
        SAFEHMPDF(_get_Nphi(d, Nphi));
    }
    if (phi != NULL)
    {
        SAFEALLOC(*phi, malloc(d->n->Nphi * sizeof(double)));
        SAFEHMPDF(_get_phi(d, *phi));
    }
    if (phiweights != NULL)
    {
        SAFEALLOC(*phiweights, malloc(d->n->Nphi * sizeof(double)));
        SAFEHMPDF(_get_phiweights(d, *phiweights));
    }
    if (corr_diagn != NULL)
    {
        SAFEALLOC(*corr_diagn, malloc(d->n->Nphi * sizeof(double)));
        SAFEHMPDF(_get_corr_diagn(d, *corr_diagn));
    }

    ENDFCT
}//}}}

int
hmpdf_get_cov_diagnostics1(hmpdf_obj *d, double *phi, double *phiweights, double *corr_diagn)
{//{{{
    STARTFCT

    SAFEHMPDF(prepare_cov(d));

    if (phi != NULL)
    {
        SAFEHMPDF(_get_phi(d, phi));
    }
    if (phiweights != NULL)
    {
        SAFEHMPDF(_get_phiweights(d, phiweights));
    }
    if (corr_diagn != NULL)
    {
        SAFEHMPDF(_get_corr_diagn(d, corr_diagn));
    }
    
    ENDFCT
}//}}}
