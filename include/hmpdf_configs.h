/*! \file hmpdf_configs.h */
#ifndef HMPDF_CONFIGS_H
#define HMPDF_CONFIGS_H

/*! Halo mass definitions.
 *  
 *  Presently only used to specify in terms of which
 *  definition the radial cut-off of the halo profiles
 *  is given.
 */
typedef enum
{
    // NOTE do not change ordering here!
    hmpdf_mdef_c, /*!< M200c (200x critial density) */
    hmpdf_mdef_v, /*!< Mvir (virial mass according to Bryan+Norman 1998) */
    hmpdf_mdef_m, /*!< M200m (200x mean matter density) */
} hmpdf_mdef_e;

/*! Signal types.
 *
 *  To specify which cosmological field the PDF should be computed for.
 */
typedef enum
{
    hmpdf_kappa, /*!< weak lensing convergence */
    hmpdf_tsz,   /*!< tSZ effect (Compton-y) */
} hmpdf_signaltype_e;

/*! Fixed point integration modes.
 *
 *  The integrals over halo mass and redshift are performed using fixed point quadratures,
 *  with modes taken from this enum.
 *  See the
 *  <a href="https://www.gnu.org/software/gsl/doc/html/integration.html#fixed-point-quadratures">
 *  GSL documentation</a>
 *  for details.
 *
 *  \attention Some of these integration modes are for infinite intervals and do not make sense
 *             for our application.
 */
typedef enum
{
    hmpdf_legendre, /*!< .*/
    hmpdf_chebyshev, /*!< .*/
    hmpdf_gegenbauer, /*!< .*/
    hmpdf_jacobi, /*!< .*/
    hmpdf_laguerre, /*!< .*/
    hmpdf_hermite, /*!< .*/
    hmpdf_exponential, /*!< .*/
    hmpdf_rational, /*!< .*/
    hmpdf_chebyshev2, /*!< .*/
} hmpdf_integr_mode_e;

/*! Ordering of the Arico+2020 parameters for the BCM */
typedef enum
{
    hmpdf_Arico20_M_c, /*!<.*/
    hmpdf_Arico20_M_1_z0_cen, /*!<.*/
    hmpdf_Arico20_eta, /*!<.*/
    hmpdf_Arico20_beta, /*!<.*/
#ifdef ARICO20
    hmpdf_Arico20_theta_inn, /*!<.*/
    hmpdf_Arico20_theta_out, /*!<.*/
    hmpdf_Arico20_M_inn, /*!<.*/
    hmpdf_Arico20_M_r, /*!<.*/
#endif
    hmpdf_Arico20_Nparams, /*!<. Internal use only. */
} hmpdf_Arico20_params_e;

/*! Options to hmpdf_init().
 *
 *  The variadic argument list in hmpdf_init() can be used to pass non-default options.
 *  [syntax is explained in the documentation for hmpdf_init()].
 *
 *  There is a large number of options, many of which the typical user will not need to use.
 *  Here is a brief synopsis:
 *
 *  Frequently used options:
 *      + pixelization: #hmpdf_pixel_side
 *      + ell space filters: #hmpdf_tophat_radius, #hmpdf_gaussian_fwhm,
 *                           #hmpdf_custom_ell_filter (and #hmpdf_custom_ell_filter_params)
 *      + Gaussian noise: #hmpdf_noise_pwr (and #hmpdf_noise_pwr_params)
 *      + multithreading: #hmpdf_N_threads
 *  
 *  Less frequently used options:
 *      + verbosity: #hmpdf_verbosity
 *      + behaviour when unusual states are encountered: #hmpdf_warn_is_err
 *      + halo model fit parameters: #hmpdf_Duffy08_conc_params,
 *                                   #hmpdf_Tinker10_hmf_params,
 *                                   #hmpdf_Battaglia12_tsz_params
 *      + k space filter: #hmpdf_custom_k_filter (and #hmpdf_custom_k_filter_params)
 *      + PDF internal sampling points: #hmpdf_N_signal, #hmpdf_signal_min, #hmpdf_signal_max
 *      + settings for simplified simulations: #hmpdf_map_fsky,
 *                                             #hmpdf_map_pixelgrid,
 *                                             #hmpdf_map_poisson
 *  
 *  Integration grids:
 *      + redshift integration: #hmpdf_N_z, #hmpdf_z_min, #hmpdf_z_max,
 *                              #hmpdf_zintegr_type, #hmpdf_zintegr_alpha, #hmpdf_zintegr_beta
 *      + halo mass integration: #hmpdf_N_M, #hmpdf_M_min, #hmpdf_M_max,
 *                               #hmpdf_Mintegr_type, #hmpdf_Mintegr_alpha, #hmpdf_Mintegr_beta
 *      + halo profile angular integration: #hmpdf_N_theta, #hmpdf_rout_scale, #hmpdf_rout_rdef
 *
 *  Covariance matrix calculation:
 *      + useful to improve numerical stability: #hmpdf_N_phi
 *      + integration/summation grid: #hmpdf_phi_max, #hmpdf_pixelexact_max, #hmpdf_phi_jitter,
 *                                    #hmpdf_phi_pwr
 */
typedef enum
{
    hmpdf_N_threads, /*!< number of threads to use in multithreaded parts of the code.
                    *   \par
                    *   Type: int. Default: 1.
                    *   \remark only applicable if code compiled with OpenMP.
                    *   \remark this does *not* control the multithreading behaviour of CLASS
                    */
    hmpdf_verbosity, /*!< larger values yield more detailed print output.
                      *   If 0 (default), only error messages are printed.
                      *   \par
                      *   Type: int. Default: 0.
                      */
    hmpdf_warn_is_err, /*!< Controls the behaviour when a warning is encountered.
                        *   Positive values will treat warnings as error and return,
                        *   zero value will continue execution but print a message to
                        *   stderr,
                        *   and negative values will ignore the warning.
                        *   \par
                        *   Type: int. Default: 1.
                        */
    hmpdf_class_pre, /*!< optional CLASS precision file.
                      *   \par
                      *   Type: char *. Default: None (use CLASS's default precision settings).
                      */
    hmpdf_N_z, /*!< number of sample points in redshift integration.
                *   \par
                *   Type: int. Default: 65.
                *   \remark the default value is conservative.
                *           By playing with the #hmpdf_zintegr_type
                *           (as well as #hmpdf_zintegr_alpha and #hmpdf_zintegr_beta)
                *           you can gain some speed here.
                */
    hmpdf_z_min, /*!< minimum redshift.
                  *   \par
                  *   Type: double. Default: 0.
                  *   \remark non-zero values do not make much sense
                  */
    hmpdf_z_max, /*!< maximum redshift.
                  *   \par
                  *   Type: double. Default: 6 for tSZ, source redshift for weak lensing.
                  */
    hmpdf_dndz, /*!< source redshift distribution for #hmpdf_kappa.
                 *   If set to null (default), a delta distribution
                 *   is assumed at source redshift. Otherwise, source redshift should be chosen
                 *   such that the source distribution is essentially zero for larger redshifts.
                 *   \par
                 *   Type: #hmpdf_dndz_f. Default: None.
                 *   \remark it is not assumed that the distribution is normalized, the code
                 *           will do this automatically.
                 */
    hmpdf_dndz_params, /*!< additional parameters to pass to hmpdf_dndz.
                        *   \par
                        *   Type: void *. Default: None.
                        */
    hmpdf_N_M, /*!< number of sample points in halo mass integration.
                *   \par
                *   Type: int. Default: 65.
                *   \remark the default value is conservative.
                *           By playing with the #hmpdf_Mintegr_type
                *           (as well as #hmpdf_Mintegr_alpha and #hmpdf_Mintegr_beta)
                *           you can gain some speed here.
                */
    hmpdf_M_min, /*!< minimum halo mass, in M200m definition and Msun units.
                  *   \par
                  *   Type: double. Default: 1e11.
                  */
    hmpdf_M_max, /*!< maximum halo mass, in M200m definition and Msun units.
                  *   \par
                  *   Type: double. Default: 1e16.
                  */
    hmpdf_mass_z_fix, // whether to place a mass and redshift cut, which defines which halos have a fixed tSZ profile (Battaglia+12) 
                      // Type: int. Default: 0
    hmpdf_min_mass_fix, // Above this mass and below hmpdf_max_z_fix redshift, the profiles are fixed (in M200c)
                        //Type: double. Default: 8*10^14
    hmpdf_max_z_fix, //Below this redshift and above M200c hmpdf_min_mass_fix mass, the profiles are fixed
                    //Type: double. Default: 0.5
    hmpdf_N_signal, /*!< number of points on which the one-point PDF is sampled internally
                     *   \par
                     *   Type: long. Default: 1024.
                     *   \remark two-point PDF will be sampled on the square of this.
                     *   \remark Should be chosen as a power of 2 for optimal speed.
                     *   \warning Setting this to large values can trigger numerical instability
                     *            in small-phi two-point PDF, and covariance matrix computations.
                     *   \attention This is the only instance of an expected long.
                     */
    hmpdf_signal_min, /*!< minimum signal value at which PDF is sampled internally.
                       *   \par
                       *   Type: double. Default: 0.
                       *   \remark The PDF should be well converged to zero for smaller signal values.
                       *   \remark If your signal is only positive
                       *           (this is quite often true, as long as you are not applying any
                       *            special ell- or k-space filters)
                       *           you should set this to zero.
                       *   \note   This is in the internal units,
                       *           in which negative WL convergence regions do not exist.
                       *           Read the paper for a discussion.
                       */
    hmpdf_signal_max, /*!< maximum signal value at which PDF is sampled internally.
                       *   \par
                       *   Type: double. Default: 2e-4 (tSZ), 1 (weak lensing).
                       *   \remark Should be chosen such that PDF is well converged to zero to avoid ringing.
                       *   \remark Should be at least 2x larger than maximum signal value you are interested in.
                       *   \remark For weak lensing, you should change this option depending on source redshift.
                       *   \note   This is in the internal units,
                       *           in which negative WL convergence regions do not exist.
                       *           Read the paper for a discussion.
                       */
    hmpdf_N_theta, /*!< number of points on which the signal profiles are sampled.
                    *   \par
                    *   Type: int. Default: 500.
                    *   \remark reasonable precision can be reached with as few as 100 sample points.
                    *           Runtime of hmpdf_init() is in many cases quite dominated by this setting.
                    */
    hmpdf_rout_scale, /*!< radial cut-off.
                       *   \par
                       *   Type: double. Default: 2.
                       */
    hmpdf_rout_rdef, /*!< mass/radius definition in terms of which the radial cut-off is specified.
                      *   \par
                      *   Type: #hmpdf_mdef_e. Default: #hmpdf_mdef_v (virial radius).
                      */
    hmpdf_pixel_side, /*!< pixel sidelength, in arcmin.
                       *   \par
                       *   Type: double. Default: None.
                       *   \remark required setting for covariance matrix calculation.
                       *   \remark negative values have no effect (but will not trigger a warning).
                       *   \warning PDFs computed with more than 5 arcmin in this setting should be treated
                       *            with caution.
                       */
    hmpdf_tophat_radius, /*!< includes the effect of smoothing the map with a tophat of given radius (in arcmin).
                          *   \par
                          *   Type: double. Default: None.
                          *   \remark negative values have no effect (but will not trigger a warning).
                          *   \warning PDFs computed with more than 5 arcmin in this setting should be treated
                          *            with caution.
                          */
    hmpdf_gaussian_fwhm, /*!< includes the effect of smoothing the map with a Gaussian of given FWHM (in arcmin).
                          *   \par
                          *   Type: double. Default: None.
                          *   \remark negative values have no effect (but will not trigger a warning).
                          *   \warning PDFs computed with more than 5 arcmin in this setting should be treated
                          *            with caution.
                          */
    hmpdf_custom_ell_filter, /*!< pass a user-defined ell-space filter
                              *   (for example, to include effect of a Wiener filter).
                              *   Signature of this function pointer has to conform
                              *   to the typedef #hmpdf_ell_filter_f.
                              *   \par
                              *   Type: #hmpdf_ell_filter_f. Default: None.
                              */
    hmpdf_custom_ell_filter_params, /*!< pass parameters to the ell-space filter (as its last argument).
                                     *   \par
                                     *   Type: void *. Default: None.
                                     *   \attention not supported in the python wrapper.
                                     */
    hmpdf_custom_k_filter, /*!< pass a user-defined k-space filter, with possible redshift-dependence
                            *   (for example, to emulate small-scale simulation resolution issues).
                            *   Signature of this function pointer has to conform
                            *   to the typedef #hmpdf_k_filter_f.
                            *   \par
                            *   Type: #hmpdf_k_filter_f. Default: None.
                            */
    hmpdf_custom_k_filter_params, /*!< pass parameters to the k-space filter (as its last argument).
                                   *   \par
                                   *   Type: void *. Default: None.
                                   *   \attention not supported in the python wrapper.
                                   */
    hmpdf_massfunc_corr, /*!< Correction prefactor for the Tinker mass function.
                          *   \par
                          *   Type: #hmpdf_massfunc_corr_f. Default: None.
                          */
    hmpdf_massfunc_corr_params, /*!< Additional parameters to pass to the above function.
                                 *   \par
                                 *   Type: void *. Default: None.
                                 */
    hmpdf_mass_resc, /*!< Rescaling function for halo masses going into the profiles
                      *   (masses going into mass function and bias are not affected).
                      *   \par
                      *   Type: #hmpdf_mass_resc_f. Default: None.
                      */
    hmpdf_mass_resc_params, /*!< Additional parameters to pass to the above function.
                             *   \par
                             *   Type: void *. Default: None.
                             */
    hmpdf_conc_resc, /*!< Rescaling function for concentration.
                      *   \par
                      *   Type: #hmpdf_conc_resc_f. Default: None.
                      */
    hmpdf_conc_resc_params, /*!< Additional parameters to pass to the above function.
                             *   \par
                             *   Type: void *. Default: None.
                             */
    hmpdf_mass_cuts, /*!< Redshift-dependent upper mass limits.
                      *   \par
                      *   Type: #hmpdf_mass_cuts_f. Default: None.
                      */
    hmpdf_mass_cuts_params, /*!< Additional parameters to pass to the above function.
                             *   \par
                             *   Type: void *. Default: None.
                             */
    hmpdf_bias_resc, /*!< Rescaling function for bias.
                      *   \par
                      *   Type: #hmpdf_bias_resc_f. Default: None.
                      */
    hmpdf_bias_resc_params, /*!< Additional parameters to pass to the above function.
                             *   \par
                             *   Type: void *. Default: None.
                             */
    hmpdf_Arico20_Nz, /*!< The number of parameter samples (in redshift) for the BCM.
                       *   The code interpolates in redshift between the ones that are given.
                       *   \par
                       *   Type: int. Default: 1.
                       */
    hmpdf_Arico20_z, /*!< The redshifts where the parameter values are given.
                      *   Must be in increasing order.
                      *   Does not need to be passed if Nz==1.
                      *   \par
                      *   Type: double *. Default: None.
                      */
    hmpdf_Arico20_params, /*!< The parameters for the Arico+2020 BCM.
                           *   If passed, the convergence profiles will be computed using the BCM.
                           *   Only applicable if computing the kappa PDF.
                           *   The ordering is specified by #hmpdf_Arico20_params_e.
                           *   Each redshift corresponds to a contiguous block in the array.
                           *   \par
                           *   Type: double *. Default: None.
                           */
    hmpdf_profiles_N, /*!< Number of profiles to write to file.
                       *   This is a hack for that only works with the BCM.
                       *   Will print the component density profiles to disk.
                       *   \par
                       *   Type: int. Default: 0.
                       */
    hmpdf_profiles_fnames, /*!< The file names the profiles should be printed to.
                            *   \par
                            *   Type: char **. Default: None.
                            */
    hmpdf_profiles_where, /*!< The redshifts and masses for which to evaluate profiles.
                           *   In the order [z0, M0, z1, M1, z2, M2, ... ]
                           *   Masses in units of M200m, Msun/h.
                           *   \par
                           *   Type: double *. Default: None.
                           */
    hmpdf_profiles_Nr, /*!< Number of radii where to evaluate profiles.
                        *   \par
                        *   Type: int. Default: 0.
                        */
    hmpdf_profiles_r, /*!< The radii for which to evaluate profiles.
                       *   In units of R200c.
                       *   \par
                       *   Type: double *. Default: None.
                       */
    hmpdf_tot_profiles_N, /*!< The following (all tot_profiles) are for the NFW profiles,
                               everything analogous to the above.
                           */
    hmpdf_tot_profiles_fnames, /*!< see above */
    hmpdf_tot_profiles_where, /*!< see above */
    hmpdf_tot_profiles_Nr, /*!< see above */
    hmpdf_tot_profiles_r, /*!< see above */
    hmpdf_DM_conc_params, /*!< The concentration parameters for the DM component.
                           *   If passed, also those for the baryonic component must be given.
                           *   Same format as #hmpdf_Duffy08_conc_params.
                           *   The parameterization is in terms of the hydro mass.
                           *   \par
                           *   Type: double *. Default: None.
                           */
    hmpdf_bar_conc_params, /*!< The concentration parameters for the baryonic component.
                            *   If passed, also those for the DM component must be given.
                            *   Same format as #hmpdf_Duffy08_conc_params.
                            *   The parameterization is in terms of the hydro mass.
                            *   \par
                            *   Type: double *. Default: None.
                            */
    hmpdf_N_phi, /*!< Number of pixel-separation sample points in covariance matrix calculation.
                  *   \par
                  *   Type: int. Default: 1000.
                  *   \remark as long as no numerical instability at small pixel separations is encountered,
                  *           the default value is conservative.
                  *           In that case, if you want to decrease run-time, it is recommended you decrease
                  *           this value to a few hundred, and decrease #hmpdf_pixelexact_max as well.
                  *   \remark if, however, numerical instability is encountered,
                  *           results can get more accurate by increasing this value.
                  */
    hmpdf_phi_max, /*!< maximum pixel separation in covariance matrix calculation (in arcmin).
                    *   \par
                    *   Type: double. Default: 150.
                    *   \remark should be set to approximately the largest halo radius on the sky.
                    */
    hmpdf_pixelexact_max, /*!< up to which pixel separation (in units of pixel sidelength) the summation
                           *   is performed exactly (integration afterwards).
                           *   \par
                           *   Type: int. Default: 20.
                           *   \remark the default value is conservative, can be decreased to 10 without
                           *           losing much accuracy.
                           */
    hmpdf_phi_jitter, /*!< technical (for covariance matrix).
                       *   For numerical stability, at small pixel separations an average is used,
                       *   around the desired pixel separation.
                       *   This option sets how wide (in phi) the averaging is done.
                       *   \par
                       *   Type: double. Default: 0.02.
                       */
    hmpdf_phi_pwr, /*!< technical (for covariance matrix).
                    *   Increasing this number increases the sample points at small pixel separation.
                    *   \par
                    *   Type: double. Default: 2.
                    *   \remark If you encounter numerical instability (e.g., wild covariance matrix entries),
                    *           increasing this number can potentially be helpful.
                    */
    hmpdf_zintegr_type, /*!< Fixed point integration mode for redshift integration.
                         *   \par
                         *   Type: #hmpdf_integr_mode_e. Default: #hmpdf_legendre.
                         */
    hmpdf_zintegr_alpha, /*!< Fixed point integration alpha for redshift integration
                          *   \par
                          *   Type: double. Default: 0.
                          */
    hmpdf_zintegr_beta, /*!< Fixed point integration beta for redshift integration.
                         *   \par
                         *   Type: double. Default: 0.
                         */
    hmpdf_Mintegr_type, /*!< Fixed point integration mode for halo mass integration.
                         *   \par
                         *   Type: #hmpdf_integr_mode_e. Default: #hmpdf_legendre.
                         */
    hmpdf_Mintegr_alpha, /*!< Fixed point integration alpha for halo mass integration.
                          *   \par
                          *   Type: double. Default: 0.
                          */
    hmpdf_Mintegr_beta, /*!< Fixed point integration beta for halo mass integration.
                         *   \par
                         *   Type: double. Default: 0.
                         */
    hmpdf_Duffy08_conc_params, /*!< Fit parameters in the
                                *   <a href="https://arxiv.org/abs/0804.2486">Duffy+2008</a>
                                *   concentration model.
                                *   Still relevant even if #hmpdf_DM_conc_params and #hmpdf_bar_conc_params
                                *   are passed (then it is used for mass conversions).
                                *   \par
                                *   Type: double *. Default: see src/configs.c.
                                *   \remark in the python wrapper, pass a 1d numpy array
                                */
    hmpdf_Tinker10_hmf_params, /*!< Fit parameters in the
                                *   <a href="https://arxiv.org/abs/1001.3162">Tinker+2010</a>
                                *   halo mass function.
                                *   \par
                                *   Type: double[10]. Default: see src/configs.c.
                                *   \remark in the python wrapper, pass a 1d numpy array
                                */
    hmpdf_Battaglia12_tsz_params, /*!< Fit parameters in the
                                   *   <a href="https://arxiv.org/abs/1109.3711">Battaglia+2012</a>
                                   *   pressure profile model.
                                   *   \par
                                   *   Type: double[15]. Default: see src/configs.c.
                                   *   \remark in the python wrapper, pass a 1d numpy array
                                   */
    hmpdf_noise_pwr, /*!< Option to add pixel-wise Gaussian noise with this power spectrum.
                      *   \par
                      *   Type: #hmpdf_noise_pwr_f. Default: None.
                      */
    hmpdf_noise_pwr_params, /*!< To pass custom parameters to the function #hmpdf_noise_pwr.
                             *   \par
                             *   Type: void *. Default: None.
                             */
    hmpdf_map_fsky, /*!< If you want to use the stochastic map-mapmaking algorithm
                     *   (simplified simulations), this is a required setting.
                     *   It is the sky fraction spanned by the map.
                     *   \par
                     *   Type: double. Default: None.
                     */
    hmpdf_map_pixelgrid, /*!< Increasing this value yields a more accurate averaging of the
                          *   signal profiles in each pixel.
                          *   \par
                          *   Type: int. Default: 3.
                          */
    hmpdf_map_poisson, /*!< If set to zero, the number of halos will not be drawn from the Poisson
                        *   distribution but rather from a much more concentrated one.
                        *   In that case, the mean (i.e. the averaged one-point PDF) is correctly
                        *   and much faster approached, while the scatter of individual histograms
                        *   will not be correct.
                        *   \par
                        *   Type: int. Default: 1.
                        */
    hmpdf_map_seed, /*!< If this option is set, the resulting simplified simulations (maps) will
                     *   have this random seed.
                     *   This means that, given equal map dimensions, the maps will have all halos
                     *   in identical positions.
                     *   \par
                     *   Type: int. Default: None.
                     */
    hmpdf_end_configs, /*!< required last argument in hmpdf_init_fct(), the convenience macro
                        *   hmpdf_init() takes care of that.
                        */
    // keep this last
} hmpdf_configs_e;

/*! Function pointer typedef for user-defined ell-space filter.
 *  Passed to hmpdf_init() as #hmpdf_custom_ell_filter.
 *  \param ell      angular wavenumber
 *  \param p        pointer that allows user to pass other parameters.
 *                  Passed to hmpdf_init() as #hmpdf_custom_ell_filter_params.
 *  \return W(ell)  window function at ell
 */
typedef double (*hmpdf_ell_filter_f)(double,
                                     void *);

/*! Function pointer typedef for user-defined k-space filter.
 *  Passed to hmpdf_init() as #hmpdf_custom_k_filter.
 *  \param k        comoving wavenumber in 1/Mpc
 *  \param z        redshift
 *  \param p        pointer that allows user to pass other parameters.
 *                  Passed to hmpdf_init() as #hmpdf_custom_k_filter_params.
 *  \return W(k,z)  window function at k and z
 */
typedef double (*hmpdf_k_filter_f)(double,
                                   double,
                                   void *);

/*! Function pointer typedef for user-defined mass function rescaling.
 *  Passed to hmpdf_init() as #hmpdf_massfunc_corr.
 *  \param z        redshift
 *  \param M        halo mass (M200m, Msun/h)
 *  \param p        pointer for additional parameters
 *  \return         hmf_new / hmf_Tinker at z and M
 */
typedef double (*hmpdf_massfunc_corr_f)(double,
                                        double,
                                        void *);

/*! Function pointer typedef for user-defined mass rescaling.
 *  Passed to hmpdf_init() as #hmpdf_mass_resc.
 *  \param z        redshift
 *  \param M        halo mass (M200m, Msun/h)
 *  \param p        pointer for additional parameters
 *  \return         M_new / M_old
 */
typedef double (*hmpdf_mass_resc_f)(double,
                                    double,
                                    void *);

/*! Function pointer typedef for user-defined concentration rescaling.
 *  Passed to hmpdf_init() as #hmpdf_conc_resc.
 *  \param z        redshift
 *  \param M        halo mass (M200m, Msun/h)
 *  \param p        pointer for additional parameters
 *  \return         c_new / c_old
 */
typedef double (*hmpdf_conc_resc_f)(double,
                                    double,
                                    void *);

/*! Function pointer typedef to pass redshift-dependent upper limit in mass.
 *  Passed to hmpdf_init() as #hmpdf_mass_cuts.
 *  \param z        redshift
 *  \param p        pointer for additional parameters
 *  \return         upper mass limit (M200m, Msun/h)
 */
typedef double (*hmpdf_mass_cuts_f)(double,
                                    void *);

/*! Function pointer typedef for user-defined noise power spectrum.
 *  Passed to hmpdf_init() as #hmpdf_noise_pwr.
 *  \param ell       angular wavenumber
 *  \param p         pointer that allows user to pass other parameters.
 *                   Passed to hmpdf_init() as #hmpdf_noise_pwr_params.
 *  \return N(ell)   noise power at ell.
 */
typedef double (*hmpdf_noise_pwr_f)(double,
                                    void *);

/*! Function pointer typedef for user-defined bias rescaling.
 *  Passed to hmpdf_init() as #hmpdf_bias_resc.
 *  \param z        redshift
 *  \param M        halo mass (M200m, Msun/h)
 *  \param p        pointer for additional parameters
 *  \return         rescaling of the halo bias
 */
typedef double (*hmpdf_bias_resc_f)(double,
                                    double,
                                    void *);

/*! Function pointer typedef for user-defined source redshift distribution.
 *  Passed to hmpdf_iniit() as #hmpdf_dndz.
 *  \param z        redshift
 *  \param p        pointer for additional parameters
 *  \return         dn/dz, the source redshift distribution
 */
typedef double (*hmpdf_dndz_f)(double,
                               void *);

#endif
