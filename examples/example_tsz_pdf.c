/*! [compile] */
/*! gcc --std=gnu99 -I$TACC_GSL_INC -I/scratch/07833/tg871330/software_scratch/hmpdf/include -o example_tsz_pdf example_tsz_pdf.c -L/scratch/07833/tg871330/software_scratch/hmpdf/ -lhmpdf */
/*! [compile] */
#include <stdio.h>

#include "utils.h"
#include "hmpdf.h"

/*! [example_tsz_onepoint] */
int example_tsz_onepoint(void)
{
    /* construct some binedges */
    int Nbins = 100; double tszmin = 0.0; double tszmax = 4e-5;
    double binedges[Nbins+1];
    for (int ii=0; ii<=Nbins; ii++)
        binedges[ii] = tszmin + (double)(ii)*(tszmax-tszmin)/(double)(Nbins);

    /* get a new hmpdf_obj */
    hmpdf_obj *d = hmpdf_new();
    if (!(d))
        return -1;

    /* initialize with default settings */
    if (hmpdf_init(d, "example_orig.ini", hmpdf_tsz
                   /* source redshift */))
        return -1;

    /* get the one-point PDF */
    double op[Nbins];
    if (hmpdf_get_op(d, Nbins, binedges, op,
                     1/* include two-halo term */,
                     0/* don't include noise */))

        return -1;
    
      
    FILE *fp = fopen("example_tsz_pdf.bin","w");
    fwrite(binedges,sizeof(double),Nbins+1,fp);
    fwrite(op,sizeof(double),Nbins,fp);
    fclose(fp);

    /* free memory associated with the hmpdf_obj */
    if (hmpdf_delete(d))
        return -1;

    /* do something with the one-point PDF ... */

    return 0;
}
/*! [example_tsz_onepoint] */

int main(void)
{
    if (example_tsz_onepoint())
    {
        fprintf(stderr, "failed\n");
        return -1;
    }
    else
    {
        return 0;
    }
}
