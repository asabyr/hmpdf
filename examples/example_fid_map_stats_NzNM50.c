/*! [compile] */
/* export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../ */
/* gcc --std=gnu99 -I$TACC_GSL_INC -I/work2/07833/tg871330/stampede2/software/hmpdf/include -o example_fid_map_stats_NzNM50 example_fid_map_stats_NzNM50.c -L/work2/07833/tg871330/stampede2/software/hmpdf/ -lhmpdf */
/*! [compile] */
#include <stdio.h>
#include "utils.h"
#include "hmpdf.h"

int example_map(void)
{
    /* get a new hmpdf_obj */
    hmpdf_obj *d = hmpdf_new();
    if (!(d))
        return -1;

    /* initialize with default settings */
    if (hmpdf_init(d, "fid_cosmo.ini", hmpdf_tsz,
                hmpdf_z_min, 0.005,   
		hmpdf_map_fsky, 0.01,
               hmpdf_pixel_side, 0.5,
               hmpdf_verbosity, 5,
 		hmpdf_N_z, 50,
		hmpdf_N_M, 50))
        return -1;

    /* get the one-point PDF */
    double *map;
    long Nside;
    
    if (hmpdf_get_map(d, &map, &Nside, 1))
        return -1;

    FILE *fp = fopen("example_fid_map_stats.bin", "w");
    fwrite(map, sizeof(double), Nside*Nside, fp);
    fclose(fp);

    free(map);

    /* free memory associated with the hmpdf_obj */
    if (hmpdf_delete(d))
        return -1;
    return 0;
}

int main(void)
{
    if (example_map())
    {
        fprintf(stderr, "failed\n");
        return -1;
    }
    else
    {
        return 0;
    }
}

