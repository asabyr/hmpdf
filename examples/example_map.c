/*! [compile] */
/* export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../ */
/* gcc --std=gnu99 -I../include -o example_map example_map.c -L.. -lhmpdf */
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
    if (hmpdf_init(d, "example.ini", hmpdf_tsz,
                   hmpdf_map_fsky, 0.1,
               hmpdf_pixel_side, 1.0,
               hmpdf_N_M, 10, hmpdf_N_z, 10,
               hmpdf_verbosity, 5,
		hmpdf_M_min, 1.0e13))
        return -1;

    /* get the one-point PDF */
    double *map;
    long Nside;
    
    if (hmpdf_get_map(d, &map, &Nside, 1))
        return -1;

    FILE *fp = fopen("example_map.bin", "w");
    fwrite(map, sizeof(double), Nside*Nside, fp);
    // map=np.fromfile('example_map.bin')
    // Nside=int(np.sqrt(len(map)))
    // map=map.reshape(Nside,Nside)
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

