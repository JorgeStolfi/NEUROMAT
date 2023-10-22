/* Last edited on 2013-05-31 01:50:08 by stolfilocal */
/* Tests the neuromat routines to display EEG data as images. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <fftw3.h>

#include <r2.h>
#include <float_image.h>
#include <float_pnm_image.h>
#include <float_image_paint.h>
#include <frgb.h>
#include <frgb_path.h>
#include <affirm.h>
#include <jsfile.h>
#include <jspnm_image.h>
#include <jspng_image.h>
#include <jspnm.h>

#include <neuromat.h>
#include <neuromat_image.h>
#include <neuromat_filter.h>

void nse_write_image(char *prefix, char *tag, int it, float_image_t *fim, double vlo, double vhi);
  /* Writes the image {fim} to a file named "{prefix}-{tag}.png" in PNG
    format. Samples are affinely mapped from {vlo} to 0 and from {vhi]}
    to 1, then gamma-mapped and quantized. */

void nse_compute_channel_stats(int nt, int nc, double **val, double vavg[], double vvar[], double vmin[], double vmax[]);
  /* For {ic} in {0..nc-1}, stores into {vavg[ic]} the average of all
    the measurements {val[0..nt-1][ic]} of channel {ic}, and into
    {vvar[ic]} their variance about that average.  Sets both values to
    {NAN} if {nt} is zero; sets the variances to zero if {nt} is 1.
    Also saves the min and max values in {vmin[ic],vmax[ic]}. */

float_image_t *nse_colorize_image(float_image_t *fim, double vmax, int style);
  /* Returns a false-color version of image {fim}, where the range
    {[-vmax _ +vmax]} is turned into various colors with
    {frgb_path_signed(z,1,style)}. */
    
float_image_t **nse_basis(int btype, int NX, int NY, int ne, r2_t *ctr, r2_t *rad, r2_t pos[]);
  /* Computes an interpolating image basis {bas[0..ne-1]} of type
    {btype} for electrodes in the positions {pos[0..ne-1]}. Assumes that
    the electrode positions are inside the axis-aligned ellipse with
    center {ctr} and radius {rad}. */
  
void nse_write_eeg_spectra(char *prefix, char *tag, int nt, int ne, char *name[], int kfmax, double fsmp, double **pwr);
  /* Writes the power spectra {pwr[0..ne-1][0..kfmax]} of {ne} electrode
    signals to the file "{prefix}-{tag}-pwr.txt". Writes a few header
    lines and then calls {neuromat_filter_write_spectra} to write the data. */

void nse_write_eeg_signals(char *prefix, char *tag, int nt, int nc, char *name[], int ne, double fsmp, int rstep, double flo0, double flo1, double fhi1, double fhi0, double **val);
  /* Writes the signals {val[0..nt-1][0..nc-1]}, comprising {nc}
    channels sampled at {nt} times, to the file "{prefix}-{tag}.txt".
    Writes only one very {rstep} frames, beginning with frame 
    {floor(rstep/2)}.
    Writes header lines with the sampling frequency {fsmp}, the filter
    parameters {flo0,flo1,fhi0,fhi1} if not {NAN}, and the number of
    electrodes {ne}, and the channel names {name[0..nc-1]} and then
    calls {neuromat_write_eeg_signals} to write the data. */

float_image_t *nse_make_trigger_channel_plot(int NX, int NY, int nt, int nc, double **val, int ictrig, double vthr, int *bxminP, int *bxmaxP);
  /* Creates an RGB image with {NX} columns and {NY} rows containing the 
    a slider bar for {nt} data frames with {nc} channels, and paints over 
    that slide bar the tic marks corresponding to the data frames where the
    trigger channel is above {vthr}. Assumes that {val[it][ictrig]} is the 
    value of the trigger channel in data frame {it}. 
    Also returns teh range of pixel columns {*bxminP..*bxmaxP}
    that is spanned by the slider bar. */

void nse_choose_triger_dot_position(double xmin, double xmax, double ymin, double ymax, r2_t *ctr, r2_t *rad, r2_t *ctrig, double *rtrig);
  /* Chooses the center {*ctrig} and radius {*rtrig} for the dot that
    represents the trigger channel in the movie frame images. Assumes
    that the schematic plot of the head is an ellipse with center {ctr}
    and radii {rad}, within the rectangle {[xmin_xmax]×[ymin_ymax]}. The
    dot will be located outside the ellipse, near the upper left corner
    of that rectangle. */

int main(int argc, char **argv)
  {
    int nargs = 3;
    demand(argc >= nargs, "bad args");
    char *prefix = argv[1];     /* Prefix for all output filenames. */
    int nskip = atoi(argv[2]);  /* Index of first data frame to use. */
    int nt = atoi(argv[3]);     /* Number of data frames to use. */
    
    fprintf(stderr, "output files are named %s-*\n", prefix);
    fprintf(stderr, "skipping %d frames, displaying %d frames\n", nskip, nt);

    /* Image dimensions: */
    int NX = 120;
    int NY = 160;
    
    /* Area reserved for the frame cursor: */
    int NY0 = 20;
        
    /* Geometry of nominal scalp projection ellipse: */
    r2_t ctr = (r2_t){{ 0.50*NX, 0.50*(NY0 + NY) }};
    r2_t rad = (r2_t){{ 0.45*NX, 0.45*(NY - NY0) }};

    /* Get the electrode positions: */
    int ne = 20;
    r2_t *el_pos = neuromat_image_e20_schematic_2D_points(&ctr, &rad);
    
    /* Get the channel names: */
    int nc = ne + 1; /* Includes {ne} channels plus one trigger (?) channel. */
    char **ch_name = neuromat_e20_signal_names();

    /* Make the head mask: */
    float_image_t *msk = neuromat_image_schematic_head_mask(NX, NY, &ctr, &rad);
    nse_write_image(prefix, "m", 0, msk, 0.0, 1.0);
    
    /* Make the channel basis images: */
    int btype = 2;
    float_image_t **bas = nse_basis(btype, NX, NY, ne, &ctr, &rad, el_pos);
    
    /* Write the channel basis images: */
    int ie;
    for (ie = 0; ie < ne; ie++) 
      { double vmax = +1.50;      /* Expected max abs basis value. */
        frgb_t bck = frgb_Black;  /* Dot color for main channel. */
        frgb_t wht = frgb_White;  /* Dot color for other channels. */
        double drad = 0.75;       /* Dot radius. */
        int style = 2;
        float_image_t *cim = nse_colorize_image(bas[ie], vmax, style);
        neuromat_image_draw_electrodes(cim, ne, el_pos, drad, ie, bck.c, wht.c);
        nse_write_image(prefix, "b", ie, cim, 0.0, 1.0);
        float_image_free(cim);
      } 
    
    /* Assumed sampling frequency (Hz): */
    double fsmp = 600.0;
    
    /* Read the EEG data, and write it out whole: */
    double **val = neuromat_read_eeg_signals(stdin, nskip, nt, nc);
    nse_write_eeg_signals(prefix, "r", nt, nc, ch_name, ne, fsmp, 1, NAN, NAN, NAN, NAN, val);
    
    /* Compute and write the power spectra of individual electrodes: */
    int kfmax = nt/2;
    double **pwr = neuromat_filter_compute_spectra(nt, ne, val, kfmax);
    nse_write_eeg_spectra(prefix, "r", nt, ne, ch_name, kfmax, fsmp, pwr);
    for (ie = 0; ie < ne; ie++) { free(pwr[ie]); } 
    free(pwr);
    
    /* Filter the signals, and write them out subsampled: */
    double flo0 = 1.0e-15; /* Suppress the constant term... */
    double flo1 = 0.2;     /* Keep components with {T > 5sec}. */
    double fhi1 = 35.0;    /* Keep up to 35 Hz. */
    double fhi0 = 58.0;    /* Delete above 58 Hz. */
    int rstep = 6;         /* After filtering, can resample at 100 Hz. */
        
    auto double bandpass_gain(double f, double fmax);
    
    double bandpass_gain(double f, double fmax)
      { return neuromat_filter_gain_soft_bandpass(f, flo0, flo1, fhi1, fhi0, fmax); }

    neuromat_filter_apply(nt, ne, val, fsmp, bandpass_gain);
    nse_write_eeg_signals(prefix, "f", nt, nc, ch_name, ne, fsmp/rstep, rstep, flo0, flo1, fhi1, fhi0, val);
    
    /* Compute and write the power spectra of filtered electrode signals: */
    double **pwf = neuromat_filter_compute_spectra(nt, ne, val, kfmax);
    nse_write_eeg_spectra(prefix, "f", nt, ne, ch_name, kfmax, fsmp, pwf);
    for (ie = 0; ie < ne; ie++) { free(pwf[ie]); } 
    free(pwf);

    /* Compute basic statistics of filtered electrode signals: */
    double vavg[nc], vvar[nc], vmin[nc], vmax[nc];
    nse_compute_channel_stats(nt, nc, val, vavg, vvar, vmin, vmax);
    double vvar_max = 1.0e-30;  /* Maximum variance among electrode channels. */
    fprintf(stderr, "--- filtered signal statistics ---\n"); 
    int ic;
    for (ic = 0; ic < nc; ic++) 
      { fprintf(stderr, "channel %3d = %-4s", ic, ch_name[ic]); 
        fprintf(stderr, "  avg = %+10.5f  dev = %10.5f", vavg[ic], sqrt(vvar[ic])); 
        fprintf(stderr, "  min = %+10.5f  max = %10.5f\n", vmin[ic], vmax[ic]); 
        if ((ic < ne) && (vvar[ic] > vvar_max)) { vvar_max = vvar[ic]; }
      }
    double vdev = sqrt(vvar_max); /* Max standard deviation of all signals. */
    fprintf(stderr, "maximum dev = %8.5f\n", vdev); 
    fprintf(stderr, "\n");
    
    /* Build image of trigger channel, and get its column span {bxmin..bxmax}: */
    int bxmin, bxmax;
    int ictrig = nc - 1;
    double vthr = 0.5; /* Trigger chan seems to be either 0 or 100 only. */
    float_image_t *cit = nse_make_trigger_channel_plot(NX, NY0, nt, nc, val, ictrig, vthr, &bxmin, &bxmax);
    
    /* Geometry of the trigger channel dot: */
    r2_t ctrig;
    double rtrig;
    nse_choose_triger_dot_position(0, (double)NX, (double)NY0, (double)NY, &ctr, &rad, &ctrig, &rtrig);
    
    /* Paint the frame images:  */
    float_image_t *fim = float_image_new(1, NX, NY);
    
    double zmax = 2.0; /* Expected max value of interpolated potentials. */
    double vtrigmax = fmax(vmax[ictrig],-vmin[ictrig]); /* Max abs value of trigger channel. */
    int it;
    double coef[ne];
    for (it = 0; it < nt; it += rstep)
      { float_image_fill(fim, 0.0);
        /* Combine the interp basis images with the channel values as coefs. */
        /* Each channel is mean-shifted and divided by the common scale factor {vdev}. */ 
        double *valt = val[it];
        for (ie = 0; ie < ne; ie++) { coef[ie] = (valt[ie] - vavg[ie])/vdev; }
        neuromat_image_paint_potentials(ne, coef, bas, msk, 0, fim);
        
        /* Paint the trigger channel as a dot at the upperleft corner: */
        float ztrig = valt[ictrig]/vtrigmax; /* Trigger value scaled to {[-1..+1]} */
        float_image_paint_dot(fim, 0, ctrig.c[0], ctrig.c[1], rtrig, 0, TRUE, FALSE, ztrig, NAN, 4);

        /* Colorize the image yielding {cim} and draw channel positions: */
        int style = 2;
        double drad = 0.75; /* Dot radius. */
        frgb_t col = (frgb_t){{ 1.000, 0.900, 0.000 }}; /* Dot color for channels. */
        float_image_t *cim = nse_colorize_image(fim, zmax, style);
        neuromat_image_draw_electrodes(cim, ne, el_pos, drad, -1, col.c, col.c);
        float_image_assign_rectangle(cim, 0, NX-1, 0, NY0-1, cit, 0, 0);
        neuromat_image_paint_slider_dot(cim, bxmin, bxmax, 0.5*NY0, (it + 0.5), nt);

        nse_write_image(prefix, "f", nskip + it, cim, 0.0, 1.0);
        float_image_free(cim);
      }
      
    float_image_free(cit);
    float_image_free(fim);
    float_image_free(msk);
    for (ie = 0; ie < ne; ie++) { float_image_free(bas[ie]); } 
    free(bas);
    for (it = 0; it < nt; it++) { free(val[it]); } 
    free(val);
    free(ch_name);
    free(el_pos);
    
    return 0;
  }

float_image_t *nse_make_trigger_channel_plot(int NX, int NY, int nt, int nc, double **val, int ictrig, double vthr, int *bxminP, int *bxmaxP)
  {
    demand((ictrig >= 0) && (ictrig < nc), "invalid {ictrig}");
    demand(NY % 4 == 0, "trigger image height must be multiple of 4");
    demand(NY >= 8, "trigger image height too small");
    demand(NX > 3*NY/2, "trigger image aspect too small");
    
    /* Build an image with the slider bar, and get its column span {*bxminP..*bxmaxP}: */
    float_image_t *img = neuromat_image_make_slider_bar(3, NX, NY, bxminP, bxmaxP);
    
    /* Choose row span {mymin..mymax} of the trigger markers: */
    int mysize = NY - 4;
    int mymin = (NY - mysize)/2;
    int mymax = NY - 1 - mymin;
    
    /* Paint the trigger marks: */
    neuromat_image_paint_trigger_marks(img, *bxminP, *bxmaxP, mymin, mymax, nt, nc, val, ictrig, vthr);

    return img;
  }
    
void nse_choose_triger_dot_position(double xmin, double xmax, double ymin, double ymax, r2_t *ctr, r2_t *rad, r2_t *ctrig, double *rtrig)
  {
    /* Get distance {df} from upper left rect corner to ellipse edge in ellipse's coord system: */
    double hx = ctr->c[0] - xmin;
    double hy = ymax - ctr->c[1];
    double df = hypot(hx/rad->c[0], hy/rad->c[1]) - 1.0;
    /* Place dot on that line segment: */
    double drel = 0.40; /* Relative pos of center from end of diagonal to ellipse edge. */
    double dx = drel*rad->c[0]*df;
    double dy = drel*rad->c[1]*df;
    (*ctrig) = (r2_t){{ xmin + dx, ymax - dy }};
    (*rtrig) = 0.75*fmin(dx, dy);
  }

void nse_write_eeg_signals(char *prefix, char *tag, int nt, int nc, char *name[], int ne, double fsmp, int rstep, double flo0, double flo1, double fhi1, double fhi0, double **val)
  {
    char *fname = NULL;
    asprintf(&fname, "%s-%s.txt", prefix, tag);
    FILE *wr = open_write(fname, TRUE);

    fprintf(wr, "nt = %d\n", nt);
    fprintf(wr, "nc = %d\n", nc);
    fprintf(wr, "channels =");
    int ic;
    for (ic = 0; ic < nc; ic++) { fprintf(wr, " %s", name[ic]); }
    fprintf(wr, "\n");
    fprintf(wr, "ne = %d\n", ne);
    fprintf(wr, "fsmp = %14.8f\n", fsmp);
    if (! isnan(flo0)) { fprintf(wr, "flo0 = %14.8f\n", flo0); }
    if (! isnan(flo1)) { fprintf(wr, "flo1 = %14.8f\n", flo1); }
    if (! isnan(fhi1)) { fprintf(wr, "fhi1 = %14.8f\n", fhi1); }
    if (! isnan(fhi0)) { fprintf(wr, "fhi0 = %14.8f\n", fhi0); }

    neuromat_write_eeg_signals(wr, nt, nc, val, rstep/2, nt-1, rstep);
    fclose(wr);
    free(fname);
  }

void nse_write_eeg_spectra(char *prefix, char *tag, int nt, int ne, char *name[], int kfmax, double fsmp, double **pwr)
  {
    demand(2*kfmax <= nt, "{kfmax} too high");
    
    char *fname = NULL;
    asprintf(&fname, "%s-%s-pwr.txt", prefix, tag);
    FILE *wr = open_write(fname, TRUE);

    fprintf(wr, "nt = %d\n", nt);
    fprintf(wr, "ne = %d\n", ne);
    fprintf(wr, "channels =");
    int ie;
    for (ie = 0; ie < ne; ie++) { fprintf(wr, " %s", name[ie]); }
    fprintf(wr, "\n");
    fprintf(wr, "kfmax = %d\n", kfmax);
    fprintf(wr, "fsmp = %14.8f\n", fsmp);

    neuromat_filter_write_spectra(wr, nt, ne, kfmax, fsmp, pwr);
    fclose(wr);
    free(fname);
  }

float_image_t **nse_basis(int btype, int NX, int NY, int ne, r2_t *ctr, r2_t *rad, r2_t pos[])
  {
    float_image_t **bas = NULL;
    switch(btype)
      { 
        case 0:
          { bas = neuromat_image_basis_shepard(NX, NY, ne, pos);
          }
          break;
        case 1:
          { bas = neuromat_image_basis_poly(NX, NY, ne, ctr, rad, pos);
          }
          break;
        case 2:
          { double rho = 2.0;
            bas = neuromat_image_basis_radial(NX, NY, ne, ctr, rad, rho, pos);
          }
          break;
        case 3:
          { bas = neuromat_image_basis_voronoi(NX, NY, ne, ctr, rad, pos);
          }
          break;
        default:
          assert(FALSE);
      }
    return bas;
  }

void nse_write_image(char *prefix, char *tag, int it, float_image_t *fim, double vlo, double vhi)
  {
    /* Convert the float image to integer image {pim}: */
    int NC = fim->sz[0];
    int chns = NC;
    assert((chns == 1) || (chns == 3));
    pnm_sample_t maxval = 65535;
    bool_t yup = TRUE;
    bool_t verbose_cvt = FALSE;
    pnm_image_t *pim = float_image_to_pnm_image(fim, FALSE, chns, NULL, NULL, NULL, maxval, yup, verbose_cvt);
    
    /* Write {pim} to disk: */
    char *fname = NULL;
    asprintf(&fname, "%s-%s-%06d.png", prefix, tag, it);
    double gamma = 2.5;
    bool_t verbose_png = FALSE;
    jspng_image_write (fname, pim, gamma, verbose_png);
    
    /* Cleanup: */
    pnm_image_free(pim);
    free(fname);
  }

float_image_t *nse_colorize_image(float_image_t *fim, double vmax, int style)
  { 
    int NX = fim->sz[1];
    int NY = fim->sz[2];
    
    float_image_t *cim = float_image_new(3, NX, NY);

    int ix, iy;
    for (iy = 0; iy < NY; iy++)
      { for (ix = 0; ix < NX; ix++)
          { double val = float_image_get_sample(fim, 0, ix, iy); 
            if (isnan(val))
              { float_image_fill_pixel(cim, ix, iy, 0.0); }
            else
              { double z = val/vmax;
                if (z < -1) { z = -1; }
                if (z > +1) { z = +1; }
                frgb_t rgb = frgb_path_map_signed(z, 1, style);
                float_image_set_pixel(cim, ix, iy, rgb.c);
              }
          }
      }

    return cim;
  }

void nse_compute_channel_stats(int nt, int nc, double **val, double vavg[], double vvar[], double vmin[], double vmax[])
  {
    int it, ic;
    for (ic = 0; ic < nc; ic++)
      { double sumv = 0;
        double minv = +INF;
        double maxv = -INF;
        for (it = -0; it < nt; it++) 
          { double v = val[it][ic];
            sumv += v; 
            if (v < minv) { minv = v; }
            if (v > maxv) { maxv = v; }
          }
        vmin[ic] = minv;
        vmax[ic] = maxv;
        double avg = sumv/nt;
        vavg[ic] = avg;
        double sumd2 = 0;
        for (it = -0; it < nt; it++) 
          { double d = val[it][ic] - avg; 
            sumd2 += d*d;
          }
        double var = sumd2/nt;
        vvar[ic] = var;
      }
  }
