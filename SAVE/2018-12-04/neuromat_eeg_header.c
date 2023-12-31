/* See {neuromat_eeg_header.h}. */
/* Last edited on 2023-10-22 02:55:02 by stolfi */
  
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <assert.h>

#include <fget.h>
#include <affirm.h>
#include <jsmath.h>
#include <jsstring.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_source.h>

#include <neuromat_eeg_header.h>

neuromat_eeg_header_t *neuromat_eeg_header_new(void)
  {
    neuromat_eeg_header_t *h = notnull(malloc(sizeof(neuromat_eeg_header_t)), "no mem");
    (*h) = (neuromat_eeg_header_t)
      { .nt = INT_MIN,
        .nc = INT_MIN,
        .chname = NULL,
        .ne = INT_MIN,
        .type = NULL,
        .component = NULL,
        .kfmax = INT_MIN,
        .tdeg = INT_MIN,
        .tkeep = INT_MIN,
        .fsmp = NAN,
        .flo0 = NAN,
        .flo1 = NAN,
        .fhi1 = NAN,
        .fhi0 = NAN,
        .finvert = INT_MIN,
        .orig = neuromat_eeg_source_new()
      };
    return h;
  }

void neuromat_eeg_header_free(neuromat_eeg_header_t *h)
  {
    if (h->chname != NULL)
      { demand(h->nc > 0, "cannot have {chname} without {nc}");
        int ic;
        for (ic = 0; ic < h->nc; ic++) { free(h->chname[ic]); }
        free(h->chname);
      }
    neuromat_eeg_source_free(h->orig);
  }

void neuromat_eeg_header_write(FILE *wr, neuromat_eeg_header_t *h)
  {
    if (h->nt != INT_MIN) { fprintf(wr, "nt = %d\n",h->nt); }   
    if (h->nc != INT_MIN) { fprintf(wr, "nc = %d\n", h->nc); }   
    if (h->chname != NULL) 
      { demand(h->nc > 0, "cannot have {channels} without {nc}");
        fprintf(wr, "channels =");
        int ic;
        for (ic = 0; ic < h->nc; ic++) { fprintf(wr, " %s", h->chname[ic]); }
        fprintf(wr, "\n");
      }   
    if (h->kfmax != INT_MIN) { fprintf(wr, "kfmax = %d\n", h->kfmax); }   
    if (h->ne != INT_MIN) { fprintf(wr,   "ne = %d\n", h->ne); } 
    if (h->type != NULL) 
      { demand(strlen(h->type) > 0, "empty run type");
        fprintf(wr, "type = %s\n", h->type);
      }   
    if (h->component != NULL) 
      { demand(strlen(h->component) > 0, "empty component");
        fprintf(wr, "component = %s\n", h->component);
      }   
    if (! isnan(h->fsmp)) { fprintf(wr, "fsmp = %.10f\n", h->fsmp); }   
    if (h->tdeg >= 0) 
      { demand((h->tkeep == 0) || (h->tkeep == 1), "inconsistent trend preservation flag");
        fprintf(wr, "trend = %d %d\n", h->tdeg, h->tkeep); }   
    if ((! isnan(h->flo0)) || (! isnan(h->flo1)) || (! isnan(h->fhi1)) || (! isnan(h->fhi0)) || (h->finvert >= 0))
      { demand((! isnan(h->flo0)) && (! isnan(h->flo1)) && (! isnan(h->fhi1)) && (! isnan(h->fhi0)), "inconsistent band limits");
        demand((h->finvert == 0) || (h->finvert == 1), "inconsistent filter inversion flag");
        fprintf(wr, "band = %.15g %.15g %.15g %.15g %d\n", h->flo0, h->flo1, h->fhi1, h->fhi0, h->finvert);
      }   
    neuromat_eeg_source_write(wr, "orig.", h->orig);
    fflush(wr);
  }

neuromat_eeg_header_t *neuromat_eeg_header_read(FILE *rd, int neDef, double fsmpDef, int *nlP)
  {
    neuromat_eeg_header_t *h = neuromat_eeg_header_new();
    
    auto void skip_line(void);
      /* Consumes the remaining characters up to and including the EOL. 
        Bombs out if EOF before EOL. */

    int nl = (nlP == NULL ? 0 : (*nlP)); /* File line number (including comments and headers), from 1. */
    int nh = 0; /* Header lines read (excluding comments and headers). */
    while (TRUE) 
      { /* Try to read one more line: */
        int r = fgetc(rd);
        if (r == EOF)
          { break; } 
        else if (r == '#')
          { /* Comment line: */
            skip_line();
          }
        else if (((r >= 'A') && (r <= 'Z')) || ((r >= 'a') && (r <= 'z')))
          { /* Looks like a header line: */
            ungetc(r, rd);
            char *name = fget_string(rd);
            fget_skip_and_match(rd, "=");
            fget_skip_spaces(rd);         
            neuromat_eeg_header_read_field_value(rd, name, h);
            fget_eol(rd);
            nh++;
          }
        else
          { ungetc(r, rd);
            break;
          }
        nl++;
      }
    if (nlP != NULL) { (*nlP) = nl; }
    return h;
    
    void skip_line(void)
      { 
        int r;
        do { r = fgetc(rd); } while ((r != EOF) && (r != '\n'));
        if (r == EOF) { fprintf(stderr, "** EOF found while skipping line %d\n", nl+1); exit(1); } 
      }

    /* Provide essential defaults: */
    if (isnan(h->fsmp)) { h->fsmp = 600.0; }
    if (h->ne < 0) 
      { int ne = (h->nc > 0 ? h->nc - 1 : neDef);
        fprintf(stderr, "assuming ne = %d\n", ne);
        h->ne = ne;
      }
    if (h->nc < 0) 
      { int nc = h->ne + 1;
        fprintf(stderr, "assuming nc = %d\n", nc);
        h->nc = nc;
      }
    if (h->chname == NULL) 
      { fprintf(stderr, "providing standard channel names for %d electrodes (without extra events)\n", h->ne);
        int nc1;
        neuromat_eeg_get_channel_names(h->ne, 0, NULL, &nc1, &(h->chname));
        if (nc1 != h->nc)
          { fprintf(stderr, "** inconsistent channel count: h.ne = %d  h.nc = %d  nc1 = %d\n", h->ne, h->nc, nc1); 
            exit(1);
          }
      }
  }

void neuromat_eeg_header_merge_int(int *dst, int src, char *name);
void neuromat_eeg_header_merge_double(double *dst, double src, char *name);
void neuromat_eeg_header_merge_string(char **dst, char *src, char *name);
void neuromat_eeg_header_merge_strings(int n, char ***dst, char **src, char *name);

void neuromat_eeg_header_merge(neuromat_eeg_header_t *dst, neuromat_eeg_header_t *src)
  {
    /* Dataset size and content: */
    neuromat_eeg_header_merge_int(&(dst->nt), src->nt, "nt");
    neuromat_eeg_header_merge_int(&(dst->nc), src->nc, "nc");

    /* Nature of data: */
    neuromat_eeg_header_merge_double(&(dst->fsmp), src->fsmp, "fsmp");
    neuromat_eeg_header_merge_int(&(dst->ne), src->ne, "ne");
    neuromat_eeg_header_merge_string(&(dst->type), src->type, "type");
    neuromat_eeg_header_merge_strings(dst->nc, &(dst->chname), src->chname, "chname");
    neuromat_eeg_header_merge_int(&(dst->kfmax), src->kfmax, "kfmax");
    neuromat_eeg_header_merge_string(&(dst->component), src->component, "component");
    
    /* Frequency filtering data: */
    neuromat_eeg_header_merge_int(&(dst->tdeg), src->tdeg, "tdeg");     
    neuromat_eeg_header_merge_int(&(dst->tkeep), src->tdeg, "tkeep");     
    neuromat_eeg_header_merge_double(&(dst->flo0), src->flo0, "flo0");     
    neuromat_eeg_header_merge_double(&(dst->flo1), src->flo1, "flo1");     
    neuromat_eeg_header_merge_double(&(dst->fhi1), src->fhi1, "fhi1");     
    neuromat_eeg_header_merge_double(&(dst->fhi0), src->fhi0, "fhi0");     
    neuromat_eeg_header_merge_int(&(dst->finvert), src->finvert, "finvert");

    neuromat_eeg_header_merge_orig(dst->orig, src->orig);
  }

void neuromat_eeg_header_merge_orig(neuromat_eeg_source_t *dst, neuromat_eeg_source_t *src)
  { neuromat_eeg_header_merge_string(&(dst->file), src->file, "orig.file");    
    neuromat_eeg_header_merge_int(&(dst->nt), src->nt, "orig.nt");        
    neuromat_eeg_header_merge_int(&(dst->it_ini), src->it_ini, "orig.it_ini");
    neuromat_eeg_header_merge_int(&(dst->it_fin), src->it_fin, "orig.it_fin");
    neuromat_eeg_header_merge_double(&(dst->fsmp), src->fsmp, "orig.fsmp"); 
    neuromat_eeg_header_merge_int(&(dst->subject), src->subject, "orig.subject"); 
    neuromat_eeg_header_merge_int(&(dst->run), src->run, "orig.run");     
  }

void neuromat_eeg_header_merge_int(int *dst, int src, char *name)
  { 
    if (src == INT_MIN) { return; }
    if ((*dst) == INT_MIN)
      { (*dst) = src; }
    else
      { if ((*dst) != src) 
          { fprintf(stderr, "** mismatch in field %s: dst = %d  src = %d\n", name, (*dst), src);
            assert(FALSE);
          }
      }
  }

void neuromat_eeg_header_merge_double(double *dst, double src, char *name)
  { 
    if (isnan(src)) { return; }
    if (isnan(*dst))
      { (*dst) = src; }
    else
      { if ((*dst) != src) 
          { fprintf(stderr, "** mismatch in field %s: dst = %24.15e  src = %25.15e\n", name, (*dst), src);
            assert(FALSE);
          }
      }
  }

void neuromat_eeg_header_merge_string(char **dst, char *src, char *name)
  { 
    if (src == NULL) { return; }
    if ((*dst) == NULL)
      { (*dst) = txtcat(src, ""); }
    else
      { if (strcmp((*dst), src) != 0) 
          { fprintf(stderr, "** mismatch in field %s: dst = %s  src = %s\n", name, (*dst), src);
            assert(FALSE);
          }
      }
  }

void neuromat_eeg_header_merge_strings(int n, char ***dst, char **src, char *name)
  { 
    if (src == NULL) { return; }
    assert(n != INT_MIN);
    int i;
    if ((*dst) == NULL)
      { char **cop = notnull(malloc(n*sizeof(char *)), "no mem");
        for (i = 0; i < n; i++) 
          { assert(src[i] != NULL);
            cop[i] = txtcat(src[i], "");
          }
        (*dst) = cop;
      }
    else
      { char **old = (*dst);
        for (i = 0; i < n; i++) 
        if (strcmp(old[i], src[i]) != 0) 
          { fprintf(stderr, "** mismatch in field %s[%d]: dst = %s  src = %s\n", name, i, old[i], src[i]);
            assert(FALSE);
          }
      }
  }

#define neh_NT_MAX 10000000
#define neh_NC_MAX 200
#define neh_TDEG_MAX 16

void neuromat_eeg_header_read_field_value(FILE *rd, char *name, neuromat_eeg_header_t *h)
  {
    if (strcmp(name, "nt") == 0) 
      { h->nt = fget_int(rd); 
        demand((h->nt > 0) && (h->nt <= neh_NT_MAX), "invalid {nt} in file header");
      }
    else if (strcmp(name, "nc") == 0) 
      { h->nc = fget_int(rd); 
        demand((h->nc > 0) && (h->nc <= neh_NC_MAX), "invalid {nc} in file header");
      }
    else  if (strcmp(name, "kfmax") == 0) 
      { h->kfmax = fget_int(rd); 
        demand((h->kfmax >= 0) && (h->kfmax <= neh_NT_MAX/2), "invalid {kfmax} in file header");
      }
    else if (strcmp(name, "ne") == 0) 
      { h->ne = fget_int(rd); 
        demand((h->ne > 0) && (h->ne <= neh_NC_MAX), "invalid {ne} in file header");
      }
    else if (strcmp(name, "type") == 0) 
      { h->type = fget_string(rd); }
    else if (strcmp(name, "component") == 0) 
      { h->component = fget_string(rd); }
    else if (strcmp(name, "channels") == 0) 
      { /* The number of channels {h->nc} must be known: */
        demand(h->nc != INT_MIN, "{channels} in file header before {nc}");
        h->chname = notnull(malloc(h->nc*sizeof(char*)), "no mem");
        int ic = 0;
        while(TRUE)
          { fget_skip_spaces(rd);
            int r = fgetc(rd);
            demand(r != EOF, "unexpected end-of-file");
            ungetc(r, rd); 
            if (r == '\n') { break; }
            demand(ic < h->nc, "too many channel names");
            h->chname[ic] = fget_string(rd);
            /* !!! should check valid syntax. !!! */
            ic++;
          }
        demand(ic == h->nc, "too few channel names");
      }
    else if (strcmp(name, "fsmp") == 0) 
      { h->fsmp = fget_double(rd);
        demand((! isnan(h->fsmp)) && (h->fsmp > 0), "invalid sampling frequency");
      }
    else if (strcmp(name, "trend") == 0) 
      { h->tdeg = fget_int(rd);
        h->tkeep = fget_int(rd);
        demand((h->tdeg >= -1) && (h->tdeg <= neh_TDEG_MAX), "invalid trend degree");
        demand((h->tkeep == 0) || (h->tkeep == 1), "invalid trend kept flag");
      }
    else if (strcmp(name, "band") == 0) 
      { h->flo0 = fget_double(rd); 
        h->flo1 = fget_double(rd); 
        h->fhi1 = fget_double(rd); 
        h->fhi0 = fget_double(rd); 
        h->finvert = fget_int(rd); 
        demand((h->flo0 <= h->flo1) && (h->flo1 <= h->fhi1) && (h->fhi1 <= h->fhi0), "invalid band limits");
        demand((h->flo1 == -INF) || (h->flo0 >= 0), "invalid low cutoff");
        demand((h->fhi1 == +INF) || (h->fhi0 <= h->fsmp/2), "invalid high cutoff");
        demand((h->finvert == 0) || (h->finvert == 1), "invalid filter invert flag");
      }
    else if (isprefix("orig.", name)) 
      { char *subname = strchr(name, '.');
        subname++;
        neuromat_eeg_source_read_field(rd, subname, h->orig);
      }
    else
      { fprintf(stderr, "** name = \"%s\"\n", name);
        demand(FALSE, "invalid header field name");
      }
  }

void neuromat_eeg_header_write_field_string(FILE *wr, char *pref, char *name, char *value)
  { 
    if (value != NULL) 
      { int m = (int)strlen(value);
        demand(m > 0, "empty string field value");
        demand(strcspn(value, " \240\011\012\013\014\015") == m, "invalid char in string value");
        fprintf(wr, "%s%s = %s\n", pref, name, value);
      }
  }

void neuromat_eeg_header_write_field_int_range(FILE *wr, char *pref, char *name, int vini, int vfin, int vmin, int vmax)
  { if ((vini != INT_MIN) || (vfin != INT_MIN)) 
      { demand((vmin <= vini) && (vini <= vfin) && (vfin <= vmax), "invalid index range field value");
        fprintf(wr, "%s%s = %d %d\n", pref, name, vini, vfin); 
      }  
  }

void neuromat_eeg_header_write_field_double(FILE *wr, char *pref, char *name, double value, double vmin, double vmax)
  { if (! isnan(value))    
      { demand((vmin <= value) && (value <= vmax), "invalid double field value");
        fprintf(wr, "%s%s = %.15g\n", pref, name, value);
      }  
  }
void neuromat_eeg_header_write_field_int(FILE *wr, char *pref, char *name, int value, int vmin, int vmax)
  { if (value != INT_MIN) 
      { demand((vmin <= value) && (value <= vmax), "invalid int field value");
        fprintf(wr, "%s%s = %d\n", pref, name, value);
      }  
  }
    
