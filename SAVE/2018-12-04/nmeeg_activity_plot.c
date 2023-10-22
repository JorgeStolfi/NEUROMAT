/* Plots cortical or scalp activity as cloud in cylinder. */
/* Last edited on 2018-10-25 14:19:16 by stolfilocal */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h> 
#include <math.h>

#include <frgb.h>
#include <float_image.h>
#include <affirm.h>
#include <r3.h>
#include <ppv_array.h>

typedef struct nmeeg_activity_view_t 
  { /* Primary view data: */
    r3_t obs;     /* Observer's position. */
    double ctrZ;  /* Z-coord of center of interest on Z-axis. */
    double resXY; /* Pixels per voxel in the XY direction at center of interest. */
    double resZ;  /* Pixels per voxel in the Z direction at center of interest. */
  } nmeeg_activity_view_t;
  /* Describes the viewing parameters of the activity plot. */

float_image_t nmeeg_activity_plot
  ( ppv_array_t *act,
    double rad,
    double  ht,
    int32_t tmin, 
    int32_t tmax, 
    nmeeg_activity_view_t *view, 
    nmeeg_activity_style_t *style
  );
  /* Creates a plot of cortical activity or skull potential as stored in the
  voxel array {act}. 
  
  First, the array {act} is interpreted as block of voxels in three-space, 
  aligned with the X, Y, and Z axes.  Let {NX}, {NY}, and {NZ} be the 
  number of elements of the array in the first three directions.  
  The of the array are assumed to be the X, Y, and Z
  axes of the World coordinate system.
  
  The activity is a function {func(x,y,t)} of
  position on the cortex, represented by a point {(x,y)} on the unit disk of the
  XY plane with radius {rad}; and time, represented by an integer {t} 
  (a discrete sampling time index) on the Z axis in the range {0..ht}.
  
  The function {fun} is plotted as a cloud of splatted particles according to the given {style}.
  Only the segment from {tmin} to {tmax} is plotted.
  The cylinder defined by {rad} and {ht} is projected on the image plane according to the
  given {view}.  The procedure returns a color image object."""
  
  # Create the image array and fill it with black/transparent:
  img = create_image(rad, 0, ht, view)
  
  # Draw the back part of the cylinder's frame:
  draw_cyl_frame(img, rad, 0, ht, view)
  
  # Fill the plotted area with black:
  draw_cyl_frame(img, rad, tmin, tmax, view)
  
  # Paint the activity on the image:
  for pt in range(ht+1):
    # Select the time {t} to plot, in the proper (back to front) direction:
    if view["dirT"] < 0: t = ht - pt; else: t = pt
    # Splay back to front, left to right as seen on the image:
    for pBF in range(-rad,rad+1):
      for pLR in range(-rad,rad+1):
        (x,y) = rotate_xy_plot_coords(pBF,pLR,view)
        if (x*x + y*y < rad*rad):
          fval = f(x,y,t)
          splat_dot(img, x, y, t, fval, view, style)
  
  # Draw the front part of the cylinder's frame:
  draw_cyl_frame(img, rad, ht, view)
  return img
  ######################################################################
  
def create_image(rad, ht, view):
  """Creates an image suitable for plotting a cylinder with radius {rad}
  and height {ht} seen from the given {view}.  Chooses the size of the image
  so that the cylinder will fit with a suitable margin, and the 
  resolution at the cylinder's center will be that specified by {view["resolution"]}."""
  
  # Get the XY plot directions back-to-front and left-to-right:
  backXY = rotate_xy_plot_coords(-rad,0,view)
  leftXY = rotate_xy_plot_coords(0,-rad,view)
  
  # Determine the extreme points of the projected cylinder:
  if view["dirT"] < 0:
    # Cylinder is viewed from below:
    topXYZ = (-backXY[0], -backXY[1], ht)
    leftXYZ = (leftXY[0], leftXY[1], 0)
    rightXYZ = (-leftXY[0], -leftXY[1], 0)
    botXYZ = (backXY[0], backXY[1], 0)
  else:
    # Cylinder is viewed from above:
    topXYZ = (backXY[0], backXY[1], ht)
    leftXYZ = (leftXY[0], leftXY[1], ht)
    rightXYZ = (-leftXY[0], -leftXY[1], ht)
    botXYZ = (-backXY[0], -backXY[1], 0)
  
  # Project the points to image coordinates:
