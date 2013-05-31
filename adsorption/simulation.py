#!/usr/bin/env python
from numpy import *
from scipy import weave

class Setup:
    pass

class Nondimensional_Setup:
    pass

def findXYCollision(x,y,x_array,y_array,r,width):
    """Return the index of the particle that has the largest overlap with a 
    particle of radius r centered at (x,y).  Returns -1 if no overlaps exist.
    """
    l = len(x_array)

#    print x_array
    code = r"""
           double min_value = 1e9;
           int min_index = 0;
           double d2;
           double x_double = x, y_double=y;
           double x_i, y_i;

           for (int i=0; i<l; i++) {
               // Handle periodic boundary conditions
               // x sides
               if ((x < 2*r) && (x_array[i] > width-2*r)) x_i = x_array[i]-width;
               else if ((x > width-2*r) && (x_array[i] < 2*r)) x_i = width+x_array[i];
               else x_i = x_array[i];
               // y
               if ((y < 2*r) && (y_array[i] > width-2*r)) y_i = y_array[i]-width;
               else if ((y > width-2*r) && (y_array[i] < 2*r)) y_i = width+y_array[i];
               else y_i = y_array[i];

               d2 = pow(x_double-x_i, 2) + pow(y_double-y_i,2);
               if (d2 < min_value) {
                   min_value = d2;
                   min_index = i;
               }
           }
           if (min_value < 4.0*pow(r,2))
               return_val = min_index;
           else
               return_val = -1; 
    """
    v = weave.inline(code, ['x', 'y', 'x_array', 'y_array', 'r', 'l','width'], \
        compiler='gcc')
#    print v
    return v

def collision_detection(particle_indices, x_array, y_array, z_array, r, S,
        compiler='gcc'):
    """Given an array of particle indices and arrays of (x,y,z) coordinates,
    this method returns an array containing the indices of the particles that have the
    greatest overlap with the corresponding particle in particle_indices.  If no
    particle overlaps with that particle, that spot in the array is set to -1. 
    It is assumed that all particles have radius "r".
    """
    from numpy import ones, int32
    from scipy import weave

    collision_indices = ones(len(x_array), dtype=int32)*-1       # set all to "no collision"
    num_active_particles = len(particle_indices)

    num_particles = len(x_array)

    image_n_x  = x_array.copy()
    image_n_y  = y_array.copy() + S

    image_e_x  = x_array.copy() + S
    image_e_y  = y_array.copy()

    image_s_x  = x_array.copy()
    image_s_y  = y_array.copy() - S

    image_w_x  = x_array.copy() - S
    image_w_y  = y_array.copy()

    code = r"""
           double min_value;
           int p, min_index, i, ref;
           double d2, x_i, y_i, x_ref, y_ref, z_ref;
           double min_radius = 4.0*pow(r,2);

           for (ref=0; ref<num_active_particles; ref++) {
               p = particle_indices[ref];
               x_ref = x_array[p];
               y_ref = y_array[p];
               z_ref = z_array[p];

               min_value = 1e9;
               min_index = 0;

               for (i=0; i<p; i++) {
                   d2 = pow(x_ref - x_array[i], 2) + pow(y_ref - y_array[i],2) + pow(z_ref - z_array[i],2);
                   if (d2 < min_value) {
                       min_value = d2;
                       min_index = i;
                   }
               }
               for (i=p+1; i<num_particles; i++) {
                   d2 = pow(x_ref - x_array[i], 2) + pow(y_ref - y_array[i],2) + pow(z_ref - z_array[i],2);
                   if (d2 < min_value) {
                       min_value = d2;
                       min_index = i;
                   }
               }

               // Images on west and east sides
               if ((x_ref < 2*r) && (min_value >= min_radius)) {
                   for (i=0; i<num_particles; i++) {
                       d2 = pow(x_ref - image_w_x[i], 2) + pow(y_ref - image_w_y[i],2) + pow(z_ref - z_array[i],2);
                       if (d2 < min_value) {
                           min_value = d2;
                           min_index = i;
                       }
                   }
               }
               else if ((x_ref > (S - 2*r)) && (min_value >= min_radius)) {
                   for (i=0; i<num_particles; i++) {
                       d2 = pow(x_ref - image_e_x[i], 2) + pow(y_ref - image_e_y[i],2) + pow(z_ref - z_array[i],2);
                       if (d2 < min_value) {
                           min_value = d2;
                           min_index = i;
                       }
                   }
               }

               // Images on bottom and top
               if ((y_ref < 2*r) && (min_value >= min_radius)) {
                   for (i=0; i<num_particles; i++) {
                       d2 = pow(x_ref - image_s_x[i], 2) + pow(y_ref - image_s_y[i],2) + pow(z_ref - z_array[i],2);
                       if (d2 < min_value) {
                           min_value = d2;
                           min_index = i;
                       }

                   }
               }
               else if ((y_ref > (S - 2*r)) && (min_value >= min_radius)) {
                   for (i=0; i<num_particles; i++) {
                       d2 = pow(x_ref - image_n_x[i], 2) + pow(y_ref - image_n_y[i],2) + pow(z_ref - z_array[i],2);
                       if (d2 < min_value) {
                           min_value = d2;
                           min_index = i;
                       }
                   }
               }

               if (min_value < min_radius)
                   collision_indices[p] = min_index;
           }
    """
    weave.inline(code, ['particle_indices', 'x_array', 'y_array', 'z_array',
        'collision_indices', 'r', 'num_particles', 'num_active_particles',
        'image_n_x',  'image_n_y', 'image_e_x',  'image_e_y',
        'image_s_x',  'image_s_y', 'image_w_x',  'image_w_y',
        'S'], compiler=compiler)

    return collision_indices

def plot_adsorbed_circles(adsorbed_x,adsorbed_y,radius, width, label=""):
    import pylab
    from matplotlib.patches import Circle

    # Plot each run
    fig = pylab.figure()
    ax = fig.add_subplot(111)
    for p in range(len(adsorbed_x)):
        ax.add_patch(Circle((adsorbed_x[p], adsorbed_y[p]), radius))
        # Plot "image" particles to verify that periodic boundary conditions are working
#        if adsorbed_x[p] < radius:
#            ax.add_patch(Circle((adsorbed_x[p] + width,adsorbed_y[p]), radius, facecolor='red'))
#        elif adsorbed_x[p] > (width-radius):
#            ax.add_patch(Circle((adsorbed_x[p] - width,adsorbed_y[p]), radius, facecolor='red'))
#        if adsorbed_y[p] < radius:
#            ax.add_patch(Circle((adsorbed_x[p],adsorbed_y[p] + width), radius, facecolor='red'))
#        elif adsorbed_y[p] > (width-radius):
#            ax.add_patch(Circle((adsorbed_x[p],adsorbed_y[p] - width), radius, facecolor='red'))

    ax.set_aspect(1.0)
    pylab.axhline(y=0, color='k')
    pylab.axhline(y=width, color='k')
    pylab.axvline(x=0, color='k')
    pylab.axvline(x=width, color='k')
    pylab.axis([-0.1*width, width*1.1, -0.1*width, width*1.1])
    pylab.xlabel("non-dimensional x")
    pylab.ylabel("non-dimensional y")
    pylab.title("Adsorbed particles at theta="+label)

    return ax

def plot_arcs(arc_list, image_arc_list, intersection_arc_list, width, height):
    import matplotlib as mpl
    import pylab

    fig = pylab.figure()
    ax = fig.add_subplot(1,1,1)

    for arc in arc_list:
        center = (arc.center_x, arc.center_y)
        patch = mpl.patches.Arc(center,width=2*arc.radius,height=2*arc.radius,angle=0, \
        theta1=arc.theta_source, theta2=arc.theta_target, edgecolor='k')
        ax.add_patch(patch)

        ax.add_patch(mpl.patches.Circle( (arc.source_x,arc.source_y), 0.1, facecolor='r'))
        ax.add_patch(mpl.patches.Circle( (arc.target_x,arc.target_y), 0.1, facecolor='b'))

    for arc in image_arc_list:
        center = (arc.center_x, arc.center_y)
        patch = mpl.patches.Arc(center,width=2*arc.radius,height=2*arc.radius,angle=0, \
        theta1=arc.theta_source, theta2=arc.theta_target, edgecolor='r')
        ax.add_patch(patch)
        ax.add_patch(mpl.patches.Circle( (arc.source_x,arc.source_y), 0.1, facecolor='r'))
#        ax.add_patch(mpl.patches.Circle( (arc.target_x,arc.target_y), 0.1, facecolor='b'))

    for arc in intersection_arc_list:
        center = (arc.center_x, arc.center_y)
        patch = mpl.patches.Arc(center,width=2*arc.radius,height=2*arc.radius,angle=0, \
        theta1=arc.theta_source, theta2=arc.theta_target, edgecolor='g', linewidth=2)
        ax.add_patch(patch)
        ax.add_patch(mpl.patches.Circle( (arc.source_x,arc.source_y), 0.1, facecolor='r'))
#        ax.add_patch(mpl.patches.Circle( (arc.target_x,arc.target_y), 0.1, facecolor='b'))

    ax.set_aspect(1.0)
    pylab.axhline(y=0, color='k')
    pylab.axhline(y=height, color='k')
    pylab.axvline(x=0, color='k')
    pylab.axvline(x=width, color='k')
    pylab.axis([-0.1*width, width*1.1, -0.1*width, width*1.1])
    pylab.xlabel("non-dimensional x")
    pylab.ylabel("non-dimensional y")

    return ax
