#!/usr/bin/env python
from numpy import *

def generate_hex_circle_packing(a,width):
    """Generate a domain, filled with hexagonally packed circles, of
    the given width.  The height will be determined so that the vertical 
    boundary condition is periodic.

    Arguments:
    a       particle radius
    width   domain size, in terms of particle radius

    Returns:
    x_list  list of x coordinates
    y_list  list of y coordinates
    x_size  width of domain (equal to argument width)
    y_size  height of domain
    """
    numParticles = 0

    x_list = []
    y_list = []
    y = a
    x = a
    rowNumber = 0
    # Create a row
    while y <= width*1.01:
        # Create circles in a row
        while x < width:
            x_list.append(x)
            x = x + 2*a
            y_list.append(y)
            numParticles = numParticles + 1
        y = y + a*sqrt(3.)
        rowNumber = rowNumber + 1
        if rowNumber%2 == 0:
            x = a
        else:
            x = 0
    x_size = width
    y_size = rowNumber*a*sqrt(3)

    return array(x_list), array(y_list), x_size, y_size
