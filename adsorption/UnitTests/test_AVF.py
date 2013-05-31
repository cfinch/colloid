#!/usr/bin/env python
from numpy import *
import tables
import pylab

from Surface_Reaction_Tools.analysis import compute_AVF

def plot_particles(adsorbed_x, adsorbed_y, exclusion_radius, S):
    from matplotlib.patches import Circle

    fig = pylab.figure()
    ax = fig.add_subplot(111)
    for p in range(len(adsorbed_x)):
        ax.add_patch(Circle((adsorbed_x[p], adsorbed_y[p]), exclusion_radius, color='red'))
    #for p in range(len(plot_x)):
    #    ax.add_patch(Circle((plot_x[p], plot_y[p]), 0.1))
    pylab.axis([0.0, S, 0.0, S])

a = 1.0
S = 50.0
num_replicates = 50
random_radius = a/2

delta_x = a/8
delta_y = delta_x
delta_h = 0.1
border = delta_x/2

x_grid = arange(border, S-border+delta_y/10, delta_x)
y_grid = arange(border, S-border+delta_y/10, delta_y)
h_grid = array([0])
exclusion_radius = 2*sqrt(a**2 - (h_grid[0]/2)**2)
print "Effective exclusion radius at h=" + str(h_grid[0]) + " = " + str(exclusion_radius)

## One particle
adsorbed_x = array([S*2./3.])
adsorbed_y = array([S*2./3.])
print "One particle, no images"
AVF = compute_AVF(x_grid, y_grid, h_grid, adsorbed_x, adsorbed_y, a, \
        num_replicates, random_radius)
print "  Theoretical AVF: ", str(1.0 - pi*exclusion_radius**2/S**2), "  AVF: ", str(AVF)

## Two particles, non-overlapping exclusion zones
adsorbed_x = array([S/2, S/2+10*a])
adsorbed_y = array([S/2, S/2+10*a])
print "Two particles, no overlap, no images"
AVF = compute_AVF(x_grid, y_grid, h_grid, adsorbed_x, adsorbed_y, a, \
        num_replicates, random_radius)

print "  Theoretical AVF: ", str(1.0 - 2*(pi*exclusion_radius**2/S**2)), "  AVF: ", str(AVF)
plot_particles(adsorbed_x, adsorbed_y, exclusion_radius, S)

## Two particles, overlapping exclusion zones
adsorbed_x = array([S/2, S/2+3*a])
adsorbed_y = array([S/2, S/2])
print "Two particles, overlapping, no images"
AVF = compute_AVF(x_grid, y_grid, h_grid, adsorbed_x, adsorbed_y, a, \
        num_replicates, random_radius)

delta_x = adsorbed_x[1] - adsorbed_x[0]
theta = 2*arccos(delta_x/(2*exclusion_radius))
area = 2*pi*exclusion_radius**2 - exclusion_radius**2*(theta-sin(theta))
fraction = 1.0 - area/S**2
print "  Theoretical AVF: ", str(fraction), "  AVF: ", str(AVF)
plot_particles(adsorbed_x, adsorbed_y, exclusion_radius, S)

## Hexagonal circle packing (AVF=0.0)
numParticles = 0
# Create a row
x_list = []
y_list = []
y = a
x = a
rowNumber = 0
while y <= S+1.1*a:
    # Create circles in a row
    while x <= S+a:
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
adsorbed_x = array(x_list)
adsorbed_y = array(y_list)

print "Hexagonal circle packing:"
AVF =  compute_AVF(x_grid, y_grid, h_grid, adsorbed_x, adsorbed_y, a, \
        num_replicates, random_radius)
print "  Theoretical AVF: ", str(0.0), "  Computed AVF: ", str(AVF)

## Hexagonal circle packing, 3 particles deleted
# Delete 3 particles from center of grid
adsorbed_x[int(len(adsorbed_x)/2)] = -10.0
adsorbed_y[int(len(adsorbed_y)/2)] = -10.0
adsorbed_x[int(len(adsorbed_x)/2)+1] = -10.0
adsorbed_y[int(len(adsorbed_y)/2)+1] = -10.0
adsorbed_x[int(len(adsorbed_x)/2)+27] = -10.0
adsorbed_y[int(len(adsorbed_y)/2)+27] = -10.0

print "Hexagonal circle packing gaps:"
AVF = compute_AVF(x_grid, y_grid, h_grid, adsorbed_x, adsorbed_y, a, \
        num_replicates, random_radius)
print "  Theoretical AVF: ", str((2*a)**2*(sqrt(3)/2-pi/6)/S**2), "  Computed AVF: ", str(AVF)
plot_particles(adsorbed_x, adsorbed_y, exclusion_radius, S)

pylab.show()

