#!/usr/bin/env python
import unittest
from numpy import *
import surface
from Surface_Reaction_Tools.simulation import plot_arcs, plot_adsorbed_circles
import pylab

class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.width = 10.0
        self.radius = 1.0
        self.effective_radius = 2*self.radius
        self.delta = 1e-9
        
    def tearDown(self):
        import pylab
#        pylab.show()

    def segment_area(self, radius, distance):
        """Returns the area of a circular segment. Arguments:
        radius      obvious
        distance    distance from center of circle to segment"""
        theta = 2*arccos(distance/radius)
        return 0.5*radius**2*(theta-sin(theta))

    def overlap_area(self, radius, distance):
        """Returns the area of the intersection of two identical
        overlapping circles. Arguments:
        radius      radius of both circles
        distance    center-to-center distance"""
        theta = 2*arccos(distance/2/radius)
        return radius**2*(theta-sin(theta))

class TestCircles(BaseTestCase):
    def testOneCircle(self):
        """One circle, entirely within the domain"""
        adsorbed_x = array([5.0]);
        adsorbed_y = array([5.0]);

        arrangement = surface.surface_arrangement()
        arrangement.set_particles(adsorbed_x, adsorbed_y, self.effective_radius, self.width, self.width)
        unv_area = arrangement.unavailable_area()

        arc_list = arrangement.get_arcs()
        image_arc_list = arrangement.get_image_arcs()
        intersection_arc_list = arrangement.get_intersection_arcs()
        plot_arcs(arc_list, image_arc_list, intersection_arc_list, self.width, self.width)
        pylab.title("One circle, entirely within the domain")

        assert abs(unv_area - pi*self.effective_radius**2) < self.delta, 'Incorrect area for one particle.'

    def testOneImageCircleY(self):
        """One circle, overlapping the y edge of the domain"""
        adsorbed_x = array([5.0])
        adsorbed_y = array([self.width])

        arrangement = surface.surface_arrangement()
        arrangement.set_particles(adsorbed_x, adsorbed_y, self.effective_radius, self.width, self.width)
        unv_area = arrangement.unavailable_area()

        arc_list = arrangement.get_arcs()
        image_arc_list = arrangement.get_image_arcs()
        intersection_arc_list = arrangement.get_intersection_arcs()
        plot_arcs(arc_list, image_arc_list, intersection_arc_list, self.width, self.width)
        pylab.title("One circle, overlapping y edge")

        assert abs(unv_area - pi*self.effective_radius**2) < self.delta, \
                'Incorrect area for one particle with image.'

    def testOneImageCircleX(self):
        """One circle, overlapping the x edge of the domain"""
        adsorbed_x = array([self.width])
        adsorbed_y = array([5.0])

        arrangement = surface.surface_arrangement()
        arrangement.set_particles(adsorbed_x, adsorbed_y, self.effective_radius, self.width, self.width)
        unv_area = arrangement.unavailable_area()

        arc_list = arrangement.get_arcs()
        image_arc_list = arrangement.get_image_arcs()
        intersection_arc_list = arrangement.get_intersection_arcs()
        plot_arcs(arc_list, image_arc_list, intersection_arc_list, self.width, self.width)
        pylab.title("Once circle, overlapping x edge")

        assert abs(unv_area - pi*self.effective_radius**2) < self.delta, \
                'Incorrect area for one particle with image.'

    def testOneImageCircleXY(self):
        """One circle, overlapping the x and y edges of the domain"""
        adsorbed_x = array([10.0])
        adsorbed_y = array([0.0])

        arrangement = surface.surface_arrangement()
        arrangement.set_particles(adsorbed_x, adsorbed_y, self.effective_radius, self.width, self.width)
        unv_area = arrangement.unavailable_area()

        arc_list = arrangement.get_arcs()
        image_arc_list = arrangement.get_image_arcs()
        intersection_arc_list = arrangement.get_intersection_arcs()
        plot_arcs(arc_list, image_arc_list, intersection_arc_list, self.width, self.width)
        pylab.title("Once circle, overlapping x and y edges")

        assert abs(unv_area - pi*self.effective_radius**2) < self.delta, \
                'Incorrect area for one particle with image.'

    def testTwoTouchingCircles(self):
        """Two touching circles at the y periodic BC"""
        delta_x = 2.0
        adsorbed_x = array([5.0, 5.0])
        adsorbed_y = array([2.0, self.width-2.0])

        arrangement = surface.surface_arrangement()
        arrangement.set_particles(adsorbed_x, adsorbed_y, self.effective_radius, self.width, self.width)
        unv_area = arrangement.unavailable_area()

        arc_list = arrangement.get_arcs()
        image_arc_list = arrangement.get_image_arcs()
        intersection_arc_list = arrangement.get_intersection_arcs()
        plot_arcs(arc_list, image_arc_list, intersection_arc_list, self.width, self.width)
        pylab.title("Two touching circles at the y periodic BC")

        actual_area = 2*pi*self.effective_radius**2
        assert abs(unv_area - actual_area) < self.delta, \
            'Incorrect area for two touching circles '+str(unv_area) + "   " + str(actual_area)

    def testThreeTouchingCircles(self):
        """Three circles that just touch.  Tests the ability to exclude the area of the hole
        formed by the three circles."""
        adsorbed_x = array([3.0+self.delta, 3.0+self.effective_radius, 3.0+2*self.effective_radius])
        adsorbed_y = array([3.0, 3.0+sqrt(3)*self.effective_radius, 3.0])

        arrangement = surface.surface_arrangement()
        arrangement.set_particles(adsorbed_x, adsorbed_y, self.effective_radius, self.width, self.width)
        unv_area = arrangement.unavailable_area()

        arc_list = arrangement.get_arcs()
        image_arc_list = arrangement.get_image_arcs()
        intersection_arc_list = arrangement.get_intersection_arcs()
        plot_arcs(arc_list, image_arc_list, intersection_arc_list, self.width, self.width)
        pylab.title("Three touching circles")

        actual_area = 3*pi*self.effective_radius**2
        assert abs(unv_area - actual_area) < self.delta, \
            'Incorrect area for three touching Circles '+str(unv_area) + "   " + str(actual_area)

    def testTwoOverlappingCircles(self):
        """Two overlapping circles, entirely within the domain"""
        delta_x = 2.0
        adsorbed_x = array([5.0, 5.0])
        adsorbed_y = array([5.0, 5.0+delta_x])

        arrangement = surface.surface_arrangement()
        arrangement.set_particles(adsorbed_x, adsorbed_y, self.effective_radius, self.width, self.width)
        unv_area = arrangement.unavailable_area()

        arc_list = arrangement.get_arcs()
        image_arc_list = arrangement.get_image_arcs()
        intersection_arc_list = arrangement.get_intersection_arcs()
        plot_arcs(arc_list, image_arc_list, intersection_arc_list, self.width, self.width)
        pylab.title("Two overlapping circles")

        actual_area = 2*pi*self.effective_radius**2 - self.overlap_area(self.effective_radius, delta_x)
        assert abs(unv_area - actual_area) < self.delta, \
            'Incorrect area for two overlapping circles '+str(unv_area) + "   " + str(actual_area)

    def testTwoOverlappingImageCircles(self):
        """Circles at opposite edges of the domain, set so that they overlap the same amount
        as the previous test case."""
        delta_x = 2.0
        adsorbed_x = array([5., 5.])
        adsorbed_y = array([delta_x/2, self.width-delta_x/2])

        arrangement = surface.surface_arrangement()
        arrangement.set_particles(adsorbed_x, adsorbed_y, self.effective_radius, self.width, self.width)
        unv_area = arrangement.unavailable_area()

        arc_list = arrangement.get_arcs()
        image_arc_list = arrangement.get_image_arcs()
        intersection_arc_list = arrangement.get_intersection_arcs()
        plot_arcs(arc_list, image_arc_list, intersection_arc_list, self.width, self.width)
        pylab.title("Two overlapping circles at y edges")

        actual_area = 2*pi*self.effective_radius**2 - self.overlap_area(self.effective_radius, delta_x)
        assert abs(unv_area - actual_area) < self.delta, \
            'Incorrect area for two overlapping image circles '+str(unv_area) + "   " + str(actual_area)

    def testThreeOverlappingImageCircles(self):
        """Three overlapping circles with images."""

        delta_x = 2*sqrt(3)
        delta_y = 2.0
        adsorbed_x = array([5., 5-delta_x/2, 5+delta_x/2])
        adsorbed_y = array([delta_y/2+1, self.width-delta_y/2, self.width-delta_y/2])

        arrangement = surface.surface_arrangement()
        arrangement.set_particles(adsorbed_x, adsorbed_y, self.effective_radius, self.width, self.width)
        unv_area = arrangement.unavailable_area()

        arc_list = arrangement.get_arcs()
        image_arc_list = arrangement.get_image_arcs()
        intersection_arc_list = arrangement.get_intersection_arcs()
        plot_arcs(arc_list, image_arc_list, intersection_arc_list, self.width, self.width)
        pylab.title("Three overlapping circles at y edges")

        actual_area = 3*pi*self.effective_radius**2 - 3*self.overlap_area(self.effective_radius, delta_x)

        assert abs(unv_area - actual_area) < self.delta, \
            'Incorrect area for three overlapping image circles '+str(unv_area) + "   " + str(actual_area)

    def testHexagonalCirclePackedParticles(self):
        from Surface_Reaction_Tools.arrangements import generate_hex_circle_packing
        adsorbed_x, adsorbed_y, x_size, y_size = generate_hex_circle_packing(self.radius, self.width)

        arrangement = surface.surface_arrangement()
        arrangement.set_particles(adsorbed_x, adsorbed_y, self.effective_radius, x_size, y_size)
        unv_area = arrangement.unavailable_area()

        arc_list = arrangement.get_arcs()
        image_arc_list = arrangement.get_image_arcs()
        intersection_arc_list = arrangement.get_intersection_arcs()
#        plot_adsorbed_circles(adsorbed_x, adsorbed_y, self.effective_radius, self.width)
        plot_arcs(arc_list, image_arc_list, intersection_arc_list, x_size, y_size)
        pylab.title("Hex Packed Circles")

        actual_area = x_size*y_size
        assert abs(unv_area - actual_area) < self.delta, \
            'Incorrect area for hex packed particles.\n' \
            + 'x size:' + str(x_size) + '  y size:' + str(y_size) + '\n' \
            + 'Calculated area:' + str(unv_area) + " Actual area:" + str(actual_area)
    ####
####

if __name__ == "__main__":
    unittest.main()

