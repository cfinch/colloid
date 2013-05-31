#include "surface_arrangement.h"

// Return a list of polygon vertices to Python for plotting
boost::python::list surface_arrangement::get_poly_vertices() {
  boost::python::list poly_list;

  for (int i=0; i<image_intersection_poly.size(); i++) {
      poly_list.append(image_intersection_poly[i]);
  }
  return poly_list;
}

// Return a list of circular arcs to Python for plotting
boost::python::list surface_arrangement::get_arcs() {
  std::list<General_polygon_with_holes_2> joined_polygon_list;    // List of polygons    
  // Copy polygons to a list
  joined_particles.polygons_with_holes(std::back_inserter (joined_polygon_list));

  return sets_to_arcs(joined_polygon_list);
}

// Return a list of image arcs to Python for plotting
boost::python::list surface_arrangement::get_image_arcs() {
  std::list<General_polygon_with_holes_2> joined_polygon_list;    // List of polygons    
  // Copy polygons to a list
  joined_images.polygons_with_holes(std::back_inserter (joined_polygon_list));

  return sets_to_arcs(joined_polygon_list);
}

// Return a list of intersection arcs to Python for plotting
boost::python::list surface_arrangement::get_intersection_arcs() {
  return sets_to_arcs(image_intersections);
}

// Convert a list of general polygons (with holes) to a list of arcs
boost::python::list surface_arrangement::sets_to_arcs(
        std::list<General_polygon_with_holes_2> joined_polygon_list) {
  Visualizer vis;
  std::list<General_polygon_with_holes_2>::iterator gp2;

  for (gp2 = joined_polygon_list.begin(); gp2 != joined_polygon_list.end(); ++gp2) {

    General_polygon_2 gp = gp2->outer_boundary();
    General_polygon_2::Curve_iterator curve;      // iterates over X_monotone_curve_2

    for (curve=gp.curves_begin(); curve!=gp.curves_end(); ++curve) {
        vis << (*curve);
    }

    std::list<General_polygon_2>::iterator hole;
    for (hole=gp2->holes_begin(); hole!=gp2->holes_end(); ++hole) {

      General_polygon_2::Curve_iterator hole_curve;
      hole->curves_begin();
      for (hole_curve=hole->curves_begin(); hole_curve!=hole->curves_end(); ++hole_curve) {
        vis << (*hole_curve);
      }
    }
  }
  return vis.get_curves();
}

// Print the configuration of adsorbed particles
void surface_arrangement::print_particles() {
    std::cout << "Adsorbed particles:" << std::endl;
    for (int i=0; i<adsorbed.size(); i++) {
        std::cout << i << ":" << adsorbed[i] << std::endl;
    }
}

// Set debug output level
void surface_arrangement::set_debug_flag(bool debug_output) {
    debug = debug_output;
}

