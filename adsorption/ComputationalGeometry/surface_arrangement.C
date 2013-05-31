#include "surface_arrangement.h"

// Constructor
surface_arrangement::surface_arrangement() {
    debug = false;
}

// Construct a polygon from a circle.
General_polygon_2 surface_arrangement::construct_polygon (const Circle_2& circle)
{
  // Subdivide the circle into two x-monotone arcs.
  Traits_2 traits;
  Curve_2 curve (circle);
  std::list<CGAL::Object> objects;
  traits.make_x_monotone_2_object() (curve, std::back_inserter(objects));
  CGAL_assertion (objects.size() == 2);

  // Construct the polygon.
  General_polygon_2 pgn;
  X_monotone_curve_2 arc;
  std::list<CGAL::Object>::iterator iter;

  for (iter = objects.begin(); iter != objects.end(); ++iter) {
    CGAL::assign (arc, *iter);
    pgn.push_back (arc);
  }

  return pgn;
}


// Compute the intersections between particles and images
void surface_arrangement::intersect_images() {
    std::list<General_polygon_with_holes_2>::iterator it, im;
    std::list<General_polygon_with_holes_2> joined_polygon_list, joined_images_list; // Union of polygons

    image_intersections.clear();

    // Copy polygons to a list
    joined_particles.polygons_with_holes(std::back_inserter (joined_polygon_list));
    joined_images.polygons_with_holes(std::back_inserter (joined_images_list));

    for (it = joined_polygon_list.begin(); it != joined_polygon_list.end(); ++it) {
      // intersect with image particles
      for (im = joined_images_list.begin(); im != joined_images_list.end(); ++im) {
        CGAL::intersection(it->outer_boundary(), im->outer_boundary(), std::back_inserter(image_intersections));
      }
    }
}


// Create a single particle, and one or two image particles if needed.
// Corner particles need two images
void surface_arrangement::create_particle(double x, double y, double radius) {
  double radius_squared = std::pow(radius,2);
  General_polygon_2 gp, image;

  // Add particle to list
  Circle_2 circle(Point_2(x,y), radius_squared);
  gp = construct_polygon(circle);
  adsorbed.push_back(gp);
  // Join to other adsorbed circles
  joined_particles.join(gp);

  // Add images
  if (x < 2*radius) {
    Circle_2 circle(Point_2(x+domain_width,y), radius_squared);
    image = construct_polygon(circle);
    joined_images.join(image);
    num_images++;
  }
  if (y < 2*radius) {
    Circle_2 circle(Point_2(x,y+domain_height), radius_squared);
    image = construct_polygon(circle);
    joined_images.join(image);
    num_images++;
  }
  if ((x < 2*radius) and (y < 2*radius)) {
    // Particle at bottom left has image at top right
    Circle_2 circle(Point_2(x+domain_width,y+domain_height), radius_squared);
    image = construct_polygon(circle);
    joined_images.join(image);
    num_images++;

  }
  if ((x > domain_width - 2*radius) and (y < 2*radius)) {
    // Particle at bottom right has image at top left      
    Circle_2 circle(Point_2(-(domain_width-x),y+domain_height), radius_squared);
    image = construct_polygon(circle);
    joined_images.join(image);
    num_images++;
  }
}


// Add more particles to an existing domain
void surface_arrangement::add_particles(boost::python::numeric::array x_array,
    boost::python::numeric::array y_array, double radius) {
  double x,y;
  int num_new_particles=0, num_previous_images;

  num_previous_images = num_images;
  num_new_particles = boost::python::extract<int>(x_array.attr("__len__")());

  // Add particles
  for (int i=0; i<num_new_particles; i++) {
    x = boost::python::extract<double>(x_array[i]);
    y = boost::python::extract<double>(y_array[i]);
    create_particle(x,y,radius);
  }
  num_adsorbed += num_new_particles;

  if (debug) std::cout << "Added " << num_new_particles << " particles and " 
    << num_images - num_previous_images << " images." << std::endl;
}


// Initialize domain and set the initial particle configuration
// Maybe want to make this a constructor which calls add_particles() in order 
// to eliminate duplicate code
void surface_arrangement::set_particles(boost::python::numeric::array x_array, 
    boost::python::numeric::array y_array, double radius, double width, double height) {

  radius_squared = std::pow(radius,2);
  double x,y;

  num_adsorbed = boost::python::extract<int>(x_array.attr("__len__")());
  num_images = 0;
  domain_width = width;     domain_height = height;

  // Construct list of adsorbed particles
  adsorbed.clear();
  adsorbed.resize(num_adsorbed);

  for (int i=0; i<num_adsorbed; i++) {
      x = boost::python::extract<double>(x_array[i]);
      y = boost::python::extract<double>(y_array[i]);
      create_particle(x,y,radius);
  }
  if (debug) {
      std::cout << "Created " << num_adsorbed << " particles and "
        << num_images << " images." << std::endl;
     std::cout << joined_particles.number_of_polygons_with_holes () 
        << " generalized polygons after merging." << std::endl;
  }

  // join polygons all at once
  //      CGAL::join (adsorbed.begin(), adsorbed.end(), std::back_inserter (joined_polygon_list));
}


#include <boost/python.hpp>
using namespace boost::python;

BOOST_PYTHON_MODULE(surface)
{
  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");

  class_<surface_arrangement>("surface_arrangement")
    .def("set_particles", &surface_arrangement::set_particles)
    .def("add_particles", &surface_arrangement::add_particles)
    .def("unavailable_area", &surface_arrangement::unavailable_area)
    .def("print_particles", &surface_arrangement::print_particles)
    .def("get_arcs", &surface_arrangement::get_arcs)
    .def("get_image_arcs", &surface_arrangement::get_image_arcs)
    .def("get_intersection_arcs", &surface_arrangement::get_intersection_arcs)
    .def("get_poly_vertices", &surface_arrangement::get_poly_vertices)
    .def("set_debug_flag", &surface_arrangement::set_debug_flag)
  ;

  class_<Arc>("Arc")
    .def_readonly("center_x", &Arc::center_x)
    .def_readonly("center_y", &Arc::center_y)
    .def_readonly("radius", &Arc::radius)
    .def_readonly("theta_source", &Arc::theta_source)
    .def_readonly("theta_target", &Arc::theta_target)
    .def_readonly("source_x", &Arc::source_x)
    .def_readonly("source_y", &Arc::source_y)
    .def_readonly("target_x", &Arc::target_x)
    .def_readonly("target_y", &Arc::target_y)
  ;

  class_<Vertex>("Vertex")
    .def_readonly("x", &Vertex::x)
    .def_readonly("y", &Vertex::y)
  ;

}

