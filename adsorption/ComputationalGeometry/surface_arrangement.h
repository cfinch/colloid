#include <list>
#include <cstdlib>

// CGAL
#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>

// Boost.python
#include "boost/python/extract.hpp"
#include "boost/python/numeric.hpp"

//typedef CGAL::Lazy_exact_nt<double>                Lazy_exact_nt;
typedef CGAL::Exact_circular_kernel_2                    Kernel;
//typedef CGAL::Cartesian<Lazy_exact_nt>                  Kernel;

typedef Kernel::Point_2                                 Point_2;
typedef Kernel::Vector_2                                Vector_2;
typedef Kernel::Circle_2                                Circle_2;
typedef Kernel::Segment_2                               Segment_2;
typedef Kernel::Line_2                               Line_2;
typedef Kernel::Circular_arc_point_2                    Circular_arc_point_2;

typedef CGAL::Gps_circle_segment_traits_2<Kernel>       Traits_2;
typedef CGAL::General_polygon_set_2<Traits_2>           General_polygon_set_2;
typedef CGAL::Polygon_2<Kernel>                         CGAL_polygon_2;

typedef Traits_2::Polygon_2                             Polygon_2;
typedef Traits_2::General_polygon_2                     General_polygon_2;
typedef Traits_2::Polygon_with_holes_2                  Polygon_with_holes_2;
typedef Traits_2::General_polygon_with_holes_2          General_polygon_with_holes_2;
typedef Traits_2::Curve_2                               Curve_2;
typedef Traits_2::X_monotone_curve_2                    X_monotone_curve_2;

const double pi = std::atan(1.0) * 4;

struct Arc {
  double center_x, center_y, radius, theta_source, theta_target, source_x,
         source_y, target_x, target_y;
};

struct Vertex {
    double x, y;
};

class Visualizer {
  private:
    std::list<struct Arc> curves;
    double compute_angle(double x, double y);

  public:
    void operator<<(X_monotone_curve_2 curve);
    void dummy(void);
    boost::python::list get_curves();
};

class surface_arrangement {
    private:
    bool debug;     // Controls debugging output
    std::vector<General_polygon_2> adsorbed;          // Circles representing adsorbed particles
    General_polygon_set_2 joined_particles;           // Union of polygons
    General_polygon_set_2 joined_images;              // Union of polygons
    std::list<General_polygon_with_holes_2> image_intersections;     // intersection of images with particles
    std::vector<struct Vertex> image_intersection_poly;

    int num_adsorbed, num_images;
    double radius_squared, domain_width, domain_height;

    // Construct a polygon from a circle.
    General_polygon_2 construct_polygon (const Circle_2& circle);

    // Compute the area of a general polygon, without holes, constructed from arbitrary circular arcs
    double general_polygon_area(General_polygon_2 general_polygon, bool hole);

    // Create a single particle, and one or two image particles if needed.
     void create_particle(double x, double y, double radius);

    // Compute intersection of particles and images
    void intersect_images(void);

    boost::python::list sets_to_arcs(std::list<General_polygon_with_holes_2> joined_polygon_list);

    public:
    surface_arrangement();

    // Add more particles to an existing domain     
     void add_particles(boost::python::numeric::array x_array, boost::python::numeric::array y_array, double radius);

    // Initialize domain and set the initial particle configuration
    void set_particles(boost::python::numeric::array x_array, boost::python::numeric::array y_array, double radius, double width, double height);

    // Print the configuration of adsorbed particles
    void print_particles();

    // Compute the area of the domain which is unavailable to the centers of adsorbing particles
    double unavailable_area();

    double intersection_area(General_polygon_2);

    // Return lists of arcs to Python for plotting
    boost::python::list get_arcs();
    boost::python::list get_image_arcs();
    boost::python::list get_intersection_arcs();
    boost::python::list get_poly_vertices();

    // Set debug output level
    void set_debug_flag(bool debug_output);
};


