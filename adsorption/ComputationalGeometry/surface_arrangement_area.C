#include "surface_arrangement.h"

// Compute the area of the domain which is unavailable to the centers
// of adsorbing particles
double surface_arrangement::unavailable_area() {
  double area, poly_area, hole_area, overlap_area;
  std::list<General_polygon_with_holes_2> joined_polygon_list;    // Union of polygons

  intersect_images();

  // Copy polygons to a list
  joined_particles.polygons_with_holes(std::back_inserter (joined_polygon_list));

  // Compute area
  area = 0.0;      

  int counter = 0;
  std::list<General_polygon_with_holes_2>::iterator it;
  for (it = joined_polygon_list.begin(); it != joined_polygon_list.end(); ++it) {

    // Compute total enclosed area
    poly_area = general_polygon_area(it->outer_boundary(), false);
    area = area + poly_area;
    if (debug) std::cout << "Polygon " << counter << " with area " << poly_area << std::endl;

    // Subtract area of holes
    std::list<General_polygon_2>::iterator hole;
    for (hole=it->holes_begin(); hole!=it->holes_end(); ++hole) {
      hole_area = general_polygon_area(*hole, true);
      area = area - hole_area;
      if (debug) std::cout << "  Hole with area " << hole_area << std::endl;
    }
    counter++;
  }
      
  // Subtract area of overlap between image particles and actual particles
  std::list<General_polygon_with_holes_2>::iterator overlap;
  for (overlap=image_intersections.begin(); overlap!=image_intersections.end(); ++overlap) {
    overlap_area = intersection_area(overlap->outer_boundary());
    area = area - overlap_area;
    if (debug) std::cout << "  Image overlap with area " << overlap_area << std::endl;
  }

  return area;
} // unavailable_area

// Compute the area of a general polygon, without holes, constructed from arbitrary circular arcs
double surface_arrangement::general_polygon_area(General_polygon_2 general_polygon, bool hole) {
  CGAL_polygon_2 poly;
  double squared_radius, arc_area, poly_area, arg, area = 0.0;
  Point_2 center, start, end;
  Vector_2 v;

  General_polygon_2::Curve_iterator curve;      // iterates over X_monotone_curve_2

  for (curve=general_polygon.curves_begin(); curve!=general_polygon.curves_end(); ++curve) {
      squared_radius = CGAL::to_double(curve->supporting_circle().squared_radius());
      center = curve->supporting_circle().center();

      start = Point_2(CGAL::to_double(curve->source().x()), CGAL::to_double(curve->source().y()));
      end = Point_2(CGAL::to_double(curve->target().x()), CGAL::to_double(curve->target().y()));
      v = end-start;

      // Ensure that rounding error does not result in an argument to asin that will result in NaN
      arg = std::sqrt(CGAL::to_double(v.squared_length()))/2/std::sqrt(squared_radius);
      if (arg > 1.0) arg = 1.0;
      arc_area = squared_radius*std::asin(arg);
      area+=arc_area;

      poly.push_back(start);
      poly.push_back(center);
  }

  poly_area = CGAL::to_double(poly.area());
  if (poly_area < 0) poly_area = poly_area * -1;

  if (not hole) area = area + poly_area;
  else area = poly_area - area;
  
  return area;
}

// Compute the area of the intersection between image and actual particles
double surface_arrangement::intersection_area(General_polygon_2 general_polygon) {
  CGAL_polygon_2 poly;
  double squared_radius, theta, segment_area, poly_area, arg, area = 0.0;
  double x,y;
  Point_2 center, start, end;
  Vector_2 v;
  struct Vertex vertex;

  image_intersection_poly.clear();

  General_polygon_2::Curve_iterator curve;      // iterates over X_monotone_curve_2

  for (curve=general_polygon.curves_begin(); curve!=general_polygon.curves_end(); ++curve) {
      squared_radius = CGAL::to_double(curve->supporting_circle().squared_radius());
      center = curve->supporting_circle().center();

      x = CGAL::to_double(curve->source().x());
      y = CGAL::to_double(curve->source().y());
      vertex.x = x;
      vertex.y = y;
      start = Point_2(x,y);
      image_intersection_poly.push_back(vertex);

      x = CGAL::to_double(curve->target().x());
      y = CGAL::to_double(curve->target().y());
      end = Point_2(x,y);

      v = end-start;

      // Ensure that rounding error does not result in an argument to asin that will result in NaN
      arg = std::sqrt(CGAL::to_double(v.squared_length()))/2/std::sqrt(squared_radius);
      if (arg > 1.0) {
          std::cout << "    Arg=" << arg << std::endl;
          arg = 1.0;
      }
      theta = 2*std::asin(arg);
      segment_area = 0.5*squared_radius*(theta - sin(theta));
      area+=segment_area;

      poly.push_back(start);

//      std::cout << "       theta=" << theta << "    arc area=" << segment_area << std::endl;
  }
//  std::cout << "           total circular segment area " << area << std::endl;
//  std::cout << poly << std::endl;
  poly_area = CGAL::to_double(poly.area());
//  std::cout << "           poly area " << poly_area << std::endl;

  area = area + poly_area;
  
  return area;
}

