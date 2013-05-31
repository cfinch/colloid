#include "surface_arrangement.h"

double Visualizer::compute_angle(double x, double y) {
  double theta;

  if ((x>0) and (y>=0)) theta = std::atan(y/x);
  else if ((x<=0) and (y>0)) theta = std::atan(std::abs(x)/y) + pi/2;
  else if ((x<0) and (y<=0)) theta = std::atan(y/x)+ pi;
  else theta = std::atan(x/std::abs(y)) + 3*pi/2;

  return theta*180/pi;
}
    
void Visualizer::operator<<(X_monotone_curve_2 curve) {
  struct Arc arc;
  double x,y;     // temporary storage
  
  arc.center_x = CGAL::to_double(curve.supporting_circle().center().x());
  arc.center_y = CGAL::to_double(curve.supporting_circle().center().y());
  arc.radius = std::sqrt(CGAL::to_double(curve.supporting_circle().squared_radius()));

  x = CGAL::to_double(curve.source().x());
  y = CGAL::to_double(curve.source().y());
  arc.source_x = x;
  arc.source_y = y;
  arc.theta_source = compute_angle(x-arc.center_x,y-arc.center_y);

  x = CGAL::to_double(curve.target().x());
  y = CGAL::to_double(curve.target().y());
  arc.target_x = x;
  arc.target_y = y;
  arc.theta_target = compute_angle(x-arc.center_x,y-arc.center_y);
  
  curves.push_back(arc);
}

boost::python::list Visualizer::get_curves() {
  boost::python::list python_curves;

  std::list<Arc>::iterator curve;
  for (curve=curves.begin(); curve!=curves.end(); ++curve) {
    python_curves.append(*curve);
  }
  return python_curves;
}

