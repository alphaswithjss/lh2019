#ifndef __SHAPE_HH__
#define __SHAPE_HH__

/// \class Shape
/// empty base class for the other shapes
class Shape{
public:
  Shape(){}
  virtual ~Shape(){}
  virtual double density(double v, double pt, bool gluon){ return 0.0;}
  virtual double fraction_below(double lambda, double pt, bool gluon){ return 0.0;}
};

#endif


