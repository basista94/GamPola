#ifndef GVECTOR4D_HH
#define GVECTOR4D_HH

#include <iostream>
#include <math.h>

class GVector3D;

class GVector4D {

  

  inline friend GVector4D operator*(double d,const GVector4D& v2); 
  inline friend GVector4D operator*(const GVector4D& v2,double d); 
  inline friend GVector4D operator/(const GVector4D& v2,double d); 
  inline friend double operator*(const GVector4D& v1,const GVector4D& v2); 
  inline friend GVector4D operator+(const GVector4D& v1,const GVector4D& v2); 
  inline friend GVector4D operator-(const GVector4D& v1,const GVector4D& v2); 
  
public:
  GVector4D();
  GVector4D(double e,double px,double py ,double pz);
  inline void set(int i,double d);
  inline void set(double e,double px,double py ,double pz);
  inline GVector4D& operator*=(double c);
  inline GVector4D& operator/=(double c);
  inline GVector4D& operator=(const GVector4D& v2);
  inline GVector4D& operator+=(const GVector4D& v2);
  inline GVector4D& operator-=(const GVector4D& v2);
  inline double get(int i) const;
  inline double cont(const GVector4D& v4) const;
  friend std::ostream& operator<<(std::ostream& s, const GVector4D& v);  
  double mass2() const;     
  double mass() const;
  void applyRotateEuler(double alpha,double beta,double gamma);
  void applyBoostTo(const GVector4D& p4, bool inverse = false);
  void applyBoostTo(const GVector3D& boost, bool inverse = false);
  GVector4D cross(const GVector4D& v2);
  double dot(const GVector4D& v2) const;
  double d3mag() const;

  double dotr3( const GVector4D& p1, const GVector4D& p2 ) const;
  double mag2r3( const GVector4D& p1 ) const;
  double magr3( const GVector4D& p1 ) const;


private:

  double v[4];

  inline double Square( double x ) const { return x*x; }

};

GVector4D rotateEuler(const GVector4D& rs,
                                 double alpha,double beta,double gamma);
GVector4D boostTo(const GVector4D& rs,
                     const GVector4D& p4, bool inverse = false);
GVector4D boostTo(const GVector4D& rs,
                     const GVector3D& boost, bool inverse = false);

inline GVector4D& GVector4D::operator=(const GVector4D& v2){

  v[0]=v2.v[0];
  v[1]=v2.v[1];
  v[2]=v2.v[2];
  v[3]=v2.v[3];
  
  return *this; 
}

inline GVector4D& GVector4D::operator+=(const GVector4D& v2){

  v[0]+=v2.v[0];
  v[1]+=v2.v[1];
  v[2]+=v2.v[2];
  v[3]+=v2.v[3];
  
  return *this; 
}

inline GVector4D& GVector4D::operator-=(const GVector4D& v2){

  v[0]-=v2.v[0];
  v[1]-=v2.v[1];
  v[2]-=v2.v[2];
  v[3]-=v2.v[3];
  
  return *this; 
}

inline double GVector4D::mass2() const{

  return v[0]*v[0]-v[1]*v[1]-v[2]*v[2]-v[3]*v[3];
}

inline GVector4D operator*(double c,const GVector4D& v2){
  
  return GVector4D(v2)*=c;
}

inline GVector4D operator*(const GVector4D& v2,double c){
  
  return GVector4D(v2)*=c;
}

inline GVector4D operator/(const GVector4D& v2,double c){
  
  return GVector4D(v2)/=c;
}

inline GVector4D& GVector4D::operator*=(double c){

  v[0]*=c;  
  v[1]*=c;  
  v[2]*=c;  
  v[3]*=c;  

  return *this;
}

inline GVector4D& GVector4D::operator/=(double c){

  double cinv=1.0/c;  
  v[0]*=cinv;  
  v[1]*=cinv;  
  v[2]*=cinv;  
  v[3]*=cinv;  

  return *this;
}

inline double operator*(const GVector4D& v1,const GVector4D& v2){

  return v1.v[0]*v2.v[0]-v1.v[1]*v2.v[1]-
         v1.v[2]*v2.v[2]-v1.v[3]*v2.v[3];
}

inline double GVector4D::cont(const GVector4D& v4) const {
  
  return v[0]*v4.v[0]-v[1]*v4.v[1]-
         v[2]*v4.v[2]-v[3]*v4.v[3];
}

inline GVector4D operator-(const GVector4D& v1,const GVector4D& v2){
  
  return GVector4D(v1)-=v2;
}

inline GVector4D operator+(const GVector4D& v1,const GVector4D& v2){
  
  return GVector4D(v1)+=v2;
}

inline double GVector4D::get(int i) const {
  return v[i];
}

inline void GVector4D::set(int i,double d){
  
  v[i]=d;
}

inline void GVector4D::set(double e,double p1,double p2, double p3){

  v[0]=e;
  v[1]=p1;
  v[2]=p2;
  v[3]=p3;
}

#endif

