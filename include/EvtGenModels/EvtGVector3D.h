#ifndef GVECTOR3D_HH
#define GVECTOR3D_HH

#include <iosfwd>

class GVector3D {

  friend GVector3D rotateEuler(const GVector3D& v,
                                 double phi,double theta,double ksi);

  inline friend GVector3D operator*(double c,const GVector3D& v2); 
  inline friend double operator*(const GVector3D& v1,const GVector3D& v2); 
  inline friend GVector3D operator+(const GVector3D& v1,const GVector3D& v2);
  inline friend GVector3D operator-(const GVector3D& v1,const GVector3D& v2);
  inline friend GVector3D operator*(const GVector3D& v1,double c);
  inline friend GVector3D operator/(const GVector3D& v1,double c);
  friend GVector3D cross(const GVector3D& v1,const GVector3D& v2);
  
public:
  GVector3D();
  GVector3D(double x,double y ,double z);
  virtual ~GVector3D(); 
  inline GVector3D& operator*=(const double c);
  inline GVector3D& operator/=(const double c);
  inline GVector3D& operator+=(const GVector3D& v2);
  inline GVector3D& operator-=(const GVector3D& v2);
  inline void set(int i,double d);
  inline void set(double x,double y ,double z);
  void applyRotateEuler(double phi,double theta,double ksi);
  inline double get(int i) const;
  friend std::ostream& operator<<(std::ostream& s,const GVector3D& v);
  double dot(const GVector3D& v2);
  double d3mag() const;

private:

  double v[3];

};

inline GVector3D& GVector3D::operator*=(const double c){

  v[0]*=c;
  v[1]*=c;
  v[2]*=c;
  return *this;
}

inline GVector3D& GVector3D::operator/=(const double c){

  v[0]/=c;
  v[1]/=c;
  v[2]/=c;
  return *this;
}

inline GVector3D& GVector3D::operator+=(const GVector3D& v2){

  v[0]+=v2.v[0];
  v[1]+=v2.v[1];
  v[2]+=v2.v[2];
  return *this;
}

inline GVector3D& GVector3D::operator-=(const GVector3D& v2){

  v[0]-=v2.v[0];
  v[1]-=v2.v[1];
  v[2]-=v2.v[2];
  return *this;
}

inline GVector3D operator*(double c,const GVector3D& v2){
  
  return GVector3D(v2)*=c;
}

inline GVector3D operator*(const GVector3D& v1,double c){
  
  return GVector3D(v1)*=c;
}

inline GVector3D operator/(const GVector3D& v1,double c){

  return GVector3D(v1)/=c; 
}

inline double operator*(const GVector3D& v1,const GVector3D& v2){

  return v1.v[0]*v2.v[0]+v1.v[1]*v2.v[1]+v1.v[2]*v2.v[2];
}

inline GVector3D operator+(const GVector3D& v1,const GVector3D& v2) {
  
  return GVector3D(v1)+=v2; 
}

inline GVector3D operator-(const GVector3D& v1,const GVector3D& v2) {
  
  return GVector3D(v1)-=v2; 

}

inline double GVector3D::get(int i) const {
  return v[i];
}

inline void GVector3D::set(int i,double d){
  
  v[i]=d;
}

inline void GVector3D::set(double x,double y, double z){

  v[0]=x;
  v[1]=y;
  v[2]=z;
}

#endif

