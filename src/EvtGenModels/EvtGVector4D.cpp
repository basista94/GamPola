
#include <iostream>
#include <cmath>
#include "EvtGenModels/EvtGVector4D.hh"
#include "EvtGenModels/EvtGVector3D.hh"

using std::ostream;

GVector4D::GVector4D() {
  v[0] = 0.0; v[1] = 0.0; v[2] = 0.0; v[3] = 0.0;
}

GVector4D::GVector4D(double e,double p1,double p2, double p3){
  
  v[0]=e; v[1]=p1; v[2]=p2; v[3]=p3;
}

double GVector4D::mass() const{

  double m2=v[0]*v[0]-v[1]*v[1]-v[2]*v[2]-v[3]*v[3];

  if (m2>0.0) {
    return sqrt(m2);
  }
  else{
    return 0.0;
  }
}


GVector4D rotateEuler(const GVector4D& rs,
                        double alpha,double beta,double gamma){

  GVector4D tmp(rs);
  tmp.applyRotateEuler(alpha,beta,gamma);
  return tmp;

}

GVector4D boostTo(const GVector4D& rs,
                    const GVector4D& p4, bool inverse){

  GVector4D tmp(rs);
  tmp.applyBoostTo(p4, inverse);
  return tmp;

}

GVector4D boostTo(const GVector4D& rs,
                    const GVector3D& boost, bool inverse){

  GVector4D tmp(rs);
  tmp.applyBoostTo(boost, inverse);
  return tmp;

}



void GVector4D::applyRotateEuler(double phi,double theta,double ksi){

  double sp=sin(phi);
  double st=sin(theta);
  double sk=sin(ksi);
  double cp=cos(phi);
  double ct=cos(theta);
  double ck=cos(ksi);

  double x=( ck*ct*cp-sk*sp)*v[1]+( -sk*ct*cp-ck*sp)*v[2]+st*cp*v[3];
  double y=( ck*ct*sp+sk*cp)*v[1]+(-sk*ct*sp+ck*cp)*v[2]+st*sp*v[3];
  double z=-ck*st*v[1]+sk*st*v[2]+ct*v[3];

  v[1]=x;
  v[2]=y;
  v[3]=z;
  
}

ostream& operator<<(ostream& s, const GVector4D& v){

  s<<"("<<v.v[0]<<","<<v.v[1]<<","<<v.v[2]<<","<<v.v[3]<<")";

  return s;

}

void GVector4D::applyBoostTo(const GVector4D& p4, bool inverse){

  double e=p4.get(0);

  GVector3D boost(p4.get(1)/e,p4.get(2)/e,p4.get(3)/e);

  applyBoostTo(boost, inverse);

  return;

}

void GVector4D::applyBoostTo(const GVector3D& boost, bool inverse){

  double bx,by,bz,gamma,b2;

  bx=boost.get(0);
  by=boost.get(1);
  bz=boost.get(2);

  double bxx=bx*bx;
  double byy=by*by;
  double bzz=bz*bz;

  b2=bxx+byy+bzz;

  if (b2 > 0.0 && b2 < 1.0) {

    gamma=1.0/sqrt(1.0-b2);

    double gb2=(gamma-1.0)/b2;

    double gb2xy=gb2*bx*by;
    double gb2xz=gb2*bx*bz;
    double gb2yz=gb2*by*bz;

    double gbx=gamma*bx;
    double gby=gamma*by;
    double gbz=gamma*bz;

    double e2=v[0];
    double px2=v[1];
    double py2=v[2];
    double pz2=v[3];

    if ( inverse ) {
      v[0]=gamma*e2-gbx*px2-gby*py2-gbz*pz2;
      
      v[1]=-gbx*e2+gb2*bxx*px2+px2+gb2xy*py2+gb2xz*pz2;
 
      v[2]=-gby*e2+gb2*byy*py2+py2+gb2xy*px2+gb2yz*pz2;
 
      v[3]=-gbz*e2+gb2*bzz*pz2+pz2+gb2yz*py2+gb2xz*px2;
    }
    else {
      v[0]=gamma*e2+gbx*px2+gby*py2+gbz*pz2;
      
      v[1]=gbx*e2+gb2*bxx*px2+px2+gb2xy*py2+gb2xz*pz2;
 
      v[2]=gby*e2+gb2*byy*py2+py2+gb2xy*px2+gb2yz*pz2;
 
      v[3]=gbz*e2+gb2*bzz*pz2+pz2+gb2yz*py2+gb2xz*px2;
    }
  }

}

GVector4D GVector4D::cross( const GVector4D& p2 ){

  //Calcs the cross product.  Added by djl on July 27, 1995.
  //Modified for real vectros by ryd Aug 28-96

  GVector4D temp;
  
  temp.v[0] = 0.0; 
  temp.v[1] = v[2]*p2.v[3] - v[3]*p2.v[2];
  temp.v[2] = v[3]*p2.v[1] - v[1]*p2.v[3];
  temp.v[3] = v[1]*p2.v[2] - v[2]*p2.v[1];

  return temp;
}

double GVector4D::d3mag() const

// returns the 3 momentum mag.
{
  double temp;

  temp = v[1]*v[1]+v[2]*v[2]+v[3]*v[3];

  temp = sqrt( temp );

  return temp;
} // r3mag

double GVector4D::dot ( const GVector4D& p2 )const{

  //Returns the dot product of the 3 momentum.  Added by
  //djl on July 27, 1995.  for real!!!

  double temp;

  temp = v[1]*p2.v[1];
  temp += v[2]*p2.v[2];
  temp += v[3]*p2.v[3];
 
  return temp;

} //dot

// Calculate the 3-d dot product of 4-vectors p1 and p2 in the rest frame of
// 4-vector p0
double GVector4D::dotr3( const GVector4D& p1, const GVector4D& p2 ) const
{
    return 1/mass2() * ((*this) * p1) * ((*this) * p2) - p1 * p2;
}

// Calculate the 3-d magnitude squared of 4-vector p1 in the rest frame of
// 4-vector p0
double GVector4D::mag2r3( const GVector4D& p1 ) const
{
    return Square((*this) * p1)/mass2() - p1.mass2();
}

// Calculate the 3-d magnitude 4-vector p1 in the rest frame of 4-vector p0.
double GVector4D::magr3( const GVector4D& p1 ) const
{
    return sqrt(mag2r3(p1));
}







