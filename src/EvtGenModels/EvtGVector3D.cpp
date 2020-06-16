
#include <iostream>
#include <math.h>
#include "EvtGenModels/EvtGVector3D.hh"
using std::ostream;

GVector3D::~GVector3D(){}

GVector3D::GVector3D(){
  
  v[0]=v[1]=v[2]=0.0;
}


GVector3D::GVector3D(double x,double y, double z){
  
  v[0]=x; v[1]=y; v[2]=z;
}

GVector3D rotateEuler(const GVector3D& v,
                        double alpha,double beta,double gamma){

  GVector3D tmp(v);
  tmp.applyRotateEuler(alpha,beta,gamma);
  return tmp;

}


void GVector3D::applyRotateEuler(double phi,double theta,double ksi){

  double temp[3];
  double sp,st,sk,cp,ct,ck;

  sp=sin(phi);
  st=sin(theta);
  sk=sin(ksi);
  cp=cos(phi);
  ct=cos(theta);
  ck=cos(ksi);

  temp[0]=( ck*ct*cp-sk*sp)*v[0]+( -sk*ct*cp-ck*sp)*v[1]+st*cp*v[2];
  temp[1]=( ck*ct*sp+sk*cp)*v[0]+(-sk*ct*sp+ck*cp)*v[1]+st*sp*v[2];
  temp[2]=-ck*st*v[0]+sk*st*v[1]+ct*v[2];


  v[0]=temp[0];
  v[1]=temp[1];
  v[2]=temp[2];
}

ostream& operator<<(ostream& s,const GVector3D& v){
 
  s<<"("<<v.v[0]<<","<<v.v[1]<<","<<v.v[2]<<")";

  return s;

}


GVector3D cross( const GVector3D& p1,const GVector3D& p2 ){

  return GVector3D(p1.v[1]*p2.v[2] - p1.v[2]*p2.v[1],
                     p1.v[2]*p2.v[0] - p1.v[0]*p2.v[2],
                     p1.v[0]*p2.v[1] - p1.v[1]*p2.v[0]);

}

double GVector3D::d3mag() const

// returns the 3 momentum mag.
{
  double temp;

  temp = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
  temp = sqrt( temp );

  return temp;
} // r3mag

double GVector3D::dot ( const GVector3D& p2 ){

  double temp;

  temp = v[0]*p2.v[0];
  temp += v[0]*p2.v[0];
  temp += v[0]*p2.v[0];
 
  return temp;
} //dot





