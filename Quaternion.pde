// this class isn't yet used in the most
// efficient way. Over-reliance on .getR()

// See Equations 4.32 to 4.69 in Lengyel (pp. 80-89)

class Quaternion {
  
  float s;
  float[] v;
  
  /*
    CONSTRUCTORS
  */
  
  // find look-at rotation
  Quaternion(float[] _from, float[] _to) {
     v = new float[3];
     
     float[] from = new float[3];
     float[] to = new float[3];
     
     // copy and normalize
     
     float magF = pow(this.dot(_from,_from),0.5);
     float magT = pow(this.dot(_to,_to),0.5);
     
     for (int i = 0 ; i < 3 ; i++) { 
         from[i] = _from[i]/magF;
         to[i] = _to[i]/magT;
     }
     
     this.axisAngleInit(acos(this.dot(from,to)),this.cross(from,to));
  }
  
  Quaternion() {
    s = 0;
    v = new float[] {0,1,0};
  }
  
  Quaternion(float[] _v) {
    s = 0;
    v = new float[] {_v[0],_v[1],_v[2]};
  }  
  
  // There is an opportunity for significat improvement of these constructors. See:
  // http://lolengine.net/blog/2013/09/18/beautiful-maths-quaternion-from-vectors
  
  // unit quaternion, rotates 'theta' about 'axis'
  Quaternion(float theta, float[] axis) {
    v = new float[3];
    this.axisAngleInit(theta,axis);
  }
  
  Quaternion(float theta, float p1, float p2, float p3) {
    v = new float[3];
    this.axisAngleInit(theta,new float[] {p1,p2,p3});
  }
  
  // rotation matrix to quaternion
  Quaternion(Matrix A) {
      s = 0;
      this.v = new float[] {0,1,0};
      // r11+r22+r33 == 4s^2 - 1
      // r11-r22-r33 == 4x^2 - 1
      // -r11+r22-r33 == 4y^2 - 1
      // -r11-r22+r33 == 4z^2 - 1
    if (A.M.length == 16) {
      s = pow(A.M[0]+A.M[5]+A.M[10]+1,0.5)/2;
      this.v[0] = pow(A.M[0]-A.M[5]-A.M[10]+1,0.5)/2;
      this.v[1] = pow(-A.M[0]+A.M[5]-A.M[10]+1,0.5)/2;
      this.v[2] = pow(-A.M[0]-A.M[5]+A.M[10]+1,0.5)/2;
    } else if (A.M.length == 9) {
      this.s = pow(A.M[0]+A.M[4]+A.M[8]+1,0.5)/2;
      this.v[0] = pow(A.M[0]-A.M[4]-A.M[8]+1,0.5)/2;
      this.v[1] = pow(-A.M[0]+A.M[4]-A.M[8]+1,0.5)/2;
      this.v[2] = pow(-A.M[0]-A.M[4]+A.M[8]+1,0.5)/2;
    }
  }
    
  void axisAngleInit(float theta, float[] axis) {
    this.s = cos(theta/2);
    
    float mag = 0;

    for (int i = 0 ; i < 3 ; i++) v[i] = axis[i];
    mag = pow(this.dot(axis,axis),0.5);
    
    for (int i = 0 ; i < 3 ; i++) { 
      this.v[i] /= mag;  // normalize
      this.v[i] *= sin(theta/2);
    }
  }
  
  // this*q: the final term in the multiplication
  // corresponds to the first of the rotations
  Quaternion mult(Quaternion q) {
     float _s = this.s*q.s-this.dot(this.v,q.v);  // s1*s2-v1.dot(v2)
     float[] v1Xv2 = this.cross(this.v,q.v);      // v1.cross(v2)
     
     // v1.cross(v2) + s1*v2 + s2*v1
     for (int i = 0 ; i < 3 ; i++) this.v[i] = v1Xv2[i]+this.s*q.v[i]+q.s*this.v[i];
     this.s = _s;
       return this; 
  }
  
  // P' = q*P*qinv
  Quaternion rotate(Quaternion q) {
    float magSqrd = q.s*q.s+this.dot(q.v,q.v);
    
    Quaternion qinv = new Quaternion(new float[] {-q.v[0]/magSqrd,-q.v[1]/magSqrd,-q.v[2]/magSqrd});
               qinv.s = q.s/magSqrd; // shame about this step
               
    Quaternion qcopy = new Quaternion(new float[] {q.v[0],q.v[1],q.v[2]});
               qcopy.s = q.s;        // shame about this step
    
    qcopy.mult(this).mult(qinv); // q*P*qinv
     
     // thisRotated = q*P*qinv
     this.s = qcopy.s;
     for (int i = 0 ; i < 3 ; i++) this.v[i] = qcopy.v[i];
  
      return this;
  }
  
  float[] cross(float[] a, float[] b) {
      return new float[] {a[1]*b[2]-a[2]*b[1],
                           a[2]*b[0]-a[0]*b[2],
                           a[0]*b[1]-a[1]*b[0]};
  }
  
  float dot(float[] a, float[] b) {
      return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
  }
  
  
Matrix getR(int dim) {
      float[] R;
      float[] RDim;
       
        // 3x3
        R = new float[9];
        R[0] = 1-2*this.v[1]*this.v[1]-2*this.v[2]*this.v[2];
        R[1] = 2*this.v[0]*this.v[1]+2*this.s*this.v[2];
        R[2] = 2*this.v[0]*this.v[2]-2*this.s*this.v[1];
        
        R[3] = 2*this.v[0]*this.v[1]-2*this.s*this.v[2];
        R[4] = 1-2*this.v[0]*this.v[0]-2*this.v[2]*this.v[2];
        R[5] = 2*this.v[1]*this.v[2]+2*this.s*this.v[0];
        
        R[6] = 2*this.v[0]*this.v[2]+2*this.s*this.v[1];
        R[7] = 2*this.v[1]*this.v[2]-2*this.s*this.v[0];
        R[8] = 1-2*this.v[0]*this.v[0]-2*this.v[1]*this.v[1];
        
        if (dim == 3) return new Matrix(R,3,3);
        
        // 4x4
        RDim = new float[4*4];
        for (int i = 0 ; i < 4 ; i++) {
           for (int j = 0; j < 4 ; j++) {
            RDim[j*4+i] = i < 3 && j < 3 ? R[j*3+i] : i == j ? 1 : 0;
           }
        }
        
        return new Matrix(RDim,4,4);
  }
  
  // interpolate t between q1 and q2
  Quaternion interpolate(float t, Quaternion q1, Quaternion q2) {
    
    // 'signs of the quaternions q1 and q2 are usually chosen
    // such that q1 dot q2 >= 0.' (p. 89)
    
    float dp = q1.s*q2.s+this.dot(q1.v,q2.v);
    float dpB = -q1.s*q2.s-q1.v[0]*q2.v[0]-q1.v[1]*q2.v[1]-q1.v[2]*q2.v[2];
    
           dp = dpB > dp ? dpB : dp;
    
    float theta = acos(dp);
    float sinTheta = pow(1-dp*dp,2);
    
    float term1 = sin(theta*(1-t))/sin(theta);
    float term2 = sin(theta*t)/sin(theta);
    
    // q(t) = term1*q1 + term2*q2
    
     Quaternion rval = new Quaternion(new float[] {
                          term1*q1.v[0]+term2*q2.v[0],                                       
                          term1*q1.v[1]+term2*q2.v[1],  
                          term1*q1.v[2]+term2*q2.v[2]});
                          
                rval.s = term1*q1.s+term2*q2.s;
          return rval;
  }
  
  Quaternion interpolate(float t, Quaternion q2) {
        return this.interpolate(t,this,q2);
  }
  
} // end of class
