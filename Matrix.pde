// column major

class Matrix {

  int n;  // rows
  int m;  // cols
  float[] M;  // col major
  
  // identity matrix
  Matrix(int _n, int _m) {
    this.n = _n;
    this.m = _m;
    this.M = new float[n*m];
    for(int i = 0 ; i < n*m ; i++) this.M[i] = (i+this.n)%this.n == i/this.n ? 1 : 0;
  }
   
  Matrix(float[] a, int _n, int _m) {
    this.n = _n;
    this.m = _m;
    this.M = new float[n*m];
    for (int i = 0 ; i < n*m ; i++) this.M[i] = a[i];
  }
  
  // rotation matrix
  Matrix(float angle, Matrix axis, int dim) {
   this.n = dim < 3 || dim > 4 ? 3 : dim;
   this.m = n;
   
   this.M = new float[this.m*this.n];
   
   // Equation 4.20 in Lengyel, for rotation of P about A by theta:
   // P' = PCos(theta) + (AxP)Sin(theta) + A(A.P)(1-Cos(theta))
   
   // term 1 - identity matrix
   Matrix A = new Matrix(3,3);
   
   // term 2 - expression of the cross product as a linear transformation
   Matrix B = new Matrix(new float[] {0,axis.M[2],-axis.M[1],
                                       -axis.M[2],0,axis.M[0],
                                       axis.M[1],-axis.M[0],0},3,3);
                                       

   // term 3 - outer product of axis with itself
   Matrix C = new Matrix(3,3);
      
   for (int j = 0 ; j < 3 ; j++) {
    for (int i = 0 ; i < 3 ; i++) C.M[j*C.n+i] = axis.M[i]*axis.M[j];
   }
   
   A.scale(cos(angle));
   B.scale(sin(angle));
   C.scale(1-cos(angle));
   
   Matrix D = matrixAdd(A,B);
          D = matrixAdd(D,C);
   
    if (dim == 3) { 
        for (int i = 0 ; i < 9 ; i++) {
            this.M[i] = D.M[i];
        }
    } else {
        for (int i = 0 ; i < 4 ; i++) {
           for (int j = 0; j < 4 ; j++) {
            M[j*4+i] = i < 3 && j < 3 ? D.M[j*3+i] : i == j ? 1 : 0;
           }
        }   
    }
}
 
 // transforms This with A
 Matrix mult(Matrix A) {
   float[] product;

   if (A.m != this.n) {
      println("inner dimensions don't agree");
   } else {
     product = new float[A.n*this.m];
     for (int k = 0 ; k < this.m ; k++) {
      for (int i = 0 ; i < A.n ; i++) {
       float dp = 0;
            for (int j = 0 ; j < this.n ; j++) dp += A.M[j*A.n+i]*this.M[k*this.n+j];
          product[k*A.n+i] = dp;
      }
     }
     
     this.M = new float[A.n*this.m];
     this.M = product;
   }
   
   return this;
 }
 
 // find a*b
 Matrix matrixProduct(Matrix a, Matrix b) {
   
   if (a.m != b.n) println("inner dimensions don't agree"); 
   
   Matrix rval = new Matrix(a.n,b.m);
   
   for (int k = 0 ; k < b.m ; k++) {
      for (int i = 0 ; i < a.n ; i++) {
       float dp = 0;
            for (int j = 0 ; j < a.m ; j++) dp += a.M[j*a.n+i]*b.M[k*b.n+j];
            rval.M[k*a.n+i] = dp;
      }
     }
     
     return rval;
 }
 
 // find a-b
 Matrix matrixSub(Matrix a, Matrix b) {
    Matrix rval = new Matrix(a.n < b.n ? a.n : b.n,
                             a.m < b.m ? a.m : b.m);
    
    for (int j = 0 ; j < rval.m ; j++) {
       for (int i = 0 ; i < rval.n ; i++) rval.M[j*rval.n+i] = a.M[j*a.n+i]-b.M[j*b.n+i];
    }
    
    return rval;
 }
 
 // find a+b
 Matrix matrixAdd(Matrix a, Matrix b) {
    Matrix rval = new Matrix(a.n < b.n ? a.n : b.n,
                             a.m < b.m ? a.m : b.m);
    
    for (int j = 0 ; j < rval.m ; j++) {
       for (int i = 0 ; i < rval.n ; i++) rval.M[j*rval.n+i] = a.M[j*a.n+i]+b.M[j*b.n+i];
    }
    
    return rval;
 }
 
 // this-a
 Matrix sub(Matrix a) {
    int upperI = a.n < this.n ? a.n : this.n;
    int upperJ = a.m < this.m ? a.m : this.m;
    for (int j = 0 ; j < upperJ ; j++) {
       for (int i = 0 ; i < upperI ; i++) this.M[j*this.n+i] -= a.M[j*a.n+i];
    }
    return this;
 }
 
 Matrix add(Matrix a) {
    int upperI = a.n < this.n ? a.n : this.n;
    int upperJ = a.m < this.m ? a.m : this.m;
    for (int j = 0 ; j < upperJ ; j++) {
       for (int i = 0 ; i < upperI ; i++) this.M[j*this.n+i] += a.M[j*a.n+i];
    }
    return this;
 }
 
 Matrix scale(float a) {
   for (int i = 0 ; i < this.M.length ; i++) this.M[i]*=a;
   return this;
 }
 
 Matrix cross(Matrix Q) {
   // only for 3D vectors
   if (Q.M.length == 3 && this.M.length == 3) {
    Matrix rval = new Matrix(new float[] {-this.M[2]*Q.M[1]+this.M[1]*Q.M[2],
                                          this.M[2]*Q.M[0]+-this.M[0]*Q.M[2],
                                          -this.M[1]*Q.M[0]+this.M[0]*Q.M[1]},3,1);
     return rval;
   }
   
   return this;
 }
 
 float dot(Matrix Q) {
   float dp = 0;
   for (int i = 0 ; i < this.M.length ; i++) dp += this.M[i]*Q.M[i];
   
   return dp;
 }
 
 // This outer Q
 Matrix outer(Matrix Q) {
   Matrix rval = new Matrix(this.n,Q.m);
   for (int j = 0 ; j < Q.m ; j++) {
    for (int i = 0 ; i < this.n ; i++) rval.M[j*this.n+i] = this.M[i]*Q.M[j];
   }
   return rval;
 }
 
 // only for 3-Vectors
 Matrix principalComponents(Matrix Pts) {
  // build covariance matrix
  float[] C = new float[this.n*this.n];
  
  float[] mean = new float[3];
  
  for (int j = 0 ; j < Pts.m ; j++) {
     for (int i = 0 ; i < 3 ; i++) mean[i] += Pts.M[j*3+i];
  }
  
  mean[0] /= Pts.m;
  mean[1] /= Pts.m;
  mean[2] /= Pts.m;
  
  // this, for me, is the more intuitive 
  // way to find the covariance matrix
  
  /* 
  Matrix datTrans = Pts.getTranspose();

  for (int j = 0 ; j < Pts.n ; j++) {
    for (int i = 0 ; i < Pts.m ; i++) {
      datTrans.M[j*datTrans.n+i] -= mean[j];
    }
  }
  
  Matrix dat = datTrans.getTranspose();
  
  Matrix Covar = matrixProduct(datTrans,dat);
  */
  
  float val = 0;
  
  // Equation 8.2 in Lengyel, p. 213
  for (int cj = 0 ; cj < 3 ; cj++) {
   for (int ci = 0 ; ci <= cj ; ci++) {  // ci<=cj because C is symmetric
      val = 0;
     for (int pj = 0 ; pj < Pts.m ; pj++) {
       val += (Pts.M[pj*3+ci]-mean[ci])*(Pts.M[pj*3+cj]-mean[cj]);
     }
     val /= Pts.m;
     C[cj*3+ci] = val;
     if (ci != cj) C[ci*3+cj] = val;    // because C is symmetric
   }
  }
  
 Matrix naturalAxes = new Matrix(3,4);
        naturalAxes = CalculateEigenSystem(C);
   return naturalAxes;
}
 
 // many assumptions here
 // would be much better if this func normalized
 // all of the column vectors- as though normalizing
 // vectors of a basis
 void normalize() {
    float mag = pow(this.dot(this),0.5);
    for (int i = 0 ; i < this.M.length ; i++) this.M[i] /= mag;
 }
 
 float mag() {
   return pow(this.dot(this),0.5);
 }
  
 Matrix project(int row) {
     for (int j = 0 ; j < m ; j++) {
      for (int i = 0 ; i <= row ; i++) M[j*n+i] = M[j*n+i]/M[j*n+row];
     }
   return this;
 }
 
 Matrix getTranspose() {
     Matrix rval = new Matrix(m,n);
     for (int j = 0 ; j < m ; j++) {
        for (int i = 0 ; i < n ; i++) rval.M[i*m+j] = M[j*n+i];
     }
  return rval;
 }
 
 Matrix getAverage() {
   Matrix rval = new Matrix(n,1);
   // find average
   for (int j = 0 ; j < m ; j++) {
      for (int i = 0 ; i < n ; i++) {
         rval.M[i] += M[j*n+i];
      } 
   }
   
   for (int i = 0 ; i < n ; i++) rval.M[i] /= m;
   
   return rval;
 }
 
 Matrix getInverse() {
   Matrix rval = new Matrix(n,n);
   if (n != m) {
      println("This is not a square matrix"); 
   } else {
     rval.M = invert(this.M);  
   }
   return rval;
 }
 
 Matrix getCopy() {
   Matrix rval = new Matrix(this.n,this.m);
   for (int i = 0 ; i < this.n*this.m ; i++) rval.M[i] = this.M[i];
   return rval;
 }
 
 // return rectangular portion of this matrix
 Matrix isolate(int iStart, int jStart, int iRange, int jRange) {
     Matrix rval = new Matrix(iRange,jRange);
     for (int j = 0 ; j < jRange ; j++) {
       for (int i = 0 ; i < iRange ; i++) {
          rval.M[j*iRange+i] = this.M[(jStart-1+j)*this.n+i]; 
       }
     }     
   return rval;
 }
 

// ------------------------------
// ------- MATRIX INVERSE -------
// ------------------------------

// for square matrices only
float[] invert(float[] _mat) {
  int n = (int) round(sqrt(_mat.length));
  float[] mat = new float[_mat.length];
  
  // copy
  for (int i = 0 ; i < mat.length ; i++) mat[i] = _mat[i];
  
  float[] identity = new float[mat.length];
  
  // build the identity matrix
  for (int j = 0 ; j < n ; j++) {
    for (int i = 0 ; i < n ; i++) identity[j*n+i] = i == j ? 1 : 0;
  }
  
  // Gauss-Jordan Elimination, Algorithm 3.13 in Lengyel, p. 42
  for (int j = 0 ; j < n ; j++) { // B. Set column j = 1
        int i = j;
        int sto = i;
        float stoVal = mat[j*n+i];
           
         while(i < n) {
              if (abs(mat[j*n+i]) > abs(stoVal)) {
                 stoVal = mat[j*n+i];
                 sto = i;                       
              }
              
              i++;
          }  // end of while
              
          if (stoVal == 0) {
            println("not invertible"); 
          } else {
            if (sto != j) {
               // pivot
               exchange(sto,j,mat);
               exchange(sto,j,identity);
             }
             
             stoVal = mat[j*n+j];
             
             weight(j,1/stoVal,mat);
             weight(j,1/stoVal,identity);
             
             for (int r = 0 ; r < n ; r++) {
                if (r != j) {
                  stoVal = -mat[j*n+r];

                  combine(j,stoVal,r,mat);
                  combine(j,stoVal,r,identity); 
                }
             }
           }
  } // end of loop
      
      return identity;
}  // end of invert()

//---------------------------------
//----Elementary Row Operations----
//---------------------------------

// add a multiple of one row 
// (w*r1) to another row (r2)
void combine(int r1,float w,int r2, float[] m) {
    int dim = (int) round(pow(m.length,0.5));
    float weighted = 0;
    
    for (int j = 0 ; j < dim ; j++) {
          weighted = m[j*dim+r1]*w;     // weight r1 element
          m[j*dim+r2] += weighted;      // add weighted r1 element to r2 element
    }
}

// multiply row r1 by a non-zero scalar
void weight(int r1, float w, float[] m) {
  int dim = (int) round(pow(m.length,0.5));
     for (int j = 0 ; j < dim ; j++) m[j*dim+r1] *= w; 
}

// exchange rows r1 and r2
void exchange(int r1, int r2, float[] m) {
    int dim = round(pow(m.length,0.5));
    float hold = 0;
    
    for (int j = 0 ; j < dim ; j++) {
        hold = m[j*dim+r1];          // save r1 element
          m[j*dim+r1] = m[j*dim+r2]; // replaced present r1 element with r2 element
          m[j*dim+r2] = hold;        // place saved r1 element where r2 element was
    }
}

//---------------------------------
//---------EIGEN-ANALYSIS----------
//---------------------------------

// This is effectively Code Listing 16.7 in Lengyel
// http://www.mathfor3dgameprogramming.com/code/Listing16.7.cpp
Matrix CalculateEigenSystem(float[] _m) {  // for a 3x3 symmetric matrix _m
  int maxSweeps = 32;
  float epsilon = 0.00004539992;  // check this is the val of 1e^-10
  
  float[] m = new float[9];
  float[] lambda = new float[3];   // becomes array of m's eigenvalues
  float[] r = new float[9];        // becomes matrix of m's eigenvectors (col major)
  Matrix rval = new Matrix(3,4);     // first three cols eigenvecs, final col eigenvals
  int n = 3;
  
  // copy entries of _m to m
  for (int i = 0 ; i < m.length ; i++) {
    m[i] = _m[i];
  }
  
  // set r to I
  for (int j = 0 ; j < n ; j++) {
    for (int i = 0 ; i < n ; i++) {
      r[j*n+i] = i == j ? 1 : 0; 
    }
  }
  
  // only need 6 vals because transformation 
  // M' = transpose(R)MR is symmetric also
  
  float m11 = m[0];
  float m12 = m[n];
  float m13 = m[2*n];
  float m22 = m[n+1];
  float m23 = m[2*n+1];
  float m33 = m[2*n+2];
  
  for (int a = 0 ; a < maxSweeps ; a++) {
    
    // Exit if off-diagonal entries small enough.
    if (((abs(m12) < epsilon) && (abs(m13) < epsilon)) && (abs(m23) < epsilon)) {
        println("System calculated in "+a+" sweeps");
          a = maxSweeps;
    }
    
    // Annihilate (1,2) entry.
    if (m12 != 0) {
      float u = (m22 - m11)*0.5/m12;  // Equation 16.46 in Lengyel. Derived from (16.43):
                                       // u == (M(p,p)-M(q,q)) / (2*M(p,q))
      float u2 = u*u;
      float u2p1 = u2+1;
      
      // t == sin(theta)/cos(theta), that is, t == tan(theta)
      
      // Lengyel sets tan(theta) == 1 when theta == 0, 
      // therefore shifting theta by PI/4
      
              // if u != inf ?             equation (16.49)              : equation (16.49.5)
      float t = (u2p1 != u2) ? ((u < 0 ? -1 : 1)*(pow(u2p1,0.5)-abs(u))) : 0.5/u;
      
      // find cos(theta) using t^2 + 1 = 1/cos(theta)^2
      float c = 1/pow(t*t+1,0.5); // cos(theta)
      float s = c*t;                   // sin(theta)
          
      // set main diagonal entries M'(p,p) and M'(q,q)
          
          m11 -= t*m12;  // equation (16.54)
          m22 += t*m12;  // equation (16.55)
      
      // set M'(p,q) to 0
      m12 = 0;
      
      float tempA = c*m13-s*m23;  // equation (16.40.2)
            m23 = s*m13+c*m23;     // equation (16.40.3)
            m13 = tempA;
      
      // operate on cols p and q of R
      // such that R == R_0 * R_2 * ... * R_M
      for (int i = 0 ; i < 3 ; i++) {
        float temp = c*r[i]-s*r[n+i];
              r[n+i] = s*r[i]+c*r[n+i];
              r[i] = temp;
      }
    }
    
    // Annihilate (1,3) entry.
    if (m13 != 0) {
      float u = (m33-m11)*0.5/m13;
      float u2 = u*u;
      float u2p1 = u2+1;
      float t = (u2p1 != u2) ? ((u < 0) ? -1 : 1) * (pow(u2p1,0.5) - abs(u)) : 0.5/u;
      float c = 1/pow(t*t+1,0.5);
      float s = c*t;
      
      m11 -= t*m13;
      m33 += t*m13;
      m13 = 0;
      
      float temp = c*m12-s*m23;
          m23 = s*m12+c*m23;
          m12 = temp;
      
      for (int i = 0 ; i < 3 ; i++) {
        float tempA = c*r[i]-s*r[2*n+i];
            r[2*n+i] = s*r[i]+c*r[2*n+i];
            r[i] = tempA;
      }
    }
    
    // Annihilate (2,3) entry.
    if (m23 != 0) {
      float u = (m33-m22)*0.5/m23;
      float u2 = u*u;
      float u2p1 = u2+1;
      float t = (u2p1 != u2) ? ((u < 0) ? -1 : 1)*(pow(u2p1,0.5) - abs(u)) : 0.5/u;
      float c = 1/pow(t*t+1,0.5);
      float s = c*t;
      
          m22 -= t*m23;
          m33 += t*m23;
          m23 = 0;
      
      float temp = c*m12-s*m13;
          m13 = s*m12+c*m13;
          m12 = temp;
      
      for (int i = 0 ; i < 3 ; i++) {
        float tempA = c*r[n+i]-s*r[2*n+i];
              r[2*n+i] = s*r[n+i]+c*r[2*n+i];
              r[n+i] = tempA;
      }
    }
  } // next sweep
  
  // eigenvalues
  lambda[0] = m11;
  lambda[1] = m22;
  lambda[2] = m33;
  
    // sort pcs by eigenvalues, high to low
    
    // assumes that eigenvecs are cols corresponding
    // to eigenval array positions
    while (!(lambda[0] > lambda[1] && lambda[1] > lambda[2])) {
      for (int j = 0 ; j < 2 ; j++) {
         if (lambda[j] < lambda[j+1]) {
            float hold;
            // swap principal components
            for (int i = 0 ; i < 3 ; i++) {
                hold = r[j*3+i];
                r[j*3+i] = r[(j+1)*3+i];
                r[(j+1)*3+i] = hold;
            }
            // swap eigenvalues
            hold = lambda[j+1];
            lambda[j+1] = lambda[j];
            lambda[j] = hold;
         }
      }
    }
  
  for (int i = 0 ; i < 9 ; i++) rval.M[i] = r[i];
  for (int i = 9 ; i < 12 ; i++) rval.M[i] = lambda[i-9];
  
    return rval;
}

} // end of class
