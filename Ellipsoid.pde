// This four-step method is given by Lengyel (pp. 218-220)
// 1. Perform PCA on vertices
// 2. Use Principal Components to scale the vertices so that 
//    they are bounded by a cube
// 3. Find the bounding sphere
// 4. stretch the sphere into an ellipsoid by transforming back
//    into unscaled coordinate space

class Ellipsoid {
  Matrix pc1ExtentVecs;
  Matrix Q; // centre
  float radius;
  Matrix semiAxes;
  Matrix ranges;
  Matrix scaleMat;
  Matrix pcs;
  Matrix radVec_1;
  
  Ellipsoid(Matrix valsVecs, Matrix vcs) {
    
    Matrix pcsExtents = new Matrix(new float[] {0,0,0,0,0,0},3,2); 
      /*   | min_dot_pc1 max_dot_pc1 |
           | min_dot_pc2 max_dot_pc2 |
           | min_dot_pc3 max_dot_pc3 |    */
    
    pc1ExtentVecs = new Matrix(3,2);
      /*   | pc1_low_X pc1_high_X |
           | pc1_low_Y pc1_high_Y |
           | pc1_low_Z pc1_high_Z |    */
    
    
    float[] lambda = valsVecs.isolate(1,4,3,1).M;
             pcs = valsVecs.isolate(1,1,3,3);

    float dp = 0;
     
    // establish initial extent values
    // for each principal component
    for (int k = 0 ; k < pcs.m ; k++) {
       // just project first vertex onto each principal component
       Matrix point = vcs.isolate(1,1,3,1);
              dp = point.dot(pcs.isolate(1,k+1,3,1));
         pcsExtents.M[k] = dp; 
         pcsExtents.M[k+3] = dp; 
    }
    
    // find the lengths of the extents in the 
    // principal component directions
    
    // the aim is to repopulate pcsExtents
    // and save the vertices producing the greatest
    // extents along principal component 1
    
    for (int k = 0 ; k < pcs.m ; k++) {
      Matrix prinComponent = pcs.isolate(1,k+1,3,1);
      Boolean firstComponent = k == 0 ? true : false;
      
      for (int j = 1 ; j <= vcs.m ; j++) {
         Matrix point = vcs.isolate(1,j,3,1);
                dp = point.dot(prinComponent);

         if (dp < pcsExtents.M[k]) {
           pcsExtents.M[k] = dp;    // lower
             if (firstComponent) { 
               for (int i = 0 ; i < 3 ; i++) pc1ExtentVecs.M[i] = point.M[i]; 
             }    // min in pc1 direction
         } else if (dp > pcsExtents.M[k+3]) {
           pcsExtents.M[k+3] = dp;  // upper
             if (firstComponent) { 
               for (int i = 0 ; i < 3 ; i++) pc1ExtentVecs.M[3+i] = point.M[i]; 
             }  // max in pc1 direction
         }
         
       }
    }
    
    // use these ranges as scale factors
    ranges = new Matrix(3,3);                          
      ranges.M[0] = pcsExtents.M[3]-pcsExtents.M[0];  
      ranges.M[4] = pcsExtents.M[4]-pcsExtents.M[1];
      ranges.M[8] = pcsExtents.M[5]-pcsExtents.M[2];
    
    // construct a matrix that applies non-uniform scale to
    // vertices so they can be contained in a bounding cube
    
    // Equation 8.21
    scaleMat = pcs.getTranspose();  // pcs orthogonal, so transpose == inverse
    scaleMat.mult(ranges.getInverse()).mult(pcs);
  
    // now proceed as though finding bounding sphere
    Q = new Matrix(3,1);
    
    // to start, use the first two precalculated extents along pc1
    
    // transform into the cube
    pc1ExtentVecs.mult(scaleMat);
    
    Q = Q.matrixAdd(pc1ExtentVecs.isolate(1,2,3,1),pc1ExtentVecs.isolate(1,1,3,1));
    
    Q.scale(0.5); // find the midpoint as the initial centre
    
    // construct the radius as an extent minus Q, the centre
           radVec_1 = Q.matrixSub(pc1ExtentVecs.isolate(1,1,3,1),Q);
    float radiusSqd_1 = radVec_1.dot(radVec_1);  // find the mag_sqd of the radius
    
    // surface point
    Matrix G = new Matrix(3,1);

    // test if all points are enclosed in the sphere
    for (int j = 0 ; j < vcs.m ; j++) {

       Matrix P_j = vcs.isolate(1,j+1,3,1);
              P_j.mult(scaleMat);            // transform to the cube
       Matrix radVec_2 = G.matrixSub(P_j,Q); // take difference from the centre Q

       float radiusSqd_2 = radVec_2.dot(radVec_2);
       
       if (radiusSqd_2 > radiusSqd_1) {      // compare
          // then find G
          float sc = pow(radiusSqd_1,0.5)/pow(radiusSqd_2,0.5);
          radVec_2.scale(sc);      
          
          G = G.matrixSub(Q,radVec_2);
          
          // and use G to redefine Q, the centre
          Q = G.matrixAdd(G,P_j);
          Q.scale(0.5);
                    
          // then update the radius vector
          radVec_1 = Q.matrixSub(P_j,Q);
          radiusSqd_1 = radVec_1.dot(radVec_1);
       }
    }

    // transform Q, the centre, back into unscaled coordinate 
    // space to give the ellipsoid centre
    
    Matrix sinv = scaleMat.getInverse();

    Q.mult(sinv);

    // define semi-axis lengths
    float finalRadius = pow(radiusSqd_1,0.5);
    semiAxes = new Matrix(new float[] {finalRadius*ranges.M[0],
                                        finalRadius*ranges.M[4],
                                        finalRadius*ranges.M[8]},3,1);
  }
  
  void display(Matrix transforms,Matrix display) {

    Matrix calculatedQ = new Matrix(new float[] {Q.M[0],Q.M[1],Q.M[2],1},4,1);
           calculatedQ.mult(transforms).project(2).mult(display);
    
      textAlign(CENTER,CENTER);
      stroke(#52AAF5);
      noFill();
      ellipse(calculatedQ.M[0],calculatedQ.M[1],20,20);
      fill(255);
      text("cQ",calculatedQ.M[0],calculatedQ.M[1]);
      
    // -1 to +1 
    // -0.5 to +0.5
    // -0.3 to +0.3
    
    float theta = 0;
    float phi = 0;
    float full = 2*PI;
    int steps = 20;
    float angInc = full/(steps);
    float sc1,sc2,sc3;
    float radius = pow(radVec_1.dot(radVec_1),0.5);
    
    Matrix pc1_h = new Matrix(new float[] {pcs.M[0],pcs.M[1],pcs.M[2],1},4,1);
    Matrix pc2_h = new Matrix(new float[] {pcs.M[3],pcs.M[4],pcs.M[5],1},4,1);
    Matrix pc3_h = new Matrix(new float[] {pcs.M[6],pcs.M[7],pcs.M[8],1},4,1);

    calculatedQ = new Matrix(new float[] {Q.M[0],Q.M[1],Q.M[2],1},4,1);
        
    stroke(255,50);
    noFill();
    
    for (int k = 0 ; k < steps ; k++) {
       phi = 0;
       theta += angInc;
       Matrix surfacePoint;
       beginShape();
       for (int j = 0 ; j < steps ; j++) {
          phi += angInc;
          sc1 = radius*ranges.M[0]*cos(theta)*sin(phi);
          sc2 = radius*ranges.M[4]*sin(theta)*sin(phi);
          sc3 = radius*ranges.M[8]*cos(phi);
          
          surfacePoint = pc1_h.getCopy();
          Matrix pc2_Copy = pc2_h.getCopy();
          Matrix pc3_Copy = pc3_h.getCopy();
          
          pc2_Copy.scale(sc2);
          pc3_Copy.scale(sc3);
          
          surfacePoint.scale(sc1);
          surfacePoint.add(pc2_Copy).add(pc3_Copy).add(calculatedQ);
          
          surfacePoint.M[3] = 1;
          surfacePoint.mult(transforms).project(2).mult(display);
            vertex(surfacePoint.M[0],surfacePoint.M[1]);
       }
       endShape(CLOSE);
    }

      for (int k = 0 ; k < steps ; k++) {
       phi += angInc;
       theta = 0;
       Matrix surfacePoint = new Matrix(4,1);
       beginShape();
       for (int j = 0 ; j < steps ; j++) {
          theta += angInc;
          sc1 = radius*ranges.M[0]*cos(theta)*sin(phi);
          sc2 = radius*ranges.M[4]*sin(theta)*sin(phi);
          sc3 = radius*ranges.M[8]*cos(phi);
          
          surfacePoint = pc1_h.getCopy();
          Matrix pc2_Copy = pc2_h.getCopy();
          Matrix pc3_Copy = pc3_h.getCopy();
          
          pc2_Copy.scale(sc2);
          pc3_Copy.scale(sc3);
          
          surfacePoint.scale(sc1);
          surfacePoint.add(pc2_Copy).add(pc3_Copy).add(calculatedQ);
          surfacePoint.M[3] = 1;
          surfacePoint.mult(transforms).project(2).mult(display);
            vertex(surfacePoint.M[0],surfacePoint.M[1]);
       }
       endShape(CLOSE);
    }

  }
  
} // end of class
