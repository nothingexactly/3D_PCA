class PositionTool {
  
  // This position tool assumes a very specific camera rig config.
  // The camera is attached to an arm whose pivot is *not* the centroid of
  // the data. It is the from the origin, pivotPoint, (that is, l1).
  // The position tool finds the rotation of the arm about the pivot point
  // which allows the camera to align it's lookAt with one of the principal
  // components.
  
  Matrix util;
  
  // camera's basis
  Matrix lookAt;
  Matrix frameI;
  Matrix frameJ;
  
  Matrix l1;  // camera arm pivot
  Matrix l2;  // principal component of data
  Matrix l3;  // vector from pivot to camera
  
  float l2_l1;  // theta between l2 and l1
  float l3_l2;  // theta between l3_adjust and l2
  float l3_l1;  // theta between l3_adjust and l1
  
  Quaternion camAlign;
  Quaternion armAlign;
  
  PositionTool(Matrix _lookAt,Matrix _frameI,Matrix _frameJ,
               Matrix component, Matrix pivotPoint, Matrix cameraRotationArm) {
     
     // --- 1. COPY PARAMS INTO OBJECT ---
     
     util = new Matrix(3,3);            
     
     // camera basis
     lookAt = new Matrix(new float[] {_lookAt.M[0],_lookAt.M[1],_lookAt.M[2]},3,1);
     lookAt.scale(-1);   // necessary to get the correct alignment
     
     frameI = new Matrix(new float[] {_frameI.M[0],_frameI.M[1],_frameI.M[2]},3,1);
     frameJ = new Matrix(new float[] {_frameJ.M[0],_frameJ.M[1],_frameJ.M[2]},3,1);
     
     l1 = new Matrix(new float[] {pivotPoint.M[0],pivotPoint.M[1],pivotPoint.M[2]},3,1);
     l2 = new Matrix(new float[] {component.M[0],component.M[1],component.M[2]},3,1);
     l3 = new Matrix(new float[] {cameraRotationArm.M[0],cameraRotationArm.M[1],cameraRotationArm.M[2]},3,1);
     
     // assumes orthogonal basis
     lookAt.normalize();
     frameI.normalize();
     frameJ.normalize();
    
     println("lookAt mag: "+lookAt.mag());
     println("frameJ mag: "+frameI.mag());
     println("frameK mag: "+frameJ.mag());
     println("");
       
     l2.normalize();
     
     // --- 2. FIND ROTATION FOR CAMERA'S PIVOT ARM ---
     // ---              l3 onto l3_adjust          ---
     
     // some sine rule magic to find the line
     // onto which the arm must be rotated
     
     l2_l1 = acos(l1.dot(l2)/l1.mag());
     
     l3_l2 = sin(l2_l1)*l1.mag()/l3.mag();
     l3_l2 = asin(l3_l2);
     
     l3_l1 = PI-l3_l2-l2_l1;
          
     float scale = l3.mag()*sin(l3_l1)/sin(l2_l1);
     
     Matrix intersection = l2.getCopy();  // i.e. intersection of camera origin with 
            intersection.normalize();     // an arbitrary line (principal component)
            intersection.scale(scale);    // from the centroid of the data
     
     Matrix l3_adjust = util.matrixSub(l1,intersection); // vector from pivot to intersection
     
     float adjustTheta = l3_l2;
          
    // axis about which camera arm should be rotated
    Matrix axisR = l3_adjust.cross(l3);
           axisR.normalize();
    
    // angle by which camera arm should rotate
    Matrix lla = l3.getCopy();
    Matrix llb = new Matrix(new float[] {-l3_adjust.M[0],-l3_adjust.M[1],-l3_adjust.M[2]},3,1);
      lla.normalize();
      llb.normalize();
    
    adjustTheta = acos(lla.dot(llb));
    
    // corresponding Quaternion for camera arm rotation
    armAlign = new Quaternion(adjustTheta,axisR.M);
    
    // --- 3. CHECK PROGRESS ---
    
      println("l2_l1 = "+l2_l1);
      println("l3_l2 = "+l3_l2);
      println("l3_l1 = "+l3_l1);
      println("");
      println("adjustTheta = "+adjustTheta);
      println("");
      println("l3 mag = "+l3.mag());
      println("l3 adjust mag = "+l3_adjust.mag());
     
    // --- 4. FIND ROTATION FOR CAMERA'S PIVOT ARM ---
    
     // align lookAt with component (works)
     
     Matrix lookAt_copy = lookAt.getCopy();
            lookAt_copy.normalize();
            
     Matrix l2_copy = l2.getCopy();
            l2_copy.normalize();
            
     float dAng = acos(lookAt_copy.dot(l2_copy));  // angle between lookAt and component
     Matrix rotAxis = lookAt_copy.cross(l2_copy);
     
     Quaternion rot_1 = new Quaternion(dAng,rotAxis.M);
     
     // bring frameJ into plane with normal vector targJ
     
     Matrix frameI_copy = frameI.getCopy();
            frameI_copy.normalize();
            frameI_copy.mult(rot_1.getR(3));
            
     Matrix targI = new Matrix(new float[] {0,1,0},3,1);
     Matrix planeIntersection = targI.cross(l2_copy);
            planeIntersection.normalize();
     
     // *** Find the cumulative rotation for the camera basis ***
     
     // Basis Rotation 1: align component and lookAt
     
     // angle by which basis should be rotated
     float angCorrect = acos(frameI_copy.dot(planeIntersection));
     
     // Basis Rotation 2: horizontal stabilisation of camera 
     
     // l2 is not always suitable axis of rotation
     // corAxis lies in the same line as l2 but is
     // sometimes opposite in direction. Hence need 
     // to go through standard procedure of using 
     // cross product to get this axis.
     Matrix corAxis = frameI_copy.cross(planeIntersection);
     
     // corresponding Quaternion for the two rotations to the basis
     
     camAlign = new Quaternion(angCorrect,corAxis.M); // second rotation
     camAlign.mult(rot_1);                            // cumulative rotation

  } // end of constructor

  Quaternion getCameraRotation() {
    return camAlign;
  }
 
  Quaternion getArmRotation() {
    return armAlign;
  }
  
  // transformation (as a rotation) of one
  // set of orthogonal basis vectors onto
  // another set of orthogonal basis vectors
  
  // this needs to be checked and tested
  Quaternion basisToBasis(Matrix _from, Matrix _to) {
     Matrix fromI = _from.isolate(1,1,3,1);
     Matrix fromJ = _from.isolate(1,2,3,1);
     Matrix toI = _to.isolate(1,1,3,1);
     Matrix toJ = _to.isolate(1,2,3,1);

     Quaternion r1 = new Quaternion(fromI.M,toI.M);
     
     fromJ.mult(r1.getR(3));
     Quaternion r2 = new Quaternion(fromJ.M,toJ.M);
     
     r2.mult(r1);  // the cumulative rotation
     
     return r2; 
  }
  
  Matrix getFullTransformation() {
   
    // l1 - camera arm pivot
    // l2 - principal component of data
    // l3 - vector from pivot to camera
   Matrix fullTransform = new Matrix(4,4);
   
   Matrix cTranslate = util.matrixAdd(util.matrixProduct(armAlign.getR(3),l3),l1);
  
   Matrix cRotate = camAlign.getR(3);

          
          for (int j = 0 ; j < 3 ; j++) {
           fullTransform.M[12+j] = cTranslate.M[j];
           for (int i = 0 ; i < 3 ; i++) fullTransform.M[j*4+i] = cRotate.M[j*3+i]; 
          }
          
   return fullTransform;
  }

} // end of class
