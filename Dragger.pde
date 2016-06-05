class Dragger {
  
  float[] dragStart;
  float[] dragEnd;
  Matrix screenVector;
  float dragMag;
  boolean pressed;
  boolean dragging;
  Matrix dRot;
  
  float magQuotient;
  float magUpper;
  float magLower;
  
  Dragger() {
    dragStart = new float[] {0,0};
    dragEnd = new float[] {0,0};
    dragging = pressed = false;
    dragMag = 0;
    
    magQuotient = 500;
    magUpper = 5;
    magLower = 0.05;
  }
  
  void pressed() {
      d.pressed = true;
  }
  
  void released() {
      d.dragging = d.pressed = false;
      d.dragMag = 0;
  }
  
  boolean check() {
   if (pressed && !dragging) {
      dragEnd = new float[] {pmouseX,pmouseY};
      dragging = true;
        return false;
   } else if (dragging) {
      dragStart = new float[] {dragEnd[0],dragEnd[1]};
      dragEnd = new float[] {mouseX,mouseY};
      screenVector = new Matrix(new float[] {dragEnd[1]-dragStart[1],dragEnd[0]-dragStart[0]},2,1);
      dragMag = pow(pow(screenVector.M[0],2)+pow(screenVector.M[1],2),0.5);
        return dragMag > 3 ? true : false;
   } else {
        return false; 
   }
  }
  
  Matrix reorientCamera(Matrix a) {
     float[] yAxis = new float[] {0,1,0};
     
     // first pitch about local X
     Quaternion localX_pitch = new Quaternion(screenVector.M[0]/magQuotient,a.isolate(1,1,3,1).M);
     Quaternion worldY_Pan = new Quaternion(screenVector.M[1]/magQuotient,yAxis);

     // check if local Y goes below world-XZ plane 
     Matrix xzNormal = new Matrix(yAxis,3,1);
     Matrix checkLocalY = a.isolate(1,2,3,1).mult(localX_pitch.getR(3));
     float d = xzNormal.dot(checkLocalY);
     
     if (abs(d)/d == -1) {  // safe to pitch
       return worldY_Pan.mult(localX_pitch).getR(4);
     } else {                    // neglect pitch
       return worldY_Pan.getR(4);
     } 
  }
  
} // end of class
