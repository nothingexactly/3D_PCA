class Axes {
  Matrix vcs;
  
  Axes() {
    vcs = new Matrix(4,4);
    for (int i = 3 ; i < 3*4 ; i+=4) vcs.M[i] = 1;
  }
  
  void display(Matrix R, Matrix transforms, Matrix display) {
    
    Matrix vcsCopy = new Matrix(4,3);
    Matrix origin = new Matrix(new float[] {0,0,0,1},4,1);
    
    for (int i = 0 ; i < 12 ; i++) vcsCopy.M[i] = vcs.M[i];
    
    vcsCopy.mult(R).mult(transforms).project(2).mult(display);
    origin.mult(transforms).project(2).mult(display);
    textAlign(CENTER,CENTER);
    stroke(255,0,0);  // x
    strokeWeight(2);
    fill(255,0,0);
      line(origin.M[0],origin.M[1],vcsCopy.M[0],vcsCopy.M[1]);
      text("X",vcsCopy.M[0],vcsCopy.M[1]);
    stroke(0,0,255);  // z
    fill(0,0,255);
      line(origin.M[0],origin.M[1],vcsCopy.M[8],vcsCopy.M[9]);
      text("Z",vcsCopy.M[8],vcsCopy.M[9]);   
    stroke(0,255,0);  // y
    fill(0,255,0);
      line(origin.M[0],origin.M[1],vcsCopy.M[4],vcsCopy.M[5]);
      text("Y",vcsCopy.M[4],vcsCopy.M[5]);
   strokeWeight(1);
  }
  
} // end of class
