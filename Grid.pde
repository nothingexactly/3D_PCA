class Grid {
  Matrix Line;
  Matrix l2;
  
  Grid() {
    
  }
  
  void display(Matrix transforms,Matrix display,int intensity) {
     float span = 2;
     float nlines = 10;
     float tick;
     Line = new Matrix(4,10);
     Matrix Alt;
     stroke(255,map(intensity,0,60,0,120));
     
     for (int i = 0 ; i <= nlines ; i++) {
       
       // Note translation -1 along y axis
       Alt = new Matrix(new float[] {0,-1,0,1,
                                      0,-1,0,1,
                                      0,-1,0,1,
                                      0,-1,0,1},4,4);

        // scaled by a factor of 2
        tick = -span+i*(2*span/nlines);
        
        Alt.M[0] = -span;
        Alt.M[2] = tick;
        
        Alt.M[4] = span;
        Alt.M[6] = tick;
        
        Alt.M[8] = tick;
        Alt.M[10] = -span;
        
        Alt.M[12] = tick;
        Alt.M[14] = span;
        
        Alt.mult(transforms).project(2).mult(display);
        line(Alt.M[0],Alt.M[1],Alt.M[4], Alt.M[5]);
        line(Alt.M[8],Alt.M[9],Alt.M[12], Alt.M[13]); 
     }
  }
  
} // end of class
