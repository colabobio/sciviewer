class Cell {
  String code;
  float umap2, umap1;
  float proj;
  boolean selected;
  
  float[] expression;
  
  Cell(String code, float umap1, float umap2) {
    this.code = code;
    this.umap1 = umap1;
    this.umap2 = umap2;
    
    this.selected = false;
  }
  
  void normalize(float min1, float max1, float min2, float max2) {
    umap1 = map(umap1, min1, max1, 0, 1);
    umap2 = map(umap2, min2, max2, 1, 0);
  }
  
  void initExpression(int numGenes) {
    expression = new float[numGenes];
  }
  
  void setExpression(int i, float level) {
    expression[i] = level;
  }
  
  void project(Selector sel) {
    if (selected) {
      PVector dirv = new PVector(sel.nx1 - sel.nx0, sel.ny1 - sel.ny0);
      PVector celv = new PVector(umap1 - sel.nx0, umap2 - sel.ny0);
      float a = PVector.angleBetween(dirv, celv);
      proj = cos(a) * celv.mag() / dirv.mag();
      //println(proj);
      //println(umap1, umap2, sel.nx0, sel.ny0, sel.nx1, sel.ny1);
    }    
  }
  
  void display(float x0, float y0, float w, float h) {
    float x = map(umap1, 0, 1, x0, x0 + w);
    float y = map(umap2, 0, 1, y0, y0 + h);
    
    noStroke();
    if (selGene != -1) {
      float f = constrain(map(expression[selGene], minGeneExp, maxGeneExp, 0, 1), 0, 1);
      colorMode(HSB, 360, 100, 100);
      fill((1 - f) * 170 + f * 233, 74, 93, 80);
      colorMode(RGB, 255, 255, 255);
    } else {
      if (selected) {
        fill(240, 118, 104, 80);      
      } else {
        fill(150, 80);
      }
    }
    ellipse(x, y, 5, 5);
  }  
}

class Gene {
  String name;
  int idx;
  Double r, rabs, p;
  
  Gene(String name, int idx, double r, double p) {
    this.name = name;
    this.idx = idx;
    this.r = r;
    this.rabs = Math.abs(r);
    this.p = p;
  }
}

class Selector {
  int state;
  
  int spx0, spy0;
  int spx1, spy1;
  int wx, wy;

  float nx0, ny0;
  float nx1, ny1;

  float angle;  
  float x0;
  float y0;
  float w, h;
  
  PMatrix2D tmat;
    
  Selector() {
    state = CLOSED;
    tmat = new PMatrix2D();
  }
  
  void display() {
    if (state == CLOSED) return;    
        
    if (state == SET_SPINE) {
      stroke(240, 118, 104);
      line(spx0, spy0, spx1, spy1);
    } else if (state == SET_WIDTH) {      
      stroke(240, 118, 104);
      line(spx0, spy0, spx1, spy1);
      line(spx1, spy1, wx, wy);
      displayBox();
    } else {
      displayBox();
    }
  }
  
  void apply(Cell cell, float sx0, float sy0, float sw, float sh) {
    float sx = map(cell.umap1, 0, 1, sx0, sx0 + sw);
    float sy = map(cell.umap2, 0, 1, sy0, sy0 + sh);
    
    float tx = tmat.multX(sx, sy);
    float ty = tmat.multY(sx, sy);
    
    //println(sx, sy, "|", tx, ty, "|", w, h);
    
    cell.selected = 0 <= tx && tx <= w && -h/2 <= ty && ty <= h/2;    
  }

  float transformY(float sx, float sy) {
    return tmat.multY(sx, sy);    
  }  
  
  void updateBox(int x, int y) {
    wx = x;
    wy = y;

    PVector spdir = new PVector(spx1 - spx0, spy1 - spy0);
    PVector bxdir = new PVector(wx - spx1, wy - spy1);    

    float a = PVector.angleBetween(spdir, bxdir);    
    float d = sin(a) * bxdir.mag();
    
    angle = (float)Math.atan2(spdir.y, spdir.x);
    
    h = 2 * d;
    w = spdir.mag(); 
    x0 = spx0;
    y0 = spy0;
    
    tmat.reset();
    tmat.rotate(-angle);
    tmat.translate(-x0, -y0);
    
    //println(degrees(angle));
  }
  
  void normalize(float sx0, float sy0, float sw, float sh) {
    nx0 = map(spx0, sx0, sx0 + sw, 0, 1);  
    nx1 = map(spx1, sx0, sx0 + sw, 0, 1);
    ny0 = map(spy0, sy0, sy0 + sh, 0, 1);
    ny1 = map(spy1, sy0, sy0 + sh, 0, 1);
  }
  
  void displayBox() {
    stroke(240, 118, 104);
    noFill();
    pushMatrix();
    translate(x0, y0);
    rotate(angle);
    rect(0, -h/2, w, h);
    popMatrix();    
  }  
  
  void press(int x, int y) {
    if (state == CLOSED || state == COMPLETED) {
      state = SET_SPINE;
    } else if (state == SET_SPINE) {
      state = SET_WIDTH;
    }

    if (state == SET_SPINE) {
      spx0 = spx1 = wx = x;
      spy0 = spy1 = wy = y;
      angle = 0;
      w = h = 0;
    }
  }
  
  void drag(int x, int y) {
    if (state == SET_SPINE) {
      spx1 = x;
      spy1 = y;      
    }    
  }  

  void move(int x, int y) {
    if (state == SET_WIDTH) {
      updateBox(x, y);
    }
  }

  void release(int x, int y) {
    if (state == SET_SPINE) {
      spx1 = x;
      spy1 = y;
      if (spx1 != spx0 || spy1 != spy0) {
        state = SET_WIDTH;
        wx = spx1;
        wy = spy1;
      } else {
        state = CLOSED;
      }
    } else if (state == SET_WIDTH) {
      updateBox(x, y);
      state = COMPLETED;
      requestSelection = true;
    }
  }
}
