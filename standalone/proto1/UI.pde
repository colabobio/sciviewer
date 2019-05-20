float itemHeight = 50;
float itemSpace = 10;

void initUI() {
  selector = new Selector();
  indices = new ArrayList<Integer>();
  scrollList = new ScrollableList(width/2, 0, 200, height);
  
  float w = width - (width/2 + 200);
  exportBtn = new Button(width/2 + 200 + w/2 - 75, height - 75, 100, 30, "EXPORT");  
}

class ScrollableList {  
  ScrollBar scrollbar;
  float x, y, w, h;
  ArrayList<Gene> list;
  int selItem = -1; 
  boolean dragged;
  
  ScrollableList(float x, float y, float w, float h) {
    this.x = x; 
    this.y = y; 
    this.w = w;
    this.h = h;
  }
  
  void setList(ArrayList<Gene> lst) {
    list = lst;
    scrollbar = new ScrollBar(50 * list.size(), 0.1 * w, w, h);
    selItem = -1;
  }
  
  void display() {
    if (list == null || list.size() == 0) return;
    
    pushMatrix();
    translate(x, y);
    pushMatrix();
    translate(0, scrollbar.translateY);
    noStroke();
    for (int i = 0; i < list.size(); i++) {
      fill(210);
      float rx = 20;
      float ry = i * itemHeight + itemSpace;
      float rw = w - 40;
      float rh = itemHeight - itemSpace;
      if (selItem == i) {
        stroke(240, 118, 104);
      } else {
        noStroke();
      }
      rect(rx, ry, rw, rh);
      fill(120);
      Gene gene = list.get(i);
      String text = gene.name + " " + nf(gene.r.floatValue(), 1, 2);
      text(text, rx, ry, rw, rh);
    }
    popMatrix();
    scrollbar.display();
    popMatrix();
  }
  
  void press() {
    dragged = false;
    scrollbar.open();
  }

  void drag(int my, int pmy) {
    dragged = true;
    scrollbar.update(my - pmy);
  }
  
  int release(int my) {
    scrollbar.close();
    if (!dragged) {      
      float len = my - scrollbar.translateY;
      selItem = int(len / itemHeight);      
    } 
    if (selItem != -1) {
      return list.get(selItem).idx;      
    } else {
      return -1;
    }
  }
}

class ScrollBar {
  float totalHeight;
  float translateY;
  float opacity;
  float barWidth;
  float listWidth;
  float listHeight;

  ScrollBar(float th, float bw, float lw, float lh) {
    totalHeight = th;
    barWidth = bw;
    translateY = 0;
    opacity = 0;
    
    listWidth = lw;
    listHeight = lh;
  }

  void open() {
    opacity = 150;
  }

  void close() {
    opacity = 0;
  }

  void update(float dy) {       
    if (totalHeight + translateY + dy > listHeight) {
      translateY += dy;
      if (translateY > 0) translateY = 0;      
    }
  }

  void display() {
    if (0 < opacity) {
      float frac = (listHeight / totalHeight);
      float x = listWidth - 1.5 * barWidth;
      float y = PApplet.map(translateY / totalHeight, -1, 0, listHeight, 0);
      float w = barWidth;
      float h = frac * listHeight;
      pushStyle();
      noStroke();
      fill(150, opacity);
      rect(x, y, w, h, 0.2 * w);
      popStyle();
    }
  }
}
