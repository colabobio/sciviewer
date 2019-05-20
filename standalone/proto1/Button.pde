class Button {
  float x, y, w, h;
  String label;
  
  Button(float x, float y, float w, float h, String label) {
    this.x = x;
    this.y = y;
    this.w = w;
    this.h = h;
    this.label = label;
  }
  
  void display() {
    noStroke();
    fill(120);
    rect(x, y, w, h, 15);
    
    fill(255);
    text(label, x, y, w, h);
  }
  
  boolean contains(int mx, int my) {
    return x <= mx && mx <= x + w && y <= my && my <= y + h;
  }
}
