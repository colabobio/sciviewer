import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import java.util.Collections;
import java.util.Comparator;

ArrayList<Cell> cells;
ArrayList<String> geneNames;
ArrayList<Gene> sortedGenes;
ArrayList<Integer> indices;
Selector selector;
int selGene = -1;
float minGeneExp, maxGeneExp;

ScrollableList scrollList;
boolean requestSelection;
Button exportBtn;

void setup() {
  size(1600, 800, P2D);
  textAlign(CENTER, CENTER);
  textFont(createFont("Helvetica", 14));
  frame.setTitle("EMBED VIEW");
  
  loadData();
  initUI();
}

void draw() {
  background(255);

  showUMAPScatter();
  
  if (requestSelection && 0 < indices.size()) {
    calculateGeneCorrelations();
  }
  
  selector.display();  
  scrollList.display();
    
  if (selGene != -1) {
    showGeneScatter();
    exportBtn.display();
  }
    
  //if (frameCount % 1000 == 0) println(frameRate);
}

void mousePressed() {
  if (mouseX < width/2) {
    selector.press(mouseX, mouseY);
  } else if (mouseX < width/2 + 200) {
    scrollList.press();
  }
}

void mouseDragged() {
  if (mouseX < width/2) {
    selector.drag(mouseX, mouseY);
  } else if (mouseX < width/2 + 200) {
    scrollList.drag(mouseY, pmouseY);
  }
}

void mouseMoved() {
  if (mouseX < width/2) {
    selector.move(mouseX, mouseY);
  }
}

void mouseReleased() {
  if (mouseX < width/2) {
    selector.release(mouseX, mouseY);
  } else if (mouseX < width/2 + 200) {
    int sel = scrollList.release(mouseY);
    if (sel != -1 && sel != selGene) {      
      selGene = sel;
      println("Selected gene", selGene);
      calculateGeneMinMax();
    }    
  } else if (exportBtn.contains(mouseX, mouseY)) {
    exportData();
  }
}
