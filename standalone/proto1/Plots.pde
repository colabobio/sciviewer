void showUMAPScatter() {
  float x0 = 25;
  float y0 = 25;
  float w = width/2 - 50;
  float h = height - 50;

  if (requestSelection) {
    indices.clear(); 
    selector.normalize(x0, y0, w, h);
  }
  
  for (int idx = 0; idx < cells.size(); idx++) {
    Cell cell = cells.get(idx);
    if (requestSelection) { 
      selector.apply(cell, x0, y0, w, h);
      cell.project(selector);
      if (cell.selected) {
        indices.add(idx);
      }
    }
    cell.display(x0, y0, w, h);
  }
  
  strokeWeight(2);
  stroke(120);
  noFill();
  rect(x0 - 2.5, y0 - 2.5, w + 5, h + 5);
  
  if (selGene != -1) {
    noStroke();
    for (int i = 0; i < 20; i++) {
      float f = map(i, 0, 19, 0, 1);
      colorMode(HSB, 360, 100, 100);
      fill((1 - f) * 170 + f * 233, 74, 93, 80);
      colorMode(RGB, 255, 255, 255);
      float x = map(f, 0, 1, x0 + 20, x0 + 120);
      rect(x, y0 + 20, 100.0/19, 30);
    }
    fill(130);
    text("Max exp.", x0 + 160, y0 + 35);
  }
  
  fill(130);
  text("UMAP1", x0, y0 + h + 2.5/2, w, height - y0 - h);
  pushMatrix();
  translate((x0-2.5)/2, y0 + h/2);
  rotate(-HALF_PI);
  text("UMAP2", 0, 0);
  popMatrix();  
}

void showGeneScatter() {
  float x0 = width/2 + 200 + 50;
  float w = width - x0 - 100;
  float h = w;
  float y0 = (height - h) / 2;
  
  for (int i = 0; i < indices.size(); i++) {
    int idx = indices.get(i);
    Cell cell = cells.get(idx);
    float x = map(cell.proj, 0, 1, x0 + 5, x0 + w - 5);
    float y = map(cell.expression[selGene], minGeneExp, maxGeneExp, y0 + w - 5, y0 + 5);    
    noStroke();
    fill(150, 80);
    ellipse(x, y, 10, 10);
  }

  fill(100);
  text(geneNames.get(selGene), x0, 0, w, y0);
  
  strokeWeight(2);
  stroke(120);
  noFill();
  rect(x0, y0, w, h);
  
  fill(130);
  text(nf(maxGeneExp, 1, 2), x0 - 20, y0 + 5);
  text(nf(minGeneExp, 1, 2), x0 - 20, y0 + h - 5);
  pushMatrix();
  translate(x0 - 20, y0 + h/2);
  rotate(-HALF_PI);
  text("Expression", 0, 0);
  popMatrix();
  
  text("0", x0 + 5, y0 + h + 15);
  text("1", x0 + w - 5, y0 + h + 15);
  text("Projection", x0 + 5, y0 + h + 10, w - 10, 20);
}
