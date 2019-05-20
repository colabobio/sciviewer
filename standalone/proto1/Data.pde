void loadData() {
  println("LOADING UMAP DATA...");
  Table table = loadTable("pbmc3k_umap.tsv", "header");
  //Table table = loadTable("test.tsv", "header");
  cells = new ArrayList<Cell>();
  float min1 = Float.MAX_VALUE; 
  float max1 = Float.MIN_VALUE;
  float min2 = Float.MAX_VALUE; 
  float max2 = Float.MIN_VALUE;  
  for (TableRow row : table.rows()) {
    String code = row.getString("index");
    float umap1 = row.getFloat("UMAP_1");
    float umap2 = row.getFloat("UMAP_2");
    if (Float.isNaN(umap1) || Float.isNaN(umap2)) continue;
    min1 = min(min1, umap1);
    max1 = max(max1, umap1);
    min2 = min(min2, umap2);
    max2 = max(max2, umap2);     
    cells.add(new Cell(code, umap1, umap2));
  }
  println("NORMALIZING...");
  for (Cell cell : cells) {
    cell.normalize(min1, max1, min2, max2);
  }
  println("DONE: LOADED", table.getRowCount(), "CELLS");
  
  println("LOADING GENE EXPRESSION DATA...");
  table = loadTable("pbmc3k_expression_filtered_normalized.tsv", "header");
  int numCol = table.getColumnCount();
  geneNames = new ArrayList<String>();
  for (int c = 1; c < numCol; c++) {
    geneNames.add(table.getColumnTitle(c));
  }
  for (int r = 0; r < table.getRowCount(); r++) {
    TableRow row = table.getRow(r);    
    String code = row.getString("index");
    Cell cell = cells.get(r);
    if (!cell.code.equals(code)) {
      println("ERROR: data mistmatch at row", r);
      exit();
    }
    cell.initExpression(numCol - 1);
    for (int c = 1; c < numCol; c++) {
      float level = row.getFloat(c);
      cell.setExpression(c - 1, level);
    }    
  }
  println("DONE: LOADED", numCol - 1, "GENES FOR", table.getRowCount(), "CELLS");  
}

void exportData() {
  println("SAVING DATA...");
  Table table = new Table();
  table.addColumn("index");
  table.addColumn("proj");
  for (int i = 0; i < indices.size(); i++) {
    int idx = indices.get(i);
    Cell cell = cells.get(idx);
    TableRow newRow = table.addRow();
    newRow.setString("index", cell.code);
    newRow.setFloat("proj", cell.proj);
  }
  saveTable(table, "data/cells.csv");
  
  table = new Table();
  table.addColumn("name");
  table.addColumn("R");
  table.addColumn("P");
  for (Gene gene: sortedGenes) {
    TableRow newRow = table.addRow();
    newRow.setString("name", gene.name);
    newRow.setFloat("R", gene.r.floatValue());
    newRow.setFloat("P", gene.p.floatValue());
  }  
  saveTable(table, "data/genes.csv");
   
  table = new Table();
  table.addColumn("index");
  table.addColumn("proj");
  table.addColumn("exp");
  for (int i = 0; i < indices.size(); i++) {
    int idx = indices.get(i);
    Cell cell = cells.get(idx);
    TableRow newRow = table.addRow();
    newRow.setString("index", cell.code);
    newRow.setFloat("proj", cell.proj);
    newRow.setFloat("exp", cell.expression[selGene]);    
  }
  saveTable(table, "data/" + geneNames.get(selGene) + ".csv");
  
  println("DONE.");
}
