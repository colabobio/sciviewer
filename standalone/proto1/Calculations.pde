void calculateGeneCorrelations() {
  println("Selected", indices.size(), "cells");
  
  double[][] data = new double[indices.size()][1 + geneNames.size()];
  for (int i = 0; i < indices.size(); i++) {
    int idx = indices.get(i);
    Cell cell = cells.get(idx);
    data[i][0] = cell.proj;
    for (int g = 0; g < geneNames.size(); g++) {
      data[i][1 + g] = cell.expression[g]; 
    }
  }
  
  println("Calculating Pearsons Correlation...");
  PearsonsCorrelation pearsons = new PearsonsCorrelation(data);
  RealMatrix corr = pearsons.getCorrelationMatrix();
  RealMatrix pval = pearsons.getCorrelationPValues();
  sortedGenes = new ArrayList<Gene>();
  for (int g = 1; g < corr.getColumnDimension(); g++) {
    //println(geneNames.get(g - 1), corr.getEntry(0, g), pval.getEntry(0, g));
    double r = corr.getEntry(0, g);
    double p = pval.getEntry(0, g); 
    if (pearsonsThreshold <= Math.abs(r) && p <= pvalueThreshold) {
      Gene gene = new Gene(geneNames.get(g - 1), g - 1, corr.getEntry(0, g), pval.getEntry(0, g));
      sortedGenes.add(gene);       
    }
  }
  
  Collections.sort(sortedGenes, new Comparator<Gene>() {
    @Override
    public int compare(Gene g1, Gene g2) {
      return g2.rabs.compareTo(g1.rabs);
    }
  });
    
  //for (Gene gene: sortedGenes) {
  //  println(gene.name, gene.r, gene.p);
  //}
  println("Sorted", sortedGenes.size(),"genes");
  scrollList.setList(sortedGenes);
  selGene = -1;
    
  println("Done.");
  requestSelection = false;
  //println("**********************************************");  
}

void calculateGeneMinMax() {
  minGeneExp = Float.MAX_VALUE;
  maxGeneExp = Float.MIN_VALUE;
  for (int i = 0; i < indices.size(); i++) {
    int idx = indices.get(i);
    Cell cell = cells.get(idx);
    float exp = cell.expression[selGene];
    minGeneExp = min(minGeneExp, exp);
    maxGeneExp = max(maxGeneExp, exp);    
  }
  println("Min/max expression level for gene", geneNames.get(selGene), minGeneExp, maxGeneExp);
}
