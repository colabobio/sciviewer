import pandas as pd
import numpy as np
import matplotlib as mp
import numpy.linalg as la
import scipy.stats as ss
import scipy.sparse as sp
import sys
import py5
from py5 import Sketch
from gui import ScrollableList, ScrollBar, Button, ToggleButton, Selector
from utils import angleBetween

import time

SEL_COLOR = 1
EXP_COLOR = 2
RST_COLOR = 3

GENE_WIDTH = 200
MARGIN = 50
PADDING = 50


class Gene():
    def __init__(self, n, i, r, p):
        self.name = n
        self.idx = i
        self.r = r
        self.rabs = abs(r)
        self.p = p

        
class Cell:
    def __init__(self, c, u1, u2):
        self.code = c
        self.umap1 = u1
        self.umap2 = u2
        self.proj = 0
        self.selected = False
        self.expression = []

    def normalize(self, p5obj, min1, max1, min2, max2):
        self.umap1 = p5obj.remap(self.umap1, min1, max1, 0, 1)
        self.umap2 = p5obj.remap(self.umap2, min2, max2, 1, 0)
        
    def project(self, sel):
        if self.selected:
            dirv = np.array([sel.nx1 - sel.nx0, sel.ny1 - sel.ny0])
            celv = np.array([self.umap1 - sel.nx0, self.umap2 - sel.ny0])
            a = angleBetween(dirv, celv)        
            self.proj = np.cos(a) * la.norm(celv) / la.norm(dirv)

    def createShape(self, p5obj, x0, y0, w, h):
        x = p5obj.remap(self.umap1, 0, 1, x0, x0 + w)
        y = p5obj.remap(self.umap2, 0, 1, y0, y0 + h)        
        sh = p5obj.create_shape(p5obj.ELLIPSE, x, y, 5, 5)
        sh.set_stroke(False)
        sh.set_fill(p5obj.color(150, 80))
        return sh        
       
class py5renderer(Sketch):
    def __init__(self, dataobj):
        super().__init__()
        self.data = dataobj
        self.indices = []
        self.excluded_indices = []
        self.selGene = -1
        self.selectedGene = False
        self.requestSelection = False
        self.umapShape = None
        self.scatterShape = None
        self.exportBtn = None
        self.modeBtn = None
    
    def settings(self):
        self.size(1600, 800, self.P2D)

    def setup(self):
        self.text_align(self.CENTER, self.CENTER)
        self.initUI()
        self.initUMAPShape()
        self.text_font(self.create_font("Helvetica", 14))

    def draw(self):        
        self.background(255)

        if self.selectedGene:        
            self.data.calculateGeneMinMax(self.indices, self.selGene)
            if self.modeBtn.state == 1:
                self.initScatterShape()
            else:
                self.initViolinShape()
            self.colorUMAPShape(EXP_COLOR)            
            self.selectedGene = False
                
        self.showUMAPScatter()

        if self.requestSelection and 1 < len(self.indices):
            if self.modeBtn.state == 1:
                # Correlation of expression with projection onto UMAP axis
                start = time.time()
                sortedGenes = self.data.calculateGeneCorrelations(self.indices)
                end = time.time()
                print(end - start,'seconds to calculate correlations. Sparsity: ', self.data.sparse)
            else:
                # Obtain set of unselected cells for violin plots
                setindices = set(self.indices)
                num_cells = self.data.num_cells
                self.excluded_indices = [i for i in range(num_cells) if i not in setindices]
                
                # Calculate differential expression
                start = time.time()
                sortedGenes = self.data.calculateDiffExpr(self.indices)
                end = time.time()
                print(end - start,'seconds to calculate differential expression. Sparsity: ', self.data.sparse)             

            self.scrollList.setList(sortedGenes)
            self.requestSelection = False
            self.selGene = -1
            self.colorUMAPShape(RST_COLOR)               

        self.setClip()
        self.selector.display(self)
        self.delClip()

        self.scrollList.display(self)

        if self.selGene != -1:
            if self.modeBtn.state == 1:
                self.showGeneScatter()
            else:
                # Need to make violin plot visual
                self.showGeneViolin()

        self.modeBtn.display(self)
        self.exportBtn.display(self)

    def mouse_pressed(self):
        if self.mouse_x < self.width/2:
            self.selector.press(self.mouse_x, self.mouse_y)
        elif self.scrollList.contains(self.mouse_x, self.mouse_y):
            self.scrollList.press()

    def mouse_dragged(self):
        if self.mouse_x < self.width/2:
            self.selector.drag(self.mouse_x, self.mouse_y)
        elif self.mouse_x < self.width/2 + 1.5 * GENE_WIDTH:
            self.scrollList.drag(self.mouse_x, self.pmouse_y)

    def mouse_moved(self):
        if self.mouse_x < self.width/2:
            self.selector.move(self.mouse_x, self.mouse_y)

    def mouse_released(self):
        if self.mouse_x < self.width/2:
            self.requestSelection = self.selector.release(self.mouse_x, self.mouse_y)
        elif self.exportBtn.contains(self.mouse_x, self.mouse_y):
            self.data.exportData(self.indices, self.selGene)
        elif self.modeBtn.contains(self.mouse_x, self.mouse_y):
            print("Selected mode", self.modeBtn.state)

        sel = self.scrollList.release(self.mouse_x, self.mouse_y)
        if sel != -1 and sel != self.selGene:
            self.selGene = sel
            self.selectedGene = True
            print("Selected gene", self.data.geneNames[self.selGene])

    def initUI(self):
        self.selector = Selector()
        self.scrollList = ScrollableList(self.width/2, 0, GENE_WIDTH, self.height)

        w = 100
        x0 = self.width/2 + 100 + self.width/4 - w/2
        self.exportBtn = Button(x0, self.height - 75, w, 30, "EXPORT")  

        w = 250
        x0 = self.width/2 + 100 + self.width/4 - w/2        
        self.modeBtn = ToggleButton(x0, 25, w, 30, "DIRECTIONAL", "DIFFERENTIAL")  
        
    def initUMAPShape(self):
        x0 = MARGIN/2 + PADDING
        y0 = MARGIN/2 + PADDING
        w = self.width/2 - MARGIN - 2 * PADDING
        h = self.height - MARGIN - 2 * PADDING

        self.umapShape = self.create_shape(self.GROUP)
        for cell in self.data.cells:
            sh = cell.createShape(self, x0, y0, w, h)
            self.umapShape.add_child(sh)

    def initScatterShape(self):
        x0 = self.width/2 + GENE_WIDTH + MARGIN
        w = self.width - x0 - MARGIN
        h = w
        y0 = (self.height - h) / 2        
        self.scatterShape = self.create_shape(self.GROUP)
        expr = self.data.expr[self.indices, self.selGene]
        for i in range(0, len(self.indices)):
            idx = self.indices[i]
            cell = self.data.cells[idx]
            x = self.remap(cell.proj, 0, 1, x0 + 5, x0 + w - 5)
            y = self.remap(expr[i], self.data.minGeneExp, self.data.maxGeneExp, y0 + w - 5, y0 + 5)
            sh = self.create_shape(self.ELLIPSE, x, y, 10, 10)
            sh.set_stroke(False)
            sh.set_fill(self.color(150, 80))
            self.scatterShape.add_child(sh)
           
    def initViolinShape(self, bw_method=None):
        def _kde_method(X, coords):
            if np.all(X[0] == X):
                return (X[0] == coords).astype(float)
            kde = mp.mlab.GaussianKDE(X, bw_method)
            return kde.evaluate(coords)
        
        selected_values = self.data.expr[self.indices, self.selGene]
        excluded_values = self.data.expr[self.excluded_indices, self.selGene]
        
        selected_values_nonzero = [x for x in selected_values if x != 0]
        excluded_values_nonzero = [x for x in excluded_values if x != 0]
        
        self.selected_stats = mp.cbook.violin_stats(selected_values, _kde_method)
        self.excluded_stats = mp.cbook.violin_stats(excluded_values, _kde_method)
        ## See https://matplotlib.org/stable/_modules/matplotlib/axes/_axes.html#Axes.violin for drawing the actual violins
        
        '''
        x0 = self.width/2 + 200 + 50
        w = self.width - x0 - 100
        h = w
        y0 = (self.height - h) / 2        
        self.scatterShape = self.create_shape(self.GROUP)        
        for i in range(0, len(self.indices)):
            idx = self.indices[i]
            cell = self.data.cells[idx]
            x = self.remap(cell.proj, 0, 1, x0 + 5, x0 + w - 5)
            y = self.remap(cell.expression[self.selGene], self.data.minGeneExp, self.data.maxGeneExp, y0 + w - 5, y0 + 5)
            sh = self.create_shape(self.ELLIPSE, x, y, 10, 10)
            sh.set_stroke(False)
            sh.set_fill(self.color(150, 80))
            self.scatterShape.add_child(sh)
        '''
            
            
            
            
            
            

    def colorUMAPShape(self, mode):
        if mode == RST_COLOR:
            self.umapShape.set_fill(self.color(150, 80))
        elif mode == SEL_COLOR:
            for idx in range(0, len(self.data.cells)):
                cell = self.data.cells[idx]
                sh = self.umapShape.get_child(idx)
                if cell.selected:
                    cl = self.color(240, 118, 104, 80)                
                else:
                    cl = self.color(150, 80)            
                sh.set_fill(cl)
        elif mode == EXP_COLOR:
            color_grad = self.data.expr[:, self.selGene]
            minGeneExp = self.data.minGeneExp
            maxGeneExp = self.data.maxGeneExp           
            color_grad -= minGeneExp
            color_grad /= (maxGeneExp-minGeneExp)
            self.color_mode(self.HSB, 360, 100, 100)
            print(color_grad.min(), color_grad.max())
            for idx in range(self.data.num_cells):
                sh = self.umapShape.get_child(idx)
                cl = self.color((1 - color_grad[idx]) * 170 + color_grad[idx] * 233, 74, 93, 80)
                sh.set_fill(cl)            
            self.color_mode(self.RGB, 255, 255, 255)

    def showUMAPScatter(self):
        x0 = MARGIN/2 + PADDING
        y0 = MARGIN/2 + PADDING
        w = self.width/2 - MARGIN - 2 * PADDING
        h = self.height - MARGIN - 2 * PADDING
      
        if self.requestSelection:
            self.indices = []
            self.selector.normalize(self, x0, y0, w, h)
            start = time.time()
            for idx in range(0, len(self.data.cells)):
                cell = self.data.cells[idx]
                self.selector.apply(self, cell, x0, y0, w, h)
                cell.project(self.selector)
                if cell.selected:
                    self.indices += [idx]
            end = time.time()
            print(end - start,'seconds to select and project cells')
                    
            if len(self.indices) == 0:
                print('No cells selected, please make another selection')
                self.requestSelection = False
        
        self.shape(self.umapShape)
        self.stroke_weight(2)
        self.stroke(120)
        self.no_fill()
        self.rect(x0 - PADDING, y0 - PADDING, w + 2 * PADDING, h + 2 * PADDING)

        if self.selGene != -1:
            self.no_stroke()
            for i in range(0, 20):
                f = self.remap(i, 0, 19, 0, 1)
                self.color_mode(self.HSB, 360, 100, 100)
                self.fill((1 - f) * 170 + f * 233, 74, 93, 80)
                self.color_mode(self.RGB, 255, 255, 255)
                x = self.remap(f, 0, 1, x0 + 20, x0 + 120)
                self.rect(x, y0 + 20, 100.0/19, 30)
            self.fill(130)
            self.text("Max exp.", x0 + 160, y0 + 35)

        self.fill(130)
        self.text("UMAP1", x0, y0 + h + 2.5/2, w, self.height - y0 - h)
        self.push_matrix()
        self.translate((x0-2.5)/2, y0 + h/2)
        self.rotate(-self.HALF_PI)
        self.text("UMAP2", 0, 0)
        self.pop_matrix()
        
    def showGeneScatter(self):
        x0 = self.width/2 + GENE_WIDTH + MARGIN
        w = self.width - x0 - MARGIN
        h = w
        y0 = (self.height - h) / 2

        self.shape(self.scatterShape)

        self.fill(100)
        self.text("Selected gene: " + self.data.geneNames[self.selGene], x0, 55, w, y0 - 55)

        self.stroke_weight(2)
        self.stroke(120)
        self.no_fill()
        self.rect(x0, y0, w, h)

        self.fill(130)
        self.text("{:1.2f}".format(self.data.maxGeneExp), x0 - 20, y0 + 5)
        self.text("{:1.2f}".format(self.data.minGeneExp), x0 - 20, y0 + h - 5)
        self.push_matrix()
        self.translate(x0 - 20, y0 + h/2)
        self.rotate(-self.HALF_PI)
        self.text("Expression", 0, 0)
        self.pop_matrix()

        self.text("0", x0 + 5, y0 + h + 15)
        self.text("1", x0 + w - 5, y0 + h + 15)
        self.text("Projection", x0 + 5, y0 + h + 10, w - 10, 20)

    def setClip(self):
        x0 = 25
        y0 = 25
        w = self.width/2 - MARGIN
        h = self.height - MARGIN       
        self.clip(x0 - 2.5, y0 - 2.5, w + 5, h + 5)

    def delClip(self):
        self.no_clip()        

    def showGeneViolin(self):
        pass        

class UMAPexplorer():
    
    def __init__(self, umap, expr, gene_names=None, cell_names=None):    
        if type(umap) is pd.core.frame.DataFrame: self.umap = umap.values
        elif type(umap) is np.ndarray: self.umap = umap
        else: sys.exit('umap must be a Pandas DataFrame or Numpy ndarray')

        if type(expr) is pd.core.frame.DataFrame:
            self.expr = expr.values
            self.geneNames = expr.columns.tolist()
            self.cellNames = expr.index.tolist()
        elif type(expr) in [np.ndarray, sp.csc.csc_matrix]:
            self.expr = expr
                        
            if gene_names is None: self.geneNames = [str(i) for i in np.arange(expr.shape[1])]
            else: self.geneNames = gene_names

            if cell_names is None: self.cellNames = np.arange(expr.shape[0])
            else: self.cellNames = cell_names
            
        else:
            sys.exit('Expression argument - expr - must be pd.DataFrame, np.ndarray, or sp.csc_matrix')
        end = time.time()

        if self.umap.shape[0] != self.expr.shape[0]:
            sys.exit('# of cells not equal between 2D embedding and expression matrix inputs')
        
        if type(expr) is sp.csc.csc_matrix: self.sparse = True
        else: self.sparse = False

        self.cells = []
        self.sortedGenes = []

        self.selected_cells = []
        self.significant_genes = []
        self.sortedGenes = []
        self.selected_gene_name = ''
        self.selected_gene_cell_data = ''
        
        self.gene_sum = None
        self.gene_sqsum = None

        self.pearsonsThreshold = 0.1
        self.tThreshold = 2.0
        self.pvalueThreshold = 0.05
        
        min1 = self.umap[:,0].min()
        max1 = self.umap[:,0].max()
        min2 = self.umap[:,1].min()
        max2 = self.umap[:,1].max()

        self.renderer = py5renderer(self)
        
        start = time.time()
        cells = []
        for i in range(self.umap.shape[0]):
            cell = Cell(self.cellNames[i], self.umap[i,0], self.umap[i,1])
            cell.normalize(self.renderer, min1, max1, min2, max2)
            cells.append(cell)
        self.cells = cells
        self.num_cells = len(cells)
        
    def calculateDiffExpr(self, indices):
        print("Selected", len(indices), "cells")

        print("Calculating differential expression...")
                
        # Calculate gene sums / sq-sums if not already calculated
        if self.gene_sum is None:
            start = time.time()
            self.gene_sum = self.expr.sum(axis=0)
            end = time.time()
            print(end - start,'seconds to calculate genesums. Sparsity: ', self.sparse)
            start = time.time()
            if self.sparse:
                self.gene_sum = np.array(self.gene_sum).reshape(-1)
                self.expr2 = self.expr.copy()
                self.expr2.data **= 2
                self.gene_sqsum = np.array(self.expr2.sum(axis=0)).reshape(-1)
            else:
                self.gene_sqsum = (self.expr**2).sum(axis=0)
            end = time.time()
            print(end - start,'seconds to calculate squared genesums. Sparsity: ', self.sparse)
        
        selected_N = len(indices)
        remainder_N = self.num_cells - selected_N
        
        selected_means = self.expr[indices,:].sum(axis=0)
        if self.sparse:
            selected_means = np.array(selected_means).reshape(-1)
        
        remainder_means = (self.gene_sum - selected_means) / remainder_N
        selected_means = selected_means / selected_N
        
        if self.sparse:
            selected_stds = np.array(self.expr2[indices,:].sum(axis=0)).reshape(-1)
        else:
            selected_stds = (self.expr[indices,:]**2).sum(axis=0)
            
        remainder_stds = np.sqrt((self.gene_sqsum - selected_stds - (remainder_N*remainder_means**2)) / (remainder_N -1))
        selected_stds = np.sqrt((selected_stds - selected_N*selected_means**2) / (selected_N -1))
        
        (T, P) = ss.ttest_ind_from_stats(selected_means, selected_stds, selected_N, remainder_means, remainder_stds, remainder_N, equal_var=False, alternative='two-sided')

        self.sortedGenes = []
        for g in range (0, len(self.geneNames)):
            if self.pearsonsThreshold <= abs(T[g]) and P[g] <= self.pvalueThreshold:
                gene = Gene(self.geneNames[g], g, T[g], P[g])
                self.sortedGenes += [gene]

        self.sortedGenes.sort(key=lambda x: x.r, reverse=False)
        return(self.sortedGenes)

    def calculateGeneCorrelations(self, indices):
        print("Selected", len(indices), "cells")

        print("Calculating correlations...") 

        self.sortedGenes = []
        vproj = []
        for i in indices:
            vproj.append(self.cells[i].proj)
            
        dexpr = self.expr[indices, :]
        vproj = np.array(vproj)
        n = vproj.size
        if self.sparse:
            # based on https://www.javaer101.com/en/article/18344934.html
            yy = vproj - vproj.mean()
            xm = dexpr.mean(axis=0).A.ravel()
            ys = yy / np.sqrt(np.dot(yy, yy))
            xs = np.sqrt(np.add.reduceat(dexpr.data**2, dexpr.indptr[:-1]) - n*xm*xm)
            rs = np.add.reduceat(dexpr.data * ys[dexpr.indices], dexpr.indptr[:-1]) / xs
        else:
            # based on https://github.com/ikizhvatov/efficient-columnwise-correlation/blob/master/columnwise_corrcoef_perf.py
            DO = dexpr - (np.sum(dexpr, 0) / np.double(n))
            DP = vproj - (np.sum(vproj) / np.double(n))        
            rs = np.dot(DP, DO) / np.sqrt(np.sum(DO ** 2, 0) * np.sum(DP ** 2))
        
        # calculate P-value
        T = -1*np.abs(rs * np.sqrt(n-2))/np.sqrt(1 - (rs**2))
        ps = ss.t.cdf(T, df=n-2)*2
        
        for (i,g) in enumerate(self.geneNames):
            if (self.pearsonsThreshold <= abs(rs[i])) and (ps[i] <= self.pvalueThreshold):
                gene = Gene(g, i, rs[i], ps[i])
                self.sortedGenes.append(gene)
        
        self.sortedGenes.sort(key=lambda x: x.r, reverse=False)
        return(self.sortedGenes)

    def calculateGeneMinMax(self, indices, selGene):
        self.minGeneExp = self.expr[:,selGene].min()
        self.maxGeneExp = self.expr[:,selGene].max()
        #self.minGeneExp = self.expr[indices,selGene].min()
        #self.maxGeneExp = self.expr[indices,selGene].max()  
        print("Min/max expression level for gene", self.geneNames[selGene], self.minGeneExp, self.maxGeneExp)

    def exportData(self, indices, selGene):
        print("EXPORTING DATA...")
        rows = []
        for i in range(0, len(indices)):
            idx = indices[i]
            cell = self.cells[idx]
            row = [cell.code, cell.proj]
            rows += [row]
        self.selected_cells = pd.DataFrame.from_records(rows, columns=['index', 'proj'])
        
        rows = []
        gene_names = []
        for gene in self.sortedGenes:
            row = [gene.r, gene.p]
            rows += [row]
            gene_names.append(gene.name)
        self.significant_genes = pd.DataFrame(rows, columns=['R', 'P'], index=gene_names)

        self.selected_gene_name = self.geneNames[selGene]

        rows = []
        expr = self.expr[indices, selGene]
        for i in range(0, len(indices)):
            idx = indices[i]
            cell = self.cells[idx]
            row = [cell.code, cell.proj, expr[i]]
            rows += [row]        
        self.selected_gene_cell_data = pd.DataFrame.from_records(rows, columns=['index', 'proj', 'exp'])        

        print("BYE")
        self.renderer.exit_sketch()
        

    def explore_data(self):
        self.renderer.run_sketch()
