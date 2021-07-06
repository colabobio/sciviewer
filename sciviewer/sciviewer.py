import pandas as pd
import numpy as np
import matplotlib as mp
import numpy.linalg as la
import scipy.stats as ss
import scipy.sparse as sp
import sys
import py5
from py5 import Sketch
import anndata
from .gui import ScrollableList, ScrollBar, Button, ToggleButton, Selector
from .utils import angle_between
#from gui import ScrollableList, ScrollBar, Button, ToggleButton, Selector
#from utils import angle_between

import time

SEL_COLOR = 1
EXP_COLOR = 2
RST_COLOR = 3

DEF_WIDTH = 1600
DEF_HEIGHT = 800

GENE_WIDTH = 200
MARGIN = 50
PADDING = 50

FONT_SIZE = 14

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
            a = angle_between(dirv, celv)        
            self.proj = np.cos(a) * la.norm(celv) / la.norm(dirv)

    def create_shape(self, p5obj, x0, y0, w, h):
        x = p5obj.remap(self.umap1, 0, 1, x0, x0 + w)
        y = p5obj.remap(self.umap2, 0, 1, y0, y0 + h)        
        sh = p5obj.create_shape(p5obj.ELLIPSE, x, y, 5, 5)
        sh.set_stroke(False)
        sh.set_fill(p5obj.color(150, 80))
        return sh        
       
class Py5Renderer(Sketch):
    def __init__(self, dataobj, width=DEF_WIDTH, height=DEF_HEIGHT):
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
        self.modeChange = False
        self.violingSelStats = None
        self.violingExlStats = None
        self.maxValue = 0
        self.minCoord = 0
        self.maxCoord = 0
        self.viewer_width = width
        self.viewer_height = height
        self.hscale = self.viewer_width / DEF_WIDTH
        self.vscale = self.viewer_height / DEF_HEIGHT
    
    def settings(self):
        self.size(self.viewer_width, self.viewer_height, self.P2D)

    def setup(self):
        surface = self.get_surface()
        surface.set_resizable(True)
        surface.set_title("SCIViewer")

        self.background(255)
        self.text_align(self.CENTER, self.CENTER)
        self.initUI()
        self.initUMAPshape()
        self.text_font(self.create_font("Helvetica", FONT_SIZE))

    def draw(self):
        self.background(255)

        self.viewer_width = self.width
        self.viewer_height = self.height        
        self.hscale = self.viewer_width / DEF_WIDTH
        self.vscale = self.viewer_height / DEF_HEIGHT
        self.scale(self.hscale, self.vscale)

        if self.selectedGene:        
            self.data.calculateGeneMinMax(self.indices, self.selGene)
            if self.modeBtn.state == 1:
                self.initScatterShape()
            else:
                self.initViolinShape()
            self.colorUMAPShape(EXP_COLOR)            
            self.selectedGene = False
        
        if self.modeChange and (len(self.indices) > 0) and (self.selGene != -1):
            if (self.modeBtn.state == 1):
                self.initScatterShape()
            else:
                self.initViolinShape()

        self.showUMAPScatter()

        if self.requestSelection and (len(self.indices) > 0):
            self.data.update_selected_cells(self.indices)
            if self.modeBtn.state == 1:
                # Correlation of expression with projection onto UMAP axis
                start = time.time()
                sortedGenes = self.data.calculateGeneCorrelations(self.indices)
                end = time.time()
                print(end - start,'seconds to calculate correlations. Sparsity: ', self.data.sparse)
            else:
                # Calculate differential expression
                start = time.time()
                sortedGenes = self.data.calculateDiffExpr(self.indices)
                end = time.time()
                print(end - start,'seconds to calculate differential expression. Sparsity: ', self.data.sparse)             

            self.scrollList.setList(sortedGenes,
                                    maxposgenes=self.data.maxdisplaygenes_pos,
                                    maxneggenes=self.data.maxdisplaygenes_neg)
            if len(sortedGenes) == 0:
                print('No genes pass the statistical significance / effect size thresholds.')
                print('Please make another selection or re-run with lower thresholds')
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
                self.showGeneViolinPlot()

        self.modeBtn.display(self)
        self.exportBtn.display(self)

    def mouse_pressed(self):
        if self.mouse_x < self.width/2:
            self.selector.press(self.mouse_x, self.mouse_y, self.width, self.height)
        elif self.scrollList.contains(self.mouse_x, self.mouse_y, self.hscale, self.vscale):
            self.scrollList.press()

    def mouse_dragged(self):
        if self.mouse_x < self.width/2:
            self.selector.drag(self.mouse_x, self.mouse_y, self.width, self.height)
        elif self.mouse_x < self.width/2 + 1.5 * GENE_WIDTH:
            self.scrollList.drag(self.mouse_x, self.pmouse_y, self.hscale, self.vscale)

    def mouse_moved(self):
        if self.mouse_x < self.width/2:
            self.selector.move(self.mouse_x, self.mouse_y, self.width, self.height)

    def mouse_released(self):
        if self.mouse_x < self.width/2:
            self.requestSelection = self.selector.release(self.mouse_x, self.mouse_y, self.width, self.height)
        elif self.exportBtn.contains(self.mouse_x, self.mouse_y, self.hscale, self.vscale):
            self.data.export_data(self.indices, self.selGene)
        else: 
            pstate = self.modeBtn.state
            if self.modeBtn.contains(self.mouse_x, self.mouse_y, self.hscale, self.vscale):
                self.modeChange = pstate != self.modeBtn.state

        sel = self.scrollList.release(self.mouse_x, self.mouse_y, self.hscale, self.vscale)
        if sel != -1 and sel != self.selGene:
            self.selGene = sel
            self.selectedGene = True
            print("Selected gene", self.data.geneNames[self.selGene])

    def initUI(self):
        self.selector = Selector(self.width, self.height)
        self.scrollList = ScrollableList(DEF_WIDTH/2, 0, GENE_WIDTH, DEF_HEIGHT)

        w = 200
        x0 = DEF_WIDTH/2 + 100 + DEF_WIDTH/4 - w/2
        self.exportBtn = Button(x0, DEF_HEIGHT - 55, w, 30, "EXPORT & CLOSE")  

        w = 250
        x0 = DEF_WIDTH/2 + 100 + DEF_WIDTH/4 - w/2        
        self.modeBtn = ToggleButton(x0, 25, w, 30, "DIRECTIONAL", "DIFFERENTIAL")  
        
    def initUMAPshape(self):
        x0 = MARGIN/2 + PADDING
        y0 = MARGIN/2 + PADDING
        w = DEF_WIDTH/2 - MARGIN - 2 * PADDING
        h = DEF_HEIGHT - MARGIN - 2 * PADDING

        self.umapShape = self.create_shape(self.GROUP)
        for cell in self.data.cells:
            sh = cell.create_shape(self, x0, y0, w, h)
            self.umapShape.add_child(sh)

    def initScatterShape(self):
        x0 = DEF_WIDTH/2 + GENE_WIDTH + MARGIN
        w = DEF_WIDTH - x0 - MARGIN
        h = w
        y0 = (DEF_HEIGHT - h) / 2        
        self.scatterShape = self.create_shape(self.GROUP)
        expr = self.data.expr[self.indices, self.selGene]
        if self.data.sparse:
            expr = expr.toarray().reshape(-1)
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
        if self.data.sparse:
            selected_values = selected_values.toarray().reshape(-1)
            excluded_values = excluded_values.toarray().reshape(-1)
        
        # selected_values_nonzero = [x for x in selected_values if x != 0]
        # excluded_values_nonzero = [x for x in excluded_values if x != 0]
        
        self.violingSelStats = mp.cbook.violin_stats(selected_values, _kde_method)
        self.violingExlStats = mp.cbook.violin_stats(excluded_values, _kde_method)        

        self.maxValue = 0
        self.minCoord = +100000
        self.maxCoord = -100000        
        for all_stats in [self.violingSelStats, self.violingExlStats]:
            for stats in all_stats:
                self.minCoord = min(self.minCoord, stats['coords'].min())
                self.maxCoord = max(self.maxCoord, stats['coords'].max())                
                self.maxValue = max(self.maxValue, stats['vals'].max())

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
            if self.data.sparse:
                color_grad = self.data.expr[:, self.selGene].toarray().reshape(-1)
            else:
                color_grad = self.data.expr[:, self.selGene].copy()
            min_expr = self.data.minGeneExp
            max_expr = self.data.maxGeneExp
            color_grad -= min_expr
            color_grad /= (max_expr-min_expr)
            self.color_mode(self.HSB, 360, 100, 100)
            for idx in range(self.data.num_cells):
                sh = self.umapShape.get_child(idx)
                cl = self.color(color_grad[idx] * 126 + (1 - color_grad[idx]) * 233, 85, 95, 80)
                sh.set_fill(cl)
            self.color_mode(self.RGB, 255, 255, 255)

    def showUMAPScatter(self):
        x0 = MARGIN/2 + PADDING
        y0 = MARGIN/2 + PADDING
        w = DEF_WIDTH/2 - MARGIN - 2 * PADDING
        h = DEF_HEIGHT - MARGIN - 2 * PADDING
      
        if self.requestSelection:
            self.indices = []
            self.selector.normalize(self, self.hscale * x0, self.vscale * y0, self.hscale * w, self.vscale * h)
            start = time.time()
            for idx in range(0, len(self.data.cells)):
                cell = self.data.cells[idx]
                self.selector.apply(self, cell, self.hscale * x0, self.vscale * y0, self.hscale * w, self.vscale * h)
                cell.project(self.selector)
                if cell.selected:
                    self.indices += [idx]

            setindices = set(self.indices)
            num_cells = self.data.num_cells
            self.excluded_indices = [i for i in range(num_cells) if i not in setindices]

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
            x00 = MARGIN/2
            y00 = MARGIN/2
            self.no_stroke()
            for i in range(0, 20):
                f = self.remap(i, 0, 19, 0, 1)
                self.color_mode(self.HSB, 360, 100, 100)
                self.fill(f * 126 + (1 - f) * 233, 85, 95, 80)
                self.color_mode(self.RGB, 255, 255, 255)
                x = self.remap(f, 0, 1, x00 + 20, x00 + 120)
                self.rect(x, y00 + 20, 100.0/19, 30)
            self.fill(130)
            self.text("Max exp.", x00 + 160, y00 + 35)

        self.fill(130)
        self.text("UMAP1", x0, y0 + h + 2.5/2, w, DEF_HEIGHT - y0 - h)
        self.push_matrix()
        self.translate((x0-2.5)/2, y0 + h/2)
        self.rotate(-self.HALF_PI)
        self.text("UMAP2", 0, 0)
        self.pop_matrix()
        
    def showGeneScatter(self):
        x0 = DEF_WIDTH/2 + GENE_WIDTH + MARGIN
        w = DEF_WIDTH - x0 - MARGIN
        h = w
        y0 = (DEF_HEIGHT - h) / 2        

        self.fill(100)
        self.text("Selected gene: " + self.data.geneNames[self.selGene], x0, 55, w, y0 - 55)

        self.shape(self.scatterShape)

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

    def showGeneViolinPlot(self):
        x0 = DEF_WIDTH/2 + GENE_WIDTH + MARGIN
        w = DEF_WIDTH - x0 - MARGIN
        h = w
        y0 = (DEF_HEIGHT - h) / 2

        self.fill(100)
        self.text("Selected gene: " + self.data.geneNames[self.selGene], x0, 55, w, y0 - 55)

        self.no_stroke()
        self.fill(self.color(240, 118, 104, 80))
        self.violin(self.violingSelStats, x0, y0, w/2, h)
        self.fill(self.color(150, 80))
        self.violin(self.violingExlStats, x0 + w/2, y0, w/2, h)

        self.stroke_weight(2)
        self.stroke(120)
        self.no_fill()
        self.rect(x0, y0, w, h)

        self.fill(130)
        self.text("{:1.2f}".format(self.data.maxGeneExp), x0 - 20, y0 + 5)
        self.text("{:1.2f}".format(self.data.minGeneExp), x0 - 20, y0 + h - 5)
        self.text("SELECTED CELLS", x0, y0 + h + 15, w/2, 2 * FONT_SIZE)
        self.text("EXCLUDED CELLS", x0 + w/2, y0 + h + 15, w/2, 2 * FONT_SIZE)

    # Adapted from https://matplotlib.org/stable/_modules/matplotlib/axes/_axes.html#Axes.violin for drawing the actual violins
    def violin(self, vpstats, x0, y0, w, h):
        N = len(vpstats)
        positions = range(1, N + 1)

        widths = 0.5
        widths = [widths] * N
        scalef = 0.4 * w / self.maxValue
        for stats, pos, width in zip(vpstats, positions, widths):
            self.no_stroke()
            self.begin_shape(self.QUAD_STRIP)
            for nx, ny in zip(stats['vals'], stats['coords']):
                y = self.remap(ny, self.minCoord, self.maxCoord, y0 + h, y0)
                cx = x0 + w/2
                self.vertex(cx - scalef * nx, y)
                self.vertex(cx + scalef * nx, y)
            self.end_shape()            

    def setClip(self):
        x0 = 25
        y0 = 25
        w = DEF_WIDTH/2 - MARGIN
        h = DEF_HEIGHT - MARGIN
        self.clip(self.hscale * (x0 - 2.5), self.vscale * (y0 - 2.5), self.hscale * (w + 5), self.vscale * (h + 5))

    def delClip(self):
        self.no_clip()        

class SCIViewer():
    """Invetarctive visualizer for 2D embeddings of single-cell gene expression
    
    When the visualizer is opened by the explore_data method, displays the 2D embedding with the option to select cells and a direction on the 2D plane. When the 'directional' option is toggled, it will compute the correlation of all genes with the direction for the selected cells and will display a list of the most correlated genes. Selecting a gene will show a scatter plot of each cell's position projected onto the direction, against its expression of the selected gene. Alternatively, if the 'differential' option is toggled, it will compute the T-test differential expression of the selected cells against all the non-selected cells and display those that are most differentially expressed. Selecting a gene will show a violin plot of the expression of that gene in the selected cells compared against all of the other cells. Press the Export and Close button to close the visualizer when you are done.
    
    Public Methods
    ------------------
    explore_data() - starts the visualizer
    
    Required arguments
    ------------------
    expr - (anndata.AnnData, pandas.DataFrame, numpy.ndarray, scipy.sparse.csc_matrix) - Either an NxG matrix of cells by genes or a Scanpy AnnData object containing an NxG matrix in the expr.X or expr.raw.X fields. Currently only supports expression data from Scanpy in numpy.ndarray (non-sparse) or sparse/csc_matrix format so please convert the attribute accordingly if it is using a different sprase format. The method was tested on log2(TP10K) data but other similar normalizations should work as well. If a Pandas DataFrame is provided, the indeces are used as cell names and the columns are used as gene names. If an AnnData is provided cell names are obtained from the index of the obs attribute and gene names are obtained from the index of the var attribute. Otherwise, these can be provided as the gene_names and cell_names optional arguments. If a csc_matrix is provided or the provided AnnData data matrix is sparse, all of the computations will exploit the sparsity of the data which can provide a substantial (orders of magnitude) speedup for large datasets.
    
    Optional arguments:
    ------------------
    umap (pandas.DataFrame, numpy.ndarray) - Nx2 matrix of cells by 2D embedding coordinates (can be any embedding including tsne, umap...). REQUIRED if expr is not a scanpy AnnData object, or if it doesn't include an embedding in the obsm attribute.
    gene_names (list of strings) -- gene names. (default: inferred from DataFrame columns or ['0', '1' ...'G'])
    cell_names (list of strings) -- cell names. (default: inferred from DataFrame indeces or ['0', '1' ...'N'])
    embedding_name (str) -- Only used if an AnnData provided for expr. Key value for obsm specifying embedding to use. (default: 'X_umap') 
    pearsonsThreshold (float)-- Pearson correlation threshold (absolute value) for genes to display in the scroll list for the directional analysis (default: 0.1)
    tThreshold (float) -- T-statistic association threshold (absolute value) for genes to display in the scroll list for the differential analysis (default: 2.0)
    pvalueThreshold (float) -- P-value  threshold (absolute value) for genes to display in the scroll list for both differential and directional analysis (default: 0.05)
    maxdisplaygenes_pos (int) -- maximum # of positively associated genes to show in the scroll list (by either directional or differential analysis) (default: 100)
    maxdisplaygenes_neg (int) -- maximum # of negatively associated genes to show in the scroll list (by either directional or differential analysis)  (default: 100)
    use_raw (boolean) -- if expr is AnnData, indicates to use .raw.X rather than .X for expression matrix (default: False)
    
    Attributes
    ----------
    selected_cells : pandas.DataFrame
        information about the currently selected cells.
        
    results_proj_correlation : pandas.DataFrame
        results from the most recent directional analysis (Pearson correlation and P-value for all genes). Updated whenever a new selection is made in the 'directional' mode
        
    results_diffexpr : pandas.DataFrame
        results from the most recent differential expression analysis (T-statistic and P-value for all genes). Updated whenever a new selection is made in the 'differential' mode
        
    significant_genes : pandas.DataFrame
        results from the most recent differential or directional analysis, whatever was run most recently. Only includes genes that pass the significance / effect size thresholds. Updated when the 'Export and close' button is pressed.
    
    selected_gene_cell_data : pandas.DataFrame
        includes the gene expression and projection coordinates for the selected cells at the time the 'Export and close' button is pressed.
    """
    
    def __init__(self, expr, umap=None, gene_names=None, cell_names=None,
                embedding_name='X_umap', pearsonsThreshold=0.1, tThreshold=2.0,
                pvalueThreshold=0.05, maxdisplaygenes_pos=100, maxdisplaygenes_neg=100,
                width=1600, height=800, use_raw=False):       
      
        if type(expr) is anndata._core.anndata.AnnData:
            if umap is None:
                umap = expr.obsm[embedding_name]
            
            if (gene_names is None):
                if not use_raw:
                    gene_names = expr.var.index.tolist()
                else:
                    gene_names = expr.raw.var.index.tolist()
                
            if (cell_names is None):
                cell_names = expr.obs.index.tolist()
                
            self.anndata_obj = expr
            if not use_raw: expr = expr.X
            else: expr = expr.raw.X
        
        if type(umap) is pd.core.frame.DataFrame: self.umap = umap.values
        elif type(umap) is np.ndarray: self.umap = umap
        else: sys.exit('umap must be a Pandas DataFrame or Numpy ndarray')

        if type(expr) is pd.core.frame.DataFrame:
            self.expr = expr.values
            self.geneNames = expr.columns.tolist()
            self.cellNames = expr.index.tolist()
        elif type(expr) in [np.ndarray, sp.csc.csc_matrix, sp.csr.csr_matrix]:
            self.expr = expr
                        
            if gene_names is None: self.geneNames = [str(i) for i in np.arange(expr.shape[1])]
            else: self.geneNames = gene_names

            if cell_names is None: self.cellNames = np.arange(expr.shape[0])
            else: self.cellNames = cell_names
            
        else:
            sys.exit('Expression argument - expr - must be pd.DataFrame, np.ndarray, sp.csc_matrix, or sp.csc.csr_matrix')
        end = time.time()

        if self.umap.shape[0] != self.expr.shape[0]:
            sys.exit('# of cells not equal between 2D embedding and expression matrix inputs')
        
        if type(expr) is sp.csc.csc_matrix: self.sparse = True
        else: self.sparse = False

        self.cells = []

        self.selected_cells = []
        self.significant_genes = None
        self.sortedGenes = []
        self.selected_gene_name = ''
        self.selected_gene_cell_data = None
        
        self.gene_sum = None
        self.gene_sqsum = None

        self.pearsonsThreshold = pearsonsThreshold
        self.tThreshold = tThreshold
        self.pvalueThreshold = pvalueThreshold
        self.maxdisplaygenes_pos = maxdisplaygenes_pos
        self.maxdisplaygenes_neg = maxdisplaygenes_neg
        
        self.results_diffexpr = None
        self.results_proj_correlation = None
        
        min1 = self.umap[:,0].min()
        max1 = self.umap[:,0].max()
        min2 = self.umap[:,1].min()
        max2 = self.umap[:,1].max()

        self.renderer = Py5Renderer(self, width=width, height=height)
        
        start = time.time()
        cells = []
        for i in range(self.umap.shape[0]):
            cell = Cell(self.cellNames[i], self.umap[i,0], self.umap[i,1])
            cell.normalize(self.renderer, min1, max1, min2, max2)
            cells.append(cell)
        self.cells = cells
        self.num_cells = len(cells)
        
    def update_selected_cells(self, indices):
        celldata = [[i, self.cells[i].code, self.cells[i].proj] for i in indices]
        self.selected_cells = pd.DataFrame(celldata, columns=['index', 'cell_name', 'projection'])
        
    def calculateDiffExpr(self, indices):
        print("Selected", len(indices), "cells")

        print("Calculating differential expression...")
                
        # Calculate gene sums & sq-sums if not already calculated
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
        
        self.results_diffexpr = pd.DataFrame([T, P], columns=self.geneNames,
                                                           index=['T', 'P']).T
        
        self.sortedGenes = []
        for g in range (0, len(self.geneNames)):
            if self.pearsonsThreshold <= abs(T[g]) and P[g] <= self.pvalueThreshold:
                gene = Gene(self.geneNames[g], g, T[g], P[g])
                self.sortedGenes += [gene]

        self.sortedGenes.sort(key=lambda x: x.r, reverse=True)
        return(self.sortedGenes)

    def calculateGeneCorrelations(self, indices):
        print("Selected", len(indices), "cells")

        print("Calculating correlations...") 

        self.sortedGenes = []
        vproj = []
        for i in indices:
            vproj.append(self.cells[i].proj)
            
        dexpr = self.expr[indices, :]
        ## only run correlation on genes with non-zero counts in the selected cell subset
        hascounts = np.array(dexpr.sum(axis=0)).reshape(-1)>0
        dexpr = dexpr[:,hascounts]
        vproj = np.array(vproj)
        n = vproj.size
        geneNames = np.array(self.geneNames)[hascounts]
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

        # obtain Pandas DataFrame with P-values / Corelations with nans for genes not tested
        self.results_proj_correlation = pd.DataFrame([rs, ps], columns=geneNames,
                                                           index=['R', 'P']).T
        nullvalues = pd.DataFrame(index=np.array(self.geneNames)[~hascounts], columns=['R', 'P'])
        self.results_proj_correlation = pd.concat([self.results_proj_correlation, nullvalues], axis=0)
        self.results_proj_correlation = self.results_proj_correlation.loc[self.geneNames, :]
        
        for (i,g) in enumerate(self.geneNames):
            r = self.results_proj_correlation.at[g, 'R']
            p = self.results_proj_correlation.at[g, 'P']
            if (self.pearsonsThreshold <= abs(r)) and (p <= self.pvalueThreshold):
                gene = Gene(g, i, r, p)
                self.sortedGenes.append(gene)

        self.sortedGenes.sort(key=lambda x: x.r, reverse=True)

        return(self.sortedGenes)

    def calculateGeneMinMax(self, indices, selGene):
        self.minGeneExp = self.expr[:,selGene].min()
        self.maxGeneExp = self.expr[:,selGene].max()
        #self.minGeneExp = self.expr[indices,selGene].min()
        #self.maxGeneExp = self.expr[indices,selGene].max()  
        print("Min/max expression level for gene", self.geneNames[selGene], self.minGeneExp, self.maxGeneExp)

    def export_data(self, indices, selGene):
        print("EXPORTING DATA...")        
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
