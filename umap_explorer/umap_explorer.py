import pandas as pd
import numpy as np
import numpy.linalg as la
import scipy.stats as ss
import sys
import py5
from py5 import Sketch
from .gui import ScrollableList, ScrollBar, Button, Selector, angleBetween

SEL_COLOR = 1
EXP_COLOR = 2
RST_COLOR = 3

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
  
    def initExpression(self, numGenes):
        self.expression = [0.0] * numGenes

    def setExpression(self, i, level):
        self.expression[i] = level

    def setAllExpressions(self, levels):
        self.expression = levels
        
    def project(self, sel):
        if self.selected:
            dirv = np.array([sel.nx1 - sel.nx0, sel.ny1 - sel.ny0])
            celv = np.array([self.umap1 - sel.nx0, self.umap2 - sel.ny0])
            a = angleBetween(dirv, celv)        
            self.proj = np.cos(a) * la.norm(celv) / la.norm(dirv)

    def getExprColor(self, p5obj, minGeneExp, maxGeneExp, selGene):
        p5obj.color_mode(p5obj.HSB, 360, 100, 100)
        f = p5obj.constrain(p5obj.remap(self.expression[selGene], 
                                        minGeneExp, maxGeneExp, 0, 1), 0, 1)
        cl = p5obj.color((1 - f) * 170 + f * 233, 74, 93, 80)
        p5obj.color_mode(p5obj.RGB, 255, 255, 255)
        return cl

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
        self.selGene = -1
        self.selectedGene = False
        self.requestSelection = False
        self.umapShape = None
        self.scatterShape = None
    
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
            self.initScatterShape()
            self.colorUMAPShape(EXP_COLOR)
            self.selectedGene = False

        self.showUMAPScatter()

        if self.requestSelection and 0 < len(self.indices):
            self.data.calculateGeneCorrelations(self.indices)

        self.setClip()
        self.selector.display(self)
        self.delClip()

        self.scrollList.display(self)

        if self.selGene != -1:
            self.showGeneScatter()
        self.exportBtn.display(self)

    def mouse_pressed(self):
        if self.mouse_x < self.width/2:
            self.selector.press(self.mouse_x, self.mouse_y)
        elif self.mouse_x < self.width/2 + 200:
            self.scrollList.press()

    def mouse_dragged(self):
        if self.mouse_x < self.width/2:
            self.selector.drag(self.mouse_x, self.mouse_y)
        elif self.mouse_x < self.width/2 + 200:
            self.scrollList.drag(self.mouse_y, self.pmouse_y)

    def mouse_moved(self):
        if self.mouse_x < self.width/2:
            self.selector.move(self.mouse_x, self.mouse_y)

    def mouse_released(self):
        if self.mouse_x < self.width/2:
            self.requestSelection = self.selector.release(self.mouse_x, self.mouse_y)
        elif self.mouse_x < self.width/2 + 200:
            sel = self.scrollList.release(self.mouse_y)
            if sel != -1 and sel != self.selGene:
                self.selGene = sel
                self.selectedGene = True
                print("Selected gene", self.data.geneNames[self.selGene])

        elif self.exportBtn.contains(self.mouse_x, self.mouse_y):
            self.data.exportData(self.indices, self.selGene)        
        
    def initUI(self):
        self.selector = Selector()
        self.scrollList = ScrollableList(self.width/2, 0, 200, self.height)  
        w = self.width - (self.width/2 + 200)
        self.exportBtn = Button(self.width/2 + 200 + w/2 - 75, self.height - 75, 100, 30, "EXPORT")  
        
    def initUMAPShape(self):
        x0 = 25
        y0 = 25
        w = self.width/2 - 50
        h = self.height - 50

        self.umapShape = self.create_shape(self.GROUP)
        for cell in self.data.cells:
            sh = cell.createShape(self, x0, y0, w, h)
            self.umapShape.add_child(sh)

    def initScatterShape(self):
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
            for idx in range(0, len(self.data.cells)):
                cell = self.data.cells[idx]
                sh = self.umapShape.get_child(idx)            
                sh.set_fill(cell.getExprColor(self, self.data.minGeneExp, self.data.maxGeneExp, self.selGene))

    def showUMAPScatter(self):
        x0 = 25
        y0 = 25
        w = self.width/2 - 50
        h = self.height - 50
       
        if self.requestSelection:
            self.indices = []
            self.selector.normalize(self, x0, y0, w, h)        
        
        for idx in range(0, len(self.data.cells)):
            cell = self.data.cells[idx]
            if self.requestSelection:            
                self.selector.apply(self, cell, x0, y0, w, h)
                cell.project(self.selector)
                if cell.selected:
                    self.indices += [idx]
        
        self.shape(self.umapShape)

        self.stroke_weight(2)
        self.stroke(120)
        self.no_fill()
        self.rect(x0 - 2.5, y0 - 2.5, w + 5, h + 5)

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
        x0 = self.width/2 + 200 + 50
        w = self.width - x0 - 100
        h = w
        y0 = (self.height - h) / 2

        self.shape(self.scatterShape)

        self.fill(100)
        self.text(self.data.geneNames[self.selGene], x0, 0, w, y0)

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
        w = self.width/2 - 50
        h = self.height - 50       
        self.clip(x0 - 2.5, y0 - 2.5, w + 5, h + 5)

    def delClip(self):
        self.no_clip()        
     
class UMAPexplorer():
    def __init__(self, umap, expr):
        self.umap = umap
        self.expr = expr

        self.cells = []
        self.geneNames = expr.columns.tolist()
        self.sortedGenes = []

        self.selected_cells = []
        self.significant_genes = []
        self.selected_gene_name = ''
        self.selected_gene_cell_data = ''

        self.pearsonsThreshold = 0.1
        self.pvalueThreshold = 0.05
        
        min1 = umap["UMAP_1"].min()
        max1 = umap["UMAP_1"].max()
        min2 = umap["UMAP_2"].min()
        max2 = umap["UMAP_2"].max()

        self.renderer = py5renderer(self)
        
        cells = []
        for i in umap.index:
            cell = Cell(i, umap.at[i,'UMAP_1'], umap.at[i,'UMAP_2'])
            cell.normalize(self.renderer, min1, max1, min2, max2)
            cell.setAllExpressions(expr.loc[i].tolist())
            cells += [cell]
        self.cells = cells

    def calculateGeneCorrelations(self, indices):
        print("Selected", len(indices), "cells")

        print("Calculating correlations...") 

        self.sortedGenes = []

        vproj = []
        rexpr = []
        for i in range(0, len(indices)):
            c = indices[i]
            cell = self.cells[c]
            vproj += [cell.proj]
            rexpr += [cell.expression]
        dexpr = pd.DataFrame.from_records(rexpr)
        for g in range (0, len(self.geneNames)):
            r, p = ss.pearsonr(vproj, dexpr[g])
            if self.pearsonsThreshold <= abs(r) and p <= self.pvalueThreshold:
                gene = Gene(self.geneNames[g], g, r, p)
                self.sortedGenes += [gene]

        self.sortedGenes.sort(key=lambda x: x.r, reverse=False)

        self.renderer.scrollList.setList(self.sortedGenes)

        self.renderer.requestSelection = False

        self.renderer.selGene = -1
        self.renderer.colorUMAPShape(RST_COLOR)
        print("Done")

    def calculateGeneMinMax(self, indices, selGene):
        self.minGeneExp = sys.float_info.max
        self.maxGeneExp = sys.float_info.min
        for i in range(0, len(indices)):
            idx = indices[i]
            cell = self.cells[idx]
            exp = cell.expression[selGene]
            self.minGeneExp = min(self.minGeneExp, exp)
            self.maxGeneExp = max(self.maxGeneExp, exp)   
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
        self.significant_genes = pd.DataFrame.from_records(rows, columns=['R', 'P'], index=gene_names)

        self.selected_gene_name = self.geneNames[selGene]

        rows = []
        for i in range(0, len(indices)):
            idx = indices[i]
            cell = self.cells[idx]
            row = [cell.code, cell.proj, cell.expression[selGene]]
            rows += [row]        
        self.selected_gene_cell_data = pd.DataFrame.from_records(rows, columns=['index', 'proj', 'exp'])        

        print("BYE")
        self.renderer.exit_sketch()
        

    def explore_data(self):
        self.renderer.run_sketch()
