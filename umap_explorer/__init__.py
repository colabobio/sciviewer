import pandas as pd
import numpy as np
import numpy.linalg as la
import scipy.stats as ss
import sys
import py5
from py5 import Sketch

usePShape = True

CLOSED = 0
SET_SPINE = 1
SET_WIDTH = 2
COMPLETED = 3

SEL_COLOR = 1
EXP_COLOR = 2
RST_COLOR = 3

itemHeight = 50
itemSpace = 10

def angle_between(v1, v2):
    # Angle calculation A
#     amt = np.dot(v1, v2) / (la.norm(v1) * la.norm(v2))
#     if amt <= -1:
#         return -py5.PI
#     elif amt >= 1:
#         return 0
#     return np.arccos(amt)
    # Angle calculation B
    cosang = np.dot(v1, v2)
    sinang = la.norm(np.cross(v1, v2))
    return np.arctan2(sinang, cosang)

class Gene():
    def __init__(self, n, i, r, p):
        self.name = n
        self.idx = i
        self.r = r
        self.rabs = abs(r)
        self.p = p

class ScrollableList:    
    def __init__(self, x, y, w, h):
        self.x = x 
        self.y = y 
        self.w = w
        self.h = h
        self.scrollbar = None
        self.genes = []
        self.selItem = -1
        self.dragged = False
  
    def set_list(self, genes):
        self.genes = genes
        self.scrollbar = ScrollBar(50 * len(genes), 0.1 * self.w, self.w, self.h)
        self.selItem = -1
  
    def display(self, py5obj):
        if not self.genes or len(self.genes) == 0: return
    
        py5obj.push_matrix()
        py5obj.translate(self.x, self.y)
        py5obj.push_matrix()
        py5obj.translate(0, self.scrollbar.translateY)
        py5obj.no_stroke()
        for i in range(0, len(self.genes)):
            py5obj.fill(210)
            rx = 20
            ry = i * itemHeight + itemSpace
            rw = self.w - 40
            rh = itemHeight - itemSpace
            if self.selItem == i:
                py5obj.stroke(240, 118, 104)
            else:
                py5obj.no_stroke()
            py5obj.rect(rx, ry, rw, rh)
            py5obj.fill(120)
            gene = self.genes[i]
            text = gene.name + " " + "{:1.2f}".format(gene.r)
            py5obj.text(text, rx, ry, rw, rh)
        py5obj.pop_matrix()
        self.scrollbar.display(py5obj)
        py5obj.pop_matrix()
  
    def press(self):
        self.dragged = False
        self.scrollbar.set_open()

    def drag(self, my, pmy):
        self.dragged = True
        self.scrollbar.update(pmy - my)

    def release(self, my):
        self.scrollbar.set_close()
        if not self.dragged:
            l = my - self.scrollbar.translateY
            self.selItem = int(l / itemHeight)

        if self.selItem != -1:
            return self.genes[self.selItem].idx  
        else:
            return -1

class ScrollBar:
    def __init__(self, th, bw, lw, lh):
        self.totalHeight = th
        self.barWidth = bw
        self.translateY = 0
        self.opacity = 0    
        self.listWidth = lw
        self.listHeight = lh

    def set_open(self):
        self.opacity = 150

    def set_close(self):
        self.opacity = 0

    def update(self, dy):
        if self.totalHeight + self.translateY + dy > self.listHeight:
            self.translateY += dy
            if self.translateY > 0: self.translateY = 0

    def display(self, py5obj):
        if 0 < self.opacity:
            frac = self.listHeight / self.totalHeight
            x = self.listWidth - 1.5 * self.barWidth
            y = py5obj.remap(self.translateY / self.totalHeight, -1, 0, self.listHeight, 0)
            w = self.barWidth
            h = frac * self.listHeight
            py5obj.push_style()
            py5obj.no_stroke()
            py5obj.fill(150, self.opacity)
            py5obj.rect(x, y, w, h, 0.2 * w)
            py5obj.pop_style()
            
class Button:
    def __init__(self, x, y, w, h, l):
        self.x = x
        self.y = y
        self.w = w
        self.h = h
        self.label = l
  
    def display(self, py5obj):
        py5obj.no_stroke()
        py5obj.fill(120)
        py5obj.rect(self.x, self.y, self.w, self.h, 15)
    
        py5obj.fill(255)
        py5obj.text(self.label, self.x, self.y, self.w, self.h)
  
    def contains(self, mx, my):
        return self.x <= mx and mx <= self.x + self.w and self.y <= my and my <= self.y + self.h

class Selector():
    def __init__(self):
        self.state = CLOSED
        self.spx0 = 0
        self.spy0 = 0
        self.spx1 = 0
        self.spy1 = 0
        self.wx = 0 
        self.wy = 0
        self.nx0 = 0 
        self.ny0 = 0        
        self.nx1 = 0
        self.ny1 = 0
        self.angle = 0
        self.x0 = 0
        self.y0 = 0
        self.w = 0
        self.h = 0
        self.tmat = np.array([[0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0]])

    def display(self, p5obj):
        if self.state == CLOSED: 
            return
        
        if self.state == SET_SPINE:
            p5obj.stroke(240, 118, 104)
            p5obj.line(self.spx0, self.spy0, self.spx1, self.spy1)
        elif self.state == SET_WIDTH:
            p5obj.stroke(240, 118, 104)
            p5obj.line(self.spx0, self.spy0, self.spx1, self.spy1)
            p5obj.line(self.spx1, self.spy1, self.wx, self.wy)
            self.display_box(p5obj)
        else:
            self.display_box(p5obj)
    
    
    def apply(self, p5obj, cell, sx0, sy0, sw, sh):
        sx = p5obj.remap(cell.umap1, 0, 1, sx0, sx0 + sw)
        sy = p5obj.remap(cell.umap2, 0, 1, sy0, sy0 + sh)
        
        # Transformation code A
#         s = np.array([sx, sy, 1])
#         t = np.matmul(self.tmat, s)
#         tx = t[0]
#         ty = t[1]

        # Transformation code B
        tx = self.mult_x(sx, sy)
        ty = self.mult_y(sx, sy)
        
        cell.selected = 0 <= tx and tx <= self.w and -self.h/2 <= ty and ty <= self.h/2

    def update_box(self, x, y):
        self.wx = x
        self.wy = y
        
        spdir = np.array([self.spx1 - self.spx0, self.spy1 - self.spy0])
        bxdir = np.array([self.wx - self.spx1, self.wy - self.spy1])

        a = angle_between(spdir, bxdir)
        d = np.sin(a) * la.norm(bxdir)
    
        self.angle = np.arctan2(spdir[1], spdir[0])
    
        self.h = 2 * d
        self.w = la.norm(spdir)
        self.x0 = self.spx0
        self.y0 = self.spy0
    
        # Transformation code A
#         s = np.sin(-self.angle)
#         c = np.cos(-self.angle)
#         tx = -self.x0 * c + self.y0 * s
#         ty = -self.x0 * s - self.y0 * c
#         self.tmat = np.array([[c, -s, tx],
#                               [s,  c, ty]])

        # Transformation code B
        self.reset()
        self.rotate(-self.angle)
        self.translate(-self.x0, -self.y0)

    def normalize(self, py5obj, sx0, sy0, sw, sh):
        self.nx0 = py5obj.remap(self.spx0, sx0, sx0 + sw, 0, 1)
        self.nx1 = py5obj.remap(self.spx1, sx0, sx0 + sw, 0, 1)
        self.ny0 = py5obj.remap(self.spy0, sy0, sy0 + sh, 0, 1)
        self.ny1 = py5obj.remap(self.spy1, sy0, sy0 + sh, 0, 1)
  
    def display_box(self, py5obj):
        py5obj.stroke(240, 118, 104)
        py5obj.no_fill()
        py5obj.push_matrix()
        py5obj.translate(self.x0, self.y0)
        py5obj.rotate(self.angle)
        py5obj.rect(0, -self.h/2, self.w, self.h)
        py5obj.pop_matrix()
        
    def press(self, x, y):
        if self.state == CLOSED or self.state == COMPLETED:
            self.state = SET_SPINE
        elif self.state == SET_SPINE:
            self.state = SET_WIDTH

        if self.state == SET_SPINE:
            self.spx0 = self.spx1 = self.wx = x
            self.spy0 = self.spy1 = self.wy = y
            self.angle = 0
            self.w = self.h = 0
  
    def drag(self, x, y):
        if self.state == SET_SPINE:
            self.spx1 = x
            self.spy1 = y      

    def move(self, x, y):
        if self.state == SET_WIDTH:
            self.update_box(x, y)

    def release(self, x, y):
        requestSelection = False
        if self.state == SET_SPINE:
            self.spx1 = x
            self.spy1 = y
            if self.spx1 != self.spx0 or self.spy1 != self.spy0:
                self.state = SET_WIDTH
                self.wx = self.spx1
                self.wy = self.spy1
            else:
                self.state = CLOSED
        elif self.state == SET_WIDTH:
            self.update_box(x, y)
            self.state = COMPLETED
            requestSelection = True
        return(requestSelection)
            
    def mult_x(self, x, y):
        return self.tmat[0][0] * x + self.tmat[0][1] * y + self.tmat[0][2]
            
    def mult_y(self, x, y):
        return self.tmat[1][0] * x + self.tmat[1][1] * y + self.tmat[1][2]
    
    def reset(self):
        self.tmat = np.array([[1.0, 0.0, 0.0],
                              [0.0, 1.0, 0.0]])        
            
    def rotate(self, angle):
        s = np.sin(angle)
        c = np.cos(angle)
        
        temp1 = self.tmat[0][0]
        temp2 = self.tmat[0][1]
        self.tmat[0][0] =  c * temp1 + s * temp2
        self.tmat[0][1] = -s * temp1 + c * temp2
        temp1 = self.tmat[1][0]
        temp2 = self.tmat[1][1]
        self.tmat[1][0] =  c * temp1 + s * temp2
        self.tmat[1][1] = -s * temp1 + c * temp2
        
    def translate(self, tx, ty):
        self.tmat[0][2] = tx*self.tmat[0][0] + ty*self.tmat[0][1] + self.tmat[0][2]
        self.tmat[1][2] = tx*self.tmat[1][0] + ty*self.tmat[1][1] + self.tmat[1][2]
        
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
  
    def init_expression(self, numGenes):
        self.expression = [0.0] * numGenes

    def set_expression(self, i, level):
        self.expression[i] = level

    def set_all_expressions(self, levels):
        self.expression = levels
        
    def project(self, sel):
        if self.selected:
            dirv = np.array([sel.nx1 - sel.nx0, sel.ny1 - sel.ny0])
            celv = np.array([self.umap1 - sel.nx0, self.umap2 - sel.ny0])
            a = angle_between(dirv, celv)        
            self.proj = np.cos(a) * la.norm(celv) / la.norm(dirv)
    
    def display(self, p5obj, x0, y0, w, h):
        x = p5obj.remap(self.umap1, 0, 1, x0, x0 + w)
        y = p5obj.remap(self.umap2, 0, 1, y0, y0 + h)
    
        p5obj.no_stroke()
        if p5obj.selGene != -1:
            f = p5obj.constrain(p5obj.remap(self.expression[p5obj.selGene],
                                            p5obj.minGeneExp, p5obj.maxGeneExp, 0, 1), 0, 1)
            p5obj.color_mode(p5obj.HSB, 360, 100, 100)
            p5obj.fill((1 - f) * 170 + f * 233, 74, 93, 80)
            p5obj.color_mode(p5obj.RGB, 255, 255, 255)
        else:
            if self.selected:
                p5obj.fill(240, 118, 104, 80)
            else:
                p5obj.fill(150, 80)
        p5obj.ellipse(x, y, 5, 5)

    def get_expr_color(self, p5obj, selGene):
        p5obj.color_mode(p5obj.HSB, 360, 100, 100)
        f = p5obj.constrain(p5obj.remap(self.expression[selGene], 
                                        p5obj.minGeneExp, p5obj.maxGeneExp, 0, 1), 0, 1)
        cl = p5obj.color((1 - f) * 170 + f * 233, 74, 93, 80)
        p5obj.color_mode(p5obj.RGB, 255, 255, 255)
        return cl

    def create_shape(self, p5obj, x0, y0, w, h):
        x = p5obj.remap(self.umap1, 0, 1, x0, x0 + w)
        y = p5obj.remap(self.umap2, 0, 1, y0, y0 + h)        
        sh = p5obj.create_shape(p5obj.ELLIPSE, x, y, 5, 5)
        sh.set_stroke(False)
        sh.set_fill(p5obj.color(150, 80))
        return sh
        
class UMAPexplorer(Sketch):
    def settings(self):
        self.size(1600, 800, self.P2D)

    def setup(self):
        self.text_align(self.CENTER, self.CENTER)
        self.init_UI()
        if usePShape:
            self.init_shape()
        self.text_font(self.create_font("Helvetica", 14))

    def draw(self):        
        self.background(255)
        self.show_UMAP_scatter()

        if self.requestSelection and 0 < len(self.indices):
            self.calculateGeneCorrelations()

        self.selector.display(self)
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
                print("Selected gene", self.selGene)
                self.calculateGeneMinMax()
                if usePShape:
                    self.color_shape(EXP_COLOR)
                
        elif self.exportBtn.contains(self.mouse_x, self.mouse_y):
            self.exportData()
            
    def init_UI(self):
        self.selector = Selector()
        self.scrollList = ScrollableList(self.width/2, 0, 200, self.height)  
        w = self.width - (self.width/2 + 200)
        self.exportBtn = Button(self.width/2 + 200 + w/2 - 75, self.height - 75, 100, 30, "EXPORT")   

    def init_shape(self):
        x0 = 25
        y0 = 25
        w = self.width/2 - 50
        h = self.height - 50

        self.cellShape = self.create_shape(self.GROUP)
        for cell in self.cells:
            sh = cell.create_shape(self, x0, y0, w, h)
            self.cellShape.add_child(sh)        

    def color_shape(self, mode):
        if mode == RST_COLOR:
            self.cellShape.set_fill(self.color(150, 80))
        elif mode == SEL_COLOR:
            for idx in range(0, len(self.cells)):
                cell = self.cells[idx]
                sh = self.cellShape.get_child(idx)
                if cell.selected:
                    cl = self.color(240, 118, 104, 80)                
                else:
                    cl = self.color(150, 80)            
                sh.set_fill(cl)
        elif mode == EXP_COLOR:
            for idx in range(0, len(self.cells)):
                cell = self.cells[idx]
                sh = self.cellShape.get_child(idx)            
                sh.set_fill(cell.get_expr_color(self, self.selGene))
            
    def show_UMAP_scatter(self):
        x0 = 25
        y0 = 25
        w = self.width/2 - 50
        h = self.height - 50
        
        
        if self.requestSelection:
            self.indices = []
            self.selector.normalize(self, x0, y0, w, h)        
        
        for idx in range(0, len(self.cells)):
            cell = self.cells[idx]
            if self.requestSelection:            
                self.selector.apply(self, cell, x0, y0, w, h)
                cell.project(self.selector)
                if cell.selected:
                    self.indices += [idx]
            if not usePShape:
                cell.display(self, x0, y0, w, h)
        
        if usePShape:
            self.shape(self.cellShape)        
        
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
        
        
    def calculateGeneCorrelations(self):
        print("Selected", len(self.indices), "cells")

        print("Calculating correlations...") 

        vproj = []
        rexpr = []
        for i in range(0, len(self.indices)):
            c = self.indices[i]
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

        self.scrollList.set_list(self.sortedGenes)

        self.requestSelection = False

        self.selGene = -1
        if usePShape:
            self.color_shape(RST_COLOR)
        print("Done")

    def calculateGeneMinMax(self):
        self.minGeneExp = sys.float_info.max
        self.maxGeneExp = sys.float_info.min
        for i in range(0, len(self.indices)):
            idx = self.indices[i]
            cell = self.cells[idx]
            exp = cell.expression[self.selGene]
            self.minGeneExp = min(self.minGeneExp, exp)
            self.maxGeneExp = max(self.maxGeneExp, exp)   
        print("Min/max expression level for gene", self.geneNames[self.selGene], self.minGeneExp, self.maxGeneExp)
        
    def showGeneScatter(self):
        x0 = self.width/2 + 200 + 50
        w = self.width - x0 - 100
        h = w
        y0 = (self.height - h) / 2

        for i in range(0, len(self.indices)):      
            idx = self.indices[i]
            cell = self.cells[idx]
            x = self.remap(cell.proj, 0, 1, x0 + 5, x0 + w - 5)
            y = self.remap(cell.expression[self.selGene], self.minGeneExp, self.maxGeneExp, y0 + w - 5, y0 + 5)
            self.no_stroke()
            self.fill(150, 80)
            self.ellipse(x, y, 10, 10)

        self.fill(100)
        self.text(self.geneNames[self.selGene], x0, 0, w, y0)

        self.stroke_weight(2)
        self.stroke(120)
        self.no_fill()
        self.rect(x0, y0, w, h)

        self.fill(130)
        self.text("{:1.2f}".format(self.maxGeneExp), x0 - 20, y0 + 5)
        self.text("{:1.2f}".format(self.minGeneExp), x0 - 20, y0 + h - 5)
        self.push_matrix()
        self.translate(x0 - 20, y0 + h/2)
        self.rotate(-self.HALF_PI)
        self.text("Expression", 0, 0)
        self.pop_matrix()

        self.text("0", x0 + 5, y0 + h + 15)
        self.text("1", x0 + w - 5, y0 + h + 15)
        self.text("Projection", x0 + 5, y0 + h + 10, w - 10, 20)
        
    def exportData(self):
        print("EXPORTING DATA...")
        rows = []
        for i in range(0, len(self.indices)):
            idx = self.indices[i]
            cell = self.cells[idx]
            row = [cell.code, cell.proj]
            rows += [row]
        self.selected_cells = pd.DataFrame.from_records(rows, columns=['index', 'proj'])

        rows = []
        for gene in self.sortedGenes:
            row = [gene.r, gene.p]
            rows += [row]
        self.significant_genes = pd.DataFrame.from_records(rows, columns=['R', 'P'])

        self.selected_gene_name = self.geneNames[self.selGene]

        rows = []
        for i in range(0, len(self.indices)):
            idx = self.indices[i]
            cell = self.cells[idx]
            row = [cell.code, cell.proj, cell.expression[self.selGene]]
            rows += [row]        
        self.selected_gene_cell_data = pd.DataFrame.from_records(rows, columns=['index', 'proj', 'exp'])        

        print("BYE")
        self.exit_sketch()
            
    def __init__(self, umap, expr):
        super().__init__()
        self.umap = umap
        self.expr = expr
        self.requestSelection = False
        
        self.cells = []
        self.geneNames = expr.columns.tolist()

        self.selected_cells = []
        self.significant_genes = []
        self.selected_gene_name = ''
        self.selected_gene_cell_data = ''
        self.selGene = -1
        self.sortedGenes = []
        
        self.minGeneExp = sys.float_info.max
        self.maxGeneExp = sys.float_info.min

        self.pearsonsThreshold = 0.1
        self.pvalueThreshold = 0.05
        
        min1 = umap["UMAP_1"].min()
        max1 = umap["UMAP_1"].max()
        min2 = umap["UMAP_2"].min()
        max2 = umap["UMAP_2"].max()

        cells = []
        for i in umap.index:
            cell = Cell(i, umap.at[i,'UMAP_1'], umap.at[i,'UMAP_2'])
            cell.normalize(self, min1, max1, min2, max2)
            cell.set_all_expressions(expr.loc[i].tolist())
            cells += [cell]
        self.cells = cells
        
        self.cellShape = None
    