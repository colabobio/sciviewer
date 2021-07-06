import numpy as np
import numpy.linalg as la
from .utils import angle_between
#from utils import angle_between

CLOSED = 0
SET_SPINE = 1
SET_WIDTH = 2
COMPLETED = 3

ITEM_HEIGHT = 50
ITEM_SPACE = 10

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
        self.pressed = False
        self.visible = False
  
    def setList(self, genes, maxposgenes=100, maxneggenes=100):
        ## Set a ceiling on the number of positive/negatively associated
        ## genes that can be shown in the scrollbar
        rs = np.array([g.r for g in genes])
        numpos = (rs>0).sum()
        numneg = (rs<0).sum()
        if numpos > maxposgenes:
            displaygenes = genes[:maxposgenes]
        else:
            displaygenes = genes[:numpos]
    
        if numneg > maxneggenes:
            displaygenes += genes[(-1*maxneggenes):]
        else:
            displaygenes += genes[(-1*numneg):]
        
        self.genes = displaygenes
        toth = ITEM_HEIGHT * len(self.genes)
        self.visible = self.h < toth
         
        self.scrollbar = ScrollBar(0.9 * self.w, 0.1 * self.w, 0.1 * self.h, len(self.genes))
        self.selItem = -1
  
    def display(self, py5obj):
        if not self.genes or len(self.genes) == 0 or not self.visible: return
    
        py5obj.push_matrix()
        py5obj.translate(self.x, self.y)

        py5obj.push_matrix()
        py5obj.translate(0, self.scrollbar.translateY)
        py5obj.no_stroke()
        for i in range(0, len(self.genes)):
            py5obj.fill(210)
            rx = 20
            ry = i * ITEM_HEIGHT + ITEM_SPACE
            rw = 0.85 * self.w - 20
            rh = ITEM_HEIGHT - ITEM_SPACE
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

        py5obj.fill(220)
        py5obj.no_stroke()
        py5obj.rect(0.9 * self.w, 0, 0.1 * self.w, self.h)
        self.scrollbar.display(py5obj)

        py5obj.pop_matrix()
  
    def press(self):
        self.pressed = True

    def drag(self, mx, my, sx, sy):
        if self.pressed:
            self.dragged = True
            if sy * self.y <= my and my <= sy * (self.y + self.h):
                self.scrollbar.update(my, self.h, sy)

    def contains(self, mx, my, sx, sy):
         return sx * (self.x + 0.9 * self.w) <= mx and mx <= sx * (self.x + self.w)

    def release(self, mx, my, sx, sy):
        insideItemArea = sx * (self.x + 20) <= mx and mx <= sx * (self.x + 0.85 * self.w) and sy * self.y <= my and my <= sy * (self.y + self.h)
        insideDragArea = sx * (self.x + 0.9 * self.w) <= mx and mx <= sx * (self.x + self.w) and sy * self.y <= my and my <= sy * (self.y + self.h)
        if not self.dragged and insideItemArea:
            l = my/sy - self.scrollbar.translateY
            newItem = int(l / ITEM_HEIGHT)
            if self.selItem == newItem:
                # self.selItem = -1 # deselect
                pass
            else:
                self.selItem = newItem
        elif self.pressed and insideDragArea:
            self.scrollbar.update(my, self.h, sy)

        self.pressed = False
        self.dragged = False

        if self.selItem != -1:
            return self.genes[self.selItem].idx  
        else:
            return -1

class ScrollBar:    
    def __init__(self, x0, bw, bh, n):
        self.posX = x0
        self.posY = 0
        self.barWidth = bw
        self.barHeight = bh
        self.numItems = n
        self.translateY = 0

    def update(self, y, th, sy):        
        ymax = th - self.barHeight
        y0 = y / sy
        if y0 <= ymax:
            self.posY = y0
            self.translateY = -(y0 / ymax) * (ITEM_SPACE + ITEM_HEIGHT * self.numItems - th)

    def display(self, py5obj):        
        x = self.posX
        y = self.posY
        w = self.barWidth
        h = self.barHeight
        py5obj.push_style()
        py5obj.no_stroke()
        py5obj.fill(150)
        py5obj.rect(x, y, w, h, w)
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
  
    def contains(self, mx, my, sx, sy):
        return sx * self.x <= mx and mx <= sx * (self.x + self.w) and sy * self.y <= my and my <= sy * (self.y + self.h)


class ToggleButton:
    def __init__(self, x, y, w, h, l1, l2):
        self.x = x
        self.y = y
        self.w = w
        self.h = h
        self.label1 = l1
        self.label2 = l2
        self.state = 1
  
    def display(self, py5obj):
        py5obj.no_stroke()

        if self.state == 1:
            py5obj.fill(120)
        else:
            py5obj.fill(180)
        py5obj.rect(self.x, self.y, self.w/2, self.h, 15, 0, 0, 15)
    
        py5obj.fill(255)
        py5obj.text(self.label1, self.x, self.y, self.w/2, self.h)
  
        if self.state == 1:
            py5obj.fill(180)
        else:
            py5obj.fill(120)
        py5obj.rect(self.x + self.w/2, self.y, self.w/2, self.h, 0, 15, 15, 0)
    
        py5obj.fill(255)
        py5obj.text(self.label2, self.x + self.w/2, self.y, self.w/2, self.h)

    def contains(self, mx, my, sx, sy):
        inside = sx * self.x <= mx and mx <= sx * (self.x + self.w) and sy * self.y <= my and my <= sy * (self.y + self.h)
        if inside:
            if sx * self.x <= mx and mx <= sx * (self.x + self.w/2):
                self.state = 1
            else:
                self.state = 2
        return inside

class Selector():
    def __init__(self, w, h):
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
        self.screen_width = w
        self.screen_height = h
        self.tmat = np.array([[0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0]])

    def display(self, p5obj):
        if self.state == CLOSED: 
            return
        
        p5obj.push_matrix()
        p5obj.scale(1/p5obj.hscale, 1/p5obj.vscale)

        p5obj.scale(p5obj.width/self.screen_width, p5obj.height/self.screen_height)

        if self.state == SET_SPINE:
            p5obj.stroke(240, 118, 104)
            p5obj.line(self.spx0, self.spy0, self.spx1, self.spy1)
        elif self.state == SET_WIDTH:
            p5obj.stroke(240, 118, 104)
            p5obj.line(self.spx0, self.spy0, self.spx1, self.spy1)
            p5obj.line(self.spx1, self.spy1, self.wx, self.wy)
            self.displayBox(p5obj)
        else:
            self.displayBox(p5obj)
    
        p5obj.pop_matrix()
    
    def apply(self, p5obj, cell, sx0, sy0, sw, sh):
        sx = p5obj.remap(cell.umap1, 0, 1, sx0, sx0 + sw)
        sy = p5obj.remap(cell.umap2, 0, 1, sy0, sy0 + sh)
        
        # Transformation code A
#         s = np.array([sx, sy, 1])
#         t = np.matmul(self.tmat, s)
#         tx = t[0]
#         ty = t[1]

        # Transformation code B
        tx = self.multX(sx, sy)
        ty = self.multY(sx, sy)
        
        cell.selected = 0 <= tx and tx <= self.w and -self.h/2 <= ty and ty <= self.h/2        
        
    def updateBox(self, x, y):
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
  
    def displayBox(self, py5obj):
        py5obj.stroke(240, 118, 104)
        py5obj.no_fill()
        py5obj.push_matrix()
        py5obj.translate(self.x0, self.y0)
        py5obj.rotate(self.angle)
        py5obj.rect(0, -self.h/2, self.w, self.h)
        py5obj.stroke(240, 118, 104, 120)
        py5obj.line(0, 0, self.w, 0)
        py5obj.line(0.9 * self.w, -self.h/2, self.w, 0)
        py5obj.line(0.9 * self.w, +self.h/2, self.w, 0)
        py5obj.pop_matrix()
        
    def press(self, x, y, w, h):
        self.screen_width = w
        self.screen_height = h
        if self.state == CLOSED or self.state == COMPLETED:
            self.state = SET_SPINE
        elif self.state == SET_SPINE:
            self.state = SET_WIDTH

        if self.state == SET_SPINE:
            self.spx0 = self.spx1 = self.wx = x
            self.spy0 = self.spy1 = self.wy = y
            self.angle = 0
            self.w = self.h = 0
  
    def drag(self, x, y, w, h):
        if self.state == SET_SPINE:
            self.screen_width = w
            self.screen_height = h            
            self.spx1 = x
            self.spy1 = y  

    def move(self, x, y, w, h):
        if self.state == SET_WIDTH:
            self.screen_width = w
            self.screen_height = h            
            self.updateBox(x, y)

    def release(self, x, y, w, h):
        requestSelection = False
        if self.state == SET_SPINE:
            self.screen_width = w
            self.screen_height = h
            self.spx1 = x
            self.spy1 = y
            if self.spx1 != self.spx0 or self.spy1 != self.spy0:
                self.state = SET_WIDTH
                self.wx = self.spx1
                self.wy = self.spy1
            else:
                self.state = CLOSED
        elif self.state == SET_WIDTH:
            self.screen_width = w
            self.screen_height = h
            self.updateBox(x, y)
            self.state = COMPLETED
            requestSelection = True
        return(requestSelection)
            
    def multX(self, x, y):
        return self.tmat[0][0] * x + self.tmat[0][1] * y + self.tmat[0][2]
            
    def multY(self, x, y):
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