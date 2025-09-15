import math
import numpy as np
import matplotlib.pyplot as plt
import os

alpha = math.sqrt(2)-1
star1 = (alpha/2) % 0.5
star2 = star1 + 0.5
radius = 5

def arc(start, end, resolution = 20):
    if start > end:
        end += 1
    seq = np.linspace(start, end, resolution)
    x = radius*np.cos(2*math.pi*seq)
    y = radius*np.sin(2*math.pi*seq)
    return (x,y)

def geodesic(start, end, resolution = 20):
    if start > end:
        end += 1
    theta = (end - start)
    if abs(theta - 0.5) < 1e-9:
        x = np.linspace(radius*math.cos(2*math.pi*start), radius*math.cos(2*math.pi*end), resolution)
        y = np.linspace(radius*math.sin(2*math.pi*start), radius*math.sin(2*math.pi*end), resolution)
    else:
        R = radius / math.cos(math.pi*theta)
        r = radius * math.tan(math.pi*theta)
        seq = np.linspace(start+0.75, end+0.25, resolution)
        x = R*math.cos(2*math.pi*(start+theta/2)) + r*np.cos(2*math.pi*seq)
        y = R*math.sin(2*math.pi*(start+theta/2)) + r*np.sin(2*math.pi*seq)
    return (x,y)

def itinerary(x,n):
    resur = ''
    for i in range(n):
        if x > star1 and x < star2:
            resur += 'L'
        elif x == star1 or x == star2:
            resur += '*'
        else:
            resur += 'R'
        x = (x*2) % 1
    return resur

def match(w1, w2):
    if len(w1) > len(w2):
        temp = w1
        w1 = w2
        w2 = temp

    resur = True
    for i in range(len(w1)):
        if w1[i] != '*' and w2[i] != '*':
            if w1[i] != w2[i]:
                resur = False
    return resur

def legalize(w):
    n = len(w)
    lst = list(w)
    Ia = itinerary(alpha, n)
    for i in range(n-2,-1,-1):
        if match(lst[i+1:], Ia):
            lst[i] = '*'
    return ''.join(lst)

def duplicating(w, D=2, resolution=20):
    n = len(w)
    Ia = itinerary(alpha, n)
    intervals = []
    for i in range(n):
        temp = i
        for j in range(i+D, n):
            if match(w[i:j], Ia):
                temp = j
        if temp != i:
            intervals.append((i, temp))
    ax = plt.gca()
    plt.axis('off')
    ax.set_aspect(1)
    ax.set_ylim([0.4,0.6])
    ax.text(-0.05, 0.5, r'$W:$')
    ax.text(-0.07, 0.45, r'$I^\alpha(\alpha)$:')
    for i in range(n):
        ax.text(0.03* i-0.01, 0.5, w[i], size = 12)
        ax.text(0.03* i-0.01, 0.45, Ia[i], size = 12)
    for i,j in intervals:
        x = np.linspace(0.03*i, 0.03*(j-1), resolution)
        ax.plot(x, 0.07/(0.0009*(j-i-1)**2)*(0.03*i-x)*(x-0.03*(j-1))+0.53, linewidth=2, color = 'r')
    plt.show()
    

class GluingLink:
    def __init__(self, word = None, **kwargs):
        if not word:
            self.arcs = kwargs['arcs']
            self.word = kwargs['Word']
            self.leaves = kwargs['leaves']
        else:
            temp = GluingLink(arcs = [], Word = '', leaves = [])
            n = len(word)
            for i in range(n-1,-1,-1):
                if word[i] == 'L':
                    temp = temp.left()
                elif word[i] == 'R':
                    temp = temp.right()
                else:
                    temp = temp.split()
            self.arcs = temp.arcs
            self.word = temp.word
            self.leaves = temp.leaves

    def __str__(self):
        return f'Gluing Link: {self.word}\nLeaves: {' '.join(self.leaves)}'

    def arc_contain(self, arc, x):
        temp = False
        if arc[0] < arc[1]:
            if arc[0] <= x and arc[1] >= x:
                temp = True
        else:
            if arc[0] <= x or arc[1] >= x:
                temp = True
        return temp

    def contains(self, x):
        temp = False
        for arc in self.arcs:
            if self.arc_contain(arc,x):
                temp = True
        return temp

    def left(self):
        temp = []
        for arc in self.arcs:
            temp1 = arc[0]/2
            if temp1 < star1:
                temp1 += 0.5
            if temp1 > star2:
                temp1 -= 0.5
            temp2 = arc[1]/2
            if temp2 < star1:
                temp2 += 0.5
            if temp2 > star2:
                temp2 -= 0.5
            if temp1 > temp2:
                t = temp1
                temp1 = temp2
                temp2 = t
            if self.arc_contain(arc, alpha):
                temp.append((temp2, star2))
                temp.append((star1, temp1))
            else:
                temp.append((temp1, temp2))
        if not temp:
            temp = [(star1, star2)]
        Temp = []
        for leaf in self.leaves:
            Temp.append('L'+leaf)
        if (not Temp) or self.contains(alpha):
            Temp.append('*')
        return GluingLink(arcs = temp, Word = 'L'+self.word, leaves = Temp)

    def right(self):
        temp = []
        for arc in self.arcs:
            temp1 = arc[0]/2
            temp2 = arc[1]/2
            if temp1 < star1:
                temp1 += 1
            if temp2 < star1:
                temp2 += 1
            if temp1 < star2:
                temp1 += 0.5
            if temp2 < star2:
                temp2 += 0.5
            if temp1 > temp2:
                t = temp1
                temp1 = temp2
                temp2 = t
            if self.arc_contain(arc, alpha):
                temp.append((temp2%1, star1))
                temp.append((star2, temp1%1))
            else:
                temp.append((temp1, temp2))
        if not temp:
            temp = [(star2, star1)]
        Temp = []
        for leaf in self.leaves:
            Temp.append('R'+leaf)
        if (not Temp) or self.contains(alpha):
            Temp.append('*')
        return GluingLink(arcs = temp, Word = 'R'+self.word, leaves = Temp)

    def split(self):
        temp = []
        cut = -1
        n = len(self.arcs)
        for i in range(n):
            if self.arc_contain(self.arcs[i], alpha):
                cut = i
        if cut == -1:
            return GluingLink(arcs = [], Word = '', leaves = [])
        arc = self.arcs[cut]
        temp1 = arc[0]/2
        if temp1 < star1:
            temp1 += 0.5
        if temp1 > star2:
            temp1 -= 0.5
        temp2 = arc[1]/2
        if temp2 < star1:
            temp2 += 0.5
        if temp2 > star2:
            temp2 -= 0.5
        if temp1 > temp2:
            t = temp1
            temp1 = temp2
            temp2 = t

        temp3 = arc[0]/2
        temp4 = arc[1]/2
        if temp3 < star1:
            temp3 += 1
        if temp4 < star1:
            temp4 += 1
        if temp3 < star2:
            temp3 += 0.5
        if temp4 < star2:
            temp4 += 0.5
        if temp3 > temp4:
            t = temp3
            temp3 = temp4
            temp4 = t

        start = (temp4%1, temp1)
        end = (temp2, temp3%1)
        temp.append(start)
        for i in range(1, n):
            arc = self.arcs[(cut+i) % n]
            temp1 = arc[0]/2
            if temp1 < star1:
                temp1 += 0.5
            if temp1 > star2:
                temp1 -= 0.5
            temp2 = arc[1]/2
            if temp2 < star1:
                temp2 += 0.5
            if temp2 > star2:
                temp2 -= 0.5
            if temp1 > temp2:
                t = temp1
                temp1 = temp2
                temp2 = t
            temp.append((temp1, temp2))

        temp.append(end)
        for i in range(1, n):
            arc = self.arcs[(cut+i) % n]
            temp1 = arc[0]/2
            if temp1 > star1 and temp1 < star2:
                temp1 = (temp1 + 0.5) % 1
            temp2 = arc[1]/2
            if temp2 > star1 and temp2 < star2:
                temp2 = (temp2 + 0.5) % 1
            if temp1 > temp2:
                t = temp1
                temp1 = temp2
                temp2 = t
            temp.append((temp1, temp2))
        Temp = []
        for leaf in self.leaves:
            Temp.append('L'+leaf)
            Temp.append('R'+leaf)
        return GluingLink(arcs = temp, Word = '*'+self.word, leaves = Temp)

    def points(self):
        n = len(self.arcs)
        x = []
        y = []
        for i in range(n-1):
            Arc = arc(self.arcs[i][0], self.arcs[i][1])
            Geo = geodesic(self.arcs[i][1], self.arcs[i+1][0])
            x.extend(Arc[0])
            x.extend(Geo[0])
            y.extend(Arc[1])
            y.extend(Geo[1])
        Arc = arc(self.arcs[n-1][0], self.arcs[n-1][1])
        Geo = geodesic(self.arcs[n-1][1], self.arcs[0][0])
        x.extend(Arc[0])
        x.extend(Geo[0])
        y.extend(Arc[1])
        y.extend(Geo[1])
        return (x,y)

    def plot(self, save = False):
        circle = plt.Circle((0, 0), radius, color='b', fill=False)
        ax = plt.gca()
        ax.set_aspect('equal')
        plt.axis('off')
        l = self.points()
        ax.add_patch(circle)
        ax.plot(l[0],l[1], linewidth=2, c='orange')
        ax.plot([radius*math.cos(2*math.pi*star1),radius*math.cos(2*math.pi*star2)], [radius*math.sin(2*math.pi*star1),radius*math.sin(2*math.pi*star2)], linewidth=2, linestyle = '--', c='red')
        ax.plot(radius*math.cos(2*math.pi*alpha), radius*math.sin(2*math.pi*alpha), marker='o', markersize=4, c='r')
        if save:
            plt.savefig(os.getcwd() + '/Gluing_link.png', format="png", bbox_inches = 'tight')
        plt.show()

class GluingSet:
    def __init__(self, *args):
        self.links = []
        self.word = ''
        if len(args) == 1:
            self.add(args[0])

    def __str__(self):
        lst = []
        for link in self.links:
            lst.append(link.word)
        return f'Gluing Set: {self.word}\nGluing Links: {' '.join(lst)}'

    def left(self):
        if not self.links:
            self.links.append(GluingLink('L'))
        else:
            temp = []
            while self.links:
                link = self.links.pop()
                temp.append(link.left())
            self.links = temp
        self.word = 'L' + self.word

    def right(self):
        if not self.links:
            self.links.append(GluingLink('R'))
        else:
            temp = []
            while self.links:
                link = self.links.pop()
                temp.append(link.right())
            self.links = temp
        self.word = 'R' + self.word

    def star(self):
        temp = []
        while self.links:
            link = self.links.pop()
            if link.contains(alpha):
                temp.append(link.split())
            else:
                temp.append(link.left())
                temp.append(link.right())
        self.links = temp
        self.word = '*' + self.word

    def add(self, word):
        n = len(word)
        for i in range(n-1,-1,-1):
            if word[i] == 'L':
                self.left()
            elif word[i] == 'R':
                self.right()
            else:
                self.star()
    def plot(self, save = False):
        circle = plt.Circle((0, 0), radius, color='b', fill=False)
        ax = plt.gca()
        ax.set_aspect('equal')
        plt.axis('off')
        ax.add_patch(circle)
        ax.plot(radius*math.cos(2*math.pi*alpha), radius*math.sin(2*math.pi*alpha), marker='o', markersize=4, c='r')
        for link in self.links:
            l = link.points()
            ax.plot(l[0],l[1], linewidth=2, c='orange')
        ax.plot([radius*math.cos(2*math.pi*star1),radius*math.cos(2*math.pi*star2)], [radius*math.sin(2*math.pi*star1),radius*math.sin(2*math.pi*star2)], linewidth=2, linestyle = '--', c='red')
        if save:
            plt.savefig(os.getcwd() + '/Gluing_set.png', format="png", bbox_inches = 'tight')
        plt.show()

def GCS(word = None, **kwargs):
    if 'x' in kwargs:
        word = itinerary(kwargs['x'], kwargs['n'])
    return GluingLink(legalize(word))