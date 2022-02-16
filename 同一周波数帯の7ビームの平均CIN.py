import os
import math
import itertools
import numpy as np
import scipy.special
from matplotlib import cm
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d


c = 299792458 # 光速[m/s]
Ds = 35786 # 静止衛星軌道までの距離[km]
R = 20 # 衛星の半径[m]
figfile = "figure"

btm = (-1573, -1393)
top = (1573, 1393)
spc = (500, 400)
baseX = np.linspace(btm[0], top[0], spc[0]) # 東西方向にspc[0]のポイント数だけサンプルする．
baseY = np.linspace(btm[1], top[1], spc[1]) # 南北方向にspc[1]のポイント数だけサンプルする．
freqList = list() # 周波数の値を入れるリスト
freqField = dict()# 角周波数ごとにビームを入れる辞書



def genBeam(freq,p):
    lmd = float(c) / freq # 波長lambda, 波長=光速/周波数であるから．
    x = np.linspace(btm[0] - p[0], top[0] - p[0], spc[0]) # 位置pにspc[0]のポイント数だけサンプルする．
    y = np.linspace(btm[1] - p[1], top[1] - p[1], spc[1]) # 位置pにspc[1]のポイント数だけサンプルする．
    X,Y = np.meshgrid(x, y) # np.meshgrid()は，一次元配列2個を使って，座標をメッシュ状に展開する．
    t = np.arctan2(np.sqrt(Y ** 2 + X ** 2), Ds) # 詳しくは中平先生の論文を参照．ビームのゲインを求めるために地上でのxy座標系から曲座標系に変換．
    s = (np.pi * R) / lmd * np.sin(t)
    return (2*scipy.special.jv(1, s) / s) ** 2 *1.5289*10**(-12)


def findex(freq):
    return freqList.index(freq)

def num(freq):
    return len(freqField[findex(freq)])

def addFreq(freq):
    if freq in freqList:
        print("{}[Hz] is exist".format(freq))
        return findex(freq)
    freqList.append(freq)
    return len(freqList) - 1


def dB(b):
    return 10 * np.log10(b)


def addBeam(freq, p):
    index = findex(freq)
    if freqField.get(index) is None:
        freqField[index] = list()
    gb = genBeam(freq, p)
    freqField[index].append((gb, p))
    return gb


def calcCI(freq,no):
    index = findex(freq)
    C = freqField[index][no][0] # CI比を知りたいビームのゲインを取得．
    I = np.zeros_like(C)
    for i in range(len(freqField[index])):
        if i != no:
            I += freqField[index][i][0]# CI比を知りたいビーム以外のビームをノイズとして加算しまくる．
    return dB(C/ (I+0.523783*10**(-13)))


def calcCIMean(ci, freq, no, r, plot = False):
    index = findex(freq)
    p = freqField[index][no][1]
    xmin, xmax = 0, 0
    for x in baseX:
        if x < p[0] - r:
            xmin += 1
        if x < p[0] + r:
            xmax += 1
    ymin, ymax = 0, 0
    for y in baseY:
        if y < p[1] - r:
            ymin += 1
        if y < p[1] + r:
            ymax += 1
            
    section = ci[int(ymin):int(ymax), int(xmin):int(xmax)]
    if plot:
        BX, BY = np.meshgrid(baseX,baseY)
        fig = plt.figure(figsize = (10.24, 7.68))
        ax = fig.add_subplot(111)
        plt.pcolormesh(BX, BY, ci, cmap="gist_ncar")
        pp = plt.colorbar(orientation='vertical')
        x1, x2 = p[0] - r, p[0] + r
        y1, y2 = p[1] - r, p[1] + r
        poly = plt.Polygon(((x1, y1), (x1, y2), (x2, y2), (x2, y1)), fill=False)
        ax.add_patch(poly)
        plt.xlabel("x[km]")
        plt.ylabel("y[km]")
        #plt.savefig(os.path.join(figfile, "main2.png"))
    
    return section.mean()



f1 = 2.5 * 10 ** 9
r = 112.5 # -3[dB] radius

addFreq(f1)

addBeam(f1, (0,0))
for deg in range(0, 360, 60):
    dx = r * math.cos(math.radians(deg))
    dy = r * math.sin(math.radians(deg))
    
    #if deg % 120 == 0:
    #    addBeam(f1, (dx * 2, dy * 2))
    if deg % 60 == 0:
        addBeam(f1, (dx * 2, dy * 2))
        

for i in range(num(f1)):
     ci = calcCIMean(calcCI(f1, i), f1, i, r, plot = True)
     print(i+1,"の平均C/(N+I)：", ci)