# -*- coding: utf-8 -*-
import numpy as np
import scipy.special


class Beam:
    c = 299792458 # 光速[m/s]
    Ds = 35786 # 静止衛星軌道までの距離[km]
    R = 20 # 衛星の直径[m]
    r = 180
    
    def __init__(self,btm = (-900, -900),top = (900, 900),spc = (90, 90)):
        self.btm = btm
        self.top = top
        self.spc = spc
        self.baseX = np.linspace(btm[0], top[0], spc[0])
        self.baseY = np.linspace(btm[1], top[1], spc[1])
        self.freqList = list() # 周波数の値を入れるリスト
        self.freqField = dict() # 角周波数ごとにビームを入れる辞書
        self.Gokei = 0
        self.Hensa = 0
        self.Gp = 10**(-90/10)
        


    def genBeam(self,freq,p):
        lmd = float(self.c) / freq # 波長lambda, 波長=光速/周波数であるから．
        x = np.linspace(self.btm[0] - p[0], self.top[0] - p[0], self.spc[0]) # 位置pにspc[0]のポイント数だけサンプルする．
        y = np.linspace(self.btm[1] - p[1], self.top[1] - p[1], self.spc[1]) # 位置pにspc[1]のポイント数だけサンプルする．
        self.X, self.Y = np.meshgrid(x, y) # np.meshgrid()は，一次元配列2個を使って，座標をメッシュ状に展開する．
        t = np.arctan2(np.sqrt(self.Y ** 2 + self.X ** 2), self.Ds) # 詳しくは中平先生の論文を参照．ビームのゲインを求めるために地上でのxy座標系から曲座標系に変換．
        s = (np.pi * self.R) / lmd * np.sin(t) # 詳しくは中平先生の論文を参照．ビームのゲインを求めるために，式の共通項を求める．
        return (2*scipy.special.jv(1, s) / s) ** 2 # 詳しくは中平先生の論文を参照．ビームのゲインを求めて返す．scipy.special.jv(n, s)は第nベッセル関数．中平先生の論文にはベッセル関数の添字に1があったので，第一次ベッセル関数であると判断して本記述とした．


    def findex(self, freq):
        return self.freqList.index(freq)


    def addFreq(self, freq):
        if freq in self.freqList:
            print("{}[Hz] is exist".format(freq))
            return self.findex(freq)
        else:
            return self.freqList.append(freq)
    

    # ビームのゲインをデシベルに変換する．
    def dB(self, b):
        return 10 * np.log10(b)


    def addBeam(self,freq,p):
        index = self.findex(freq)
        if self.freqField.get(index) is None:
            self.freqField[index] = list()
        gb = self.genBeam(freq, p)
        return self.freqField[index].append((gb, p))
    
    
    def BitRate(self,CNR,haba,reb):
        xk = np.array([1.4,2.7,4.6,5.6,7.4,10,11.2,13.8,15.8,19.5])
        yi = np.array([1,1.2,1.5,1.67,1.98,2.52,3,3.5,4,4.5])
        a = ((np.dot(xk, yi)-yi.sum()*xk.sum() / len(xk)) / ((xk ** 2).sum() - xk.sum() **2 / (len(xk))))
        b = (yi.sum() - a * xk.sum()) / len(xk)
        return (a * CNR + b)*(haba/reb)
    
    
    def Wa(self,Sum,n):
        if n == 1:
            self.Hensa = np.sum(Sum)
            return 0
        elif n == 2:
            return self.Hensa


    def Sowa(self,params,freq,haba,oi):
        shi = list()
        for j in self.freqList:
            co = 0
            poi = []
            index = self.findex(j)
            N = (1.38*10**(-23))*haba[index]*316.3
            for k in freq[j]:
                poi.append(params[k])
            for no in range(len(self.freqField[index])):
                p = self.freqField[index][no][1]#ビームの座標を取得
                xmin, xmax = 0, 0
                ymin, ymax = 0, 0
                for x in self.baseX:
                    if x < p[0] - self.r:
                        xmin += 1
                    if x < p[0] + self.r:
                        xmax += 1
                for y in self.baseY:
                    if y < p[1] - self.r:
                        ymin += 1
                    if y < p[1] + self.r:
                        ymax += 1
                for x in self.baseX[xmin:xmax]:
                    for y in self.baseY[ymin:ymax]:
                        if np.sqrt((x - p[0]) ** 2 + (y - p[1]) ** 2) <= self.r:
                            co += 1
            for no in range(len(self.freqField[index])):
                C = self.freqField[index][no][0]*poi[no] # CI比を知りたいビームのゲインを取得．
                I = np.zeros_like(C)
                for i in range(len(self.freqField[index])):
                    if i != no:
                        I += self.freqField[index][i][0]*poi[i]
                if oi == 1:
                    print(freq[j][no]," ",np.average(I))
                ci = C/ (I+N / self.Gp)
                #ai.append(I)
                p = self.freqField[index][no][1]#ビームの座標を取得
                xmin, xmax = 0, 0
                ymin, ymax = 0, 0
                for x in self.baseX:
                    if x < p[0] - self.r:
                        xmin += 1
                    if x < p[0] + self.r:
                        xmax += 1
                for y in self.baseY:
                    if y < p[1] - self.r:
                        ymin += 1
                    if y < p[1] + self.r:
                        ymax += 1
                for x in self.baseX[xmin:xmax]:
                    for y in self.baseY[ymin:ymax]:
                        if np.sqrt((x - p[0]) ** 2 + (y - p[1]) ** 2) <= self.r:
                            shi.append(self.BitRate(self.dB(ci[self.baseY==y,self.baseX==x]),haba[index],co))
                #print(freq[j][no]," ",np.average(go))
        if oi == 1:
            print("")
            self.Wa(shi,oi)
        return np.average(shi)