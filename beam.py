# -*- coding: utf-8 -*-
import xlrd
import numpy as np
import scipy.special
from geopy.distance import geodesic
from japanmap import pref_names
from japanmap import get_data, pref_points


class Beam:
    c = 299792458 # 光速[m/s]
    Ds = 35786 # 静止衛星軌道までの距離[km]
    R = 20 # 衛星の直径[m]
    r = 180 #ビームの半径[km]
    
    def __init__(self,btm = (-900, -900),top = (900, 900),spc = (90, 90)):
        self.btm = btm
        self.top = top
        self.spc = spc
        self.baseX = np.linspace(btm[0], top[0], spc[0])
        self.baseY = np.linspace(btm[1], top[1], spc[1])
        self.freqList = list() # 周波数の値を入れるリスト
        self.freqField = dict() # 角周波数ごとにビームを入れる辞書
        self.Hensa = 0
        self.Gp = 10**(-90/10)  #伝搬利得Gp
        
        
        self.oll = pref_names
        #self.oll.remove("_")
        qpqo = get_data()
        pnts = pref_points(qpqo)
        self.ken = dict()
        for i in self.oll:
            self.ken[i] = pnts[self.oll.index(i)]
        self.chu = [37.5,137]
        self.shin = dict()
        
        #緯度経度を距離に変換
        for i in self.ken.keys():
            x = []
            for j in range(len(self.ken[i])):
            
                idokyo = geodesic(self.chu,(self.ken[i][j][1],self.chu[1])).km
                keikyo = geodesic(self.chu,(self.chu[0],self.ken[i][j][0])).km
                if self.ken[i][j][0] < self.chu[1]:
                    if self.ken[i][j][1] < self.chu[0]:
                        x.append([-keikyo,-idokyo])
                    else:
                        x.append([-keikyo,idokyo])
                else:
                    if self.ken[i][j][1] < self.chu[0]:
                        x.append([keikyo,-idokyo])
                    else:
                        x.append([keikyo,idokyo])
            self.shin[i] = x
        

    #ビーム利得の算出
    def genBeam(self,freq,p):
        lmd = float(self.c) / freq # 波長lambda, 波長=光速/周波数であるから．
        x = np.linspace(self.btm[0] - p[0], self.top[0] - p[0], self.spc[0]) # 位置pにspc[0]のポイント数だけサンプルする．
        y = np.linspace(self.btm[1] - p[1], self.top[1] - p[1], self.spc[1]) # 位置pにspc[1]のポイント数だけサンプルする．
        self.X, self.Y = np.meshgrid(x, y) # np.meshgrid()は，一次元配列2個を使って，座標をメッシュ状に展開する．
        t = np.arctan2(np.sqrt(self.Y ** 2 + self.X ** 2), self.Ds) # 詳しくは中平先生の論文を参照．ビームのゲインを求めるために地上でのxy座標系から曲座標系に変換．
        s = (np.pi * self.R) / lmd * np.sin(t) # 詳しくは中平先生の論文を参照．ビームのゲインを求めるために，式の共通項を求める．
        return (2*scipy.special.jv(1, s) / s) ** 2 # 詳しくは中平先生の論文を参照．ビームのゲインを求めて返す．scipy.special.jv(n, s)は第nベッセル関数．中平先生の論文にはベッセル関数の添字に1があったので，第一次ベッセル関数であると判断して本記述とした．

    #周波数帯域の探索
    def findex(self, freq):
        return self.freqList.index(freq)

    #周波数帯域を追加する
    def addFreq(self, freq):
        if freq in self.freqList:
            print("{}[Hz] is exist".format(freq))
            return self.findex(freq)
        else:
            return self.freqList.append(freq)
    

    # ビームのゲインをデシベルに変換する．
    def dB(self, b):
        return 10 * np.log10(b)

    #pの座標のビームを追加する
    def addBeam(self,freq,p):
        index = self.findex(freq)
        if self.freqField.get(index) is None:
            self.freqField[index] = list()
        gb = self.genBeam(freq, p)
        return self.freqField[index].append((gb, p))
    
    
    #最小二乗法を用いて、各ビームのC/(N+I)から周波数利用効率を導き、スループットを算出
    def BitRate(self,CNR,haba,reb):
        xk = np.array([1.4,2.7,4.6,5.6,7.4,10,11.2,13.8,15.8,19.5])
        yi = np.array([1,1.2,1.5,1.67,1.98,2.52,3,3.5,4,4.5])
        a = ((np.dot(xk, yi)-yi.sum()*xk.sum() / len(xk)) / ((xk ** 2).sum() - xk.sum() **2 / (len(xk))))
        b = (yi.sum() - a * xk.sum()) / len(xk)
        return (a * CNR + b)*(haba/reb)
    
    #全ユーザーのスループットから総和を算出
    def Wa(self,Sum,n):
        if n == 1:
            self.Hensa = np.sum(Sum)
            return 0
        elif n == 2:
            return self.Hensa


    def Per(self,freq):
        por = {}
        for j in self.freqList:
            index = self.findex(j)
            for no in range(len(self.freqField[index])):
                p = self.freqField[index][no][1]#ビームの座標を取得
                di = {}
                for i in self.shin.keys():
                    ol = 0
                    x = []
                    for l in range(len(self.shin[i])):
                        if np.sqrt((self.shin[i][l][0] - p[0])**2 + (self.shin[i][l][1] - p[1]) ** 2) <= self.r:
                            x.append([self.shin[i][l][0],self.shin[i][l][1]])
                            ol = 1
                    if ol == 1:
                        di[i] = x
                        souwa = 0
                        for k in range(len(di[i])):
                            if k == len(di[i]) - 1:
                                wa = abs(di[i][k][0] * di[i][0][1] - di[i][k][1] * di[i][0][0])
                            else:
                                wa = abs(di[i][k][0]*di[i][k + 1][1] - di[i][k][1]*di[i][k + 1][0])
                            souwa += wa
                        souwa = abs(souwa / 2)
                        di[i] = souwa
                por[freq[j][no]] = di
        for i in por.keys():
            if len(por[i]) == 1:
                for l in por[i].keys():
                    po = l
            for j in por.keys():
                if (po in por[j]) and (i != j) and (len(por[j]) != 1):
                    del por[j][po]
        for i in self.shin.keys():
            de = {}
            for j in por.keys():
                if (i in por[j]) and (len(por[j]) != 1):
                    de[j] = por[j][i]
            if len(de) > 1:
                mak = max(de,key = de.get)
                for k in de.keys():
                    if k != mak:
                        del por[k][i]
                        
        wb=xlrd.open_workbook('jinko.xls')
        hsheet = wb.sheet_by_name('j')
        oo = hsheet.nrows
        todo = {}
        for i in range(oo):
            todo[hsheet.cell(i,0).value] = hsheet.cell(i,1).value
        
        
        for i in self.shin.keys():
            de = []
            for j in por.keys():
                if i in por[j]:
                    de.append(j)
            if len(de) == 1:
                por[de[0]][i] = todo[i] / 125263
            elif len(de) > 1:
                for j in de:
                    por[j][i] = todo[i] * len(de) / 125263
        for i in por.keys():
            por[i] = sum(por[i].values())
        return por

    #サービスエリアに存在するユーザーのスループットから平均を算出
    def Sowa(self,params,freq,haba,jinko,oi):
        shi = list()
        for j in self.freqList:
            #ai = []
            co = 0
            poi = []
            index = self.findex(j)
            pi = []
            N = (1.38*10**(-23))*haba[index]*316.3
            for k in freq[j]:
                poi.append(params[k])
            for no in range(len(self.freqField[index])):
                pe = 0
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
                            pe += 1
                pi.append(pe)
            for no in range(len(self.freqField[index])):
                go = []
                C = self.freqField[index][no][0]*poi[no] # CI比を知りたいビームのゲインを取得．
                I = np.zeros_like(C)
                for i in range(len(self.freqField[index])):
                    if i != no:
                        I += self.freqField[index][i][0]*poi[i]
                ci = C/ (I+N / self.Gp)
                if oi == 1:
                    print(freq[j][no]," ",self.dB(np.average(C)))
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
                            shi.append(self.BitRate(self.dB(ci[self.baseY==y,self.baseX==x]),haba[index],co) *(jinko[freq[j][no]]))
                            go.append(self.BitRate(self.dB(ci[self.baseY==y,self.baseX==x]),haba[index],co) *(jinko[freq[j][no]]))
                #print(freq[j][no]," ",np.average(go))
        if oi == 1:
            print("")
            self.Wa(shi,oi)
        return np.average(shi)

    #微分
    def df(self,params,freq,haba,dx,jinko,n):
        if n == 3:
            n = 0
            he = {}
            for ke in params.keys():
                bol = params.copy()
                bol[ke] = bol[ke] + dx
                he[ke] = (-self.Sowa(bol,freq,haba,jinko,n) - (-self.Sowa(params,freq,haba,jinko,n))) / dx
            return he
        elif n == 4:
            n = 0
            he = {}
            for ke in haba.keys():
                bol = haba.copy()
                bol[ke] = bol[ke] + dx
                he[ke] = (-self.Sowa(params,freq,bol,jinko,n) - (-self.Sowa(params,freq,haba,jinko,n))) / dx
            return he

#最適化アルゴリズム
class Adam:
    # インスタンス変数を定義
    def __init__(self, lr=0.001, beta1=0.9, beta2=0.999):
        self.lr = lr # 学習率
        self.beta1 = beta1 # mの減衰率
        self.beta2 = beta2 # vの減衰率
        self.iter = 0 # 試行回数を初期化
        self.m = None # モーメンタム
        self.v = None # 適合的な学習係数
    
    # パラメータの更新メソッドを定義
    def updateW(self, params, grads,W):
        # mとvを初期化
        if self.m is None: # 初回のみ
            self.m = {}
            self.v = {}
            for key, val in params.items():
                self.m[key] = np.zeros_like(val) # 全ての要素が0
                self.v[key] = np.zeros_like(val) # 全ての要素が0
        
        # パラメータごとに値を更新
        self.iter += 1 # 更新回数をカウント
        lr_t  = self.lr * np.sqrt(1.0 - self.beta2 ** self.iter) / (1.0 - self.beta1 ** self.iter) 
        for key in params.keys():
            self.m[key] = self.beta1 * self.m[key] + (1 - self.beta1) * grads[key]
            self.v[key] = self.beta2 * self.v[key] + (1 - self.beta2) * (grads[key] ** 2)
            ss = params[key] - lr_t * self.m[key] / (np.sqrt(self.v[key]) + 1e-7)
            if ss < 0:
                params[key] *=  0.1
            else:
                params[key] = ss
        if sum(params.values()) < W:
            am = W - sum(params.values())
            for key in params.keys():                
                params[key] +=  am / len(params.values())
        elif sum(params.values()) > W:
            sa = sum(params.values()) - W
            for key in params.keys():
                params[key] -= (params[key] / sum(params.values())) * sa