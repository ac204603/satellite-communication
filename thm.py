# -*- coding: utf-8 -*-
import xlrd
import numpy as np
from thb import Beam
from geopy.distance import geodesic


if __name__ == '__main__':
    
    chu = [37.5,137]
    W = 12  #総電力
    Hz = 35*(10**6)
    f = []
    n = 3
    if n == 1:
        f.append(2.5 * 10 ** 9)
    elif n > 1:
        for i in np.linspace((2.5 * 10 ** 9),(2.535 * 10 ** 9),n):
            f.append(i)
    b = Beam()
    for i in range(len(f)):
        b.addFreq(f[i])
    
    wb=xlrd.open_workbook('curcle.xlsx')
    hsheet = wb.sheet_by_name('180')  #左のシートを選択
    oo = hsheet.nrows -1           #座標の数を取得
    keido,ido,x,y = [],[],[],[]
    
    for i in range(oo):
        keido.append(hsheet.cell(i + 1,0).value)
        ido.append(hsheet.cell(i + 1,1).value)
        
        idokyo = geodesic(chu,(ido[i],chu[1])).km
        keikyo = geodesic(chu,(chu[0],keido[i])).km
        if keido[i] < chu[1]:
            if ido[i] < chu[0]:
                x.append(-keikyo)
                y.append(-idokyo)
            else:
                x.append(-keikyo)
                y.append(idokyo)
        else:
            if ido[i] < chu[0]:
                x.append(keikyo)
                y.append(-idokyo)
            else:
                x.append(keikyo)
                y.append(idokyo)
        for _ in range(n):
            if i % n == _:
                b.addBeam(f[_], (x[i], y[i]))
                break
            else:
                pass
    
    
    # パラメータの初期値を指定（リスト作成）
    params = {}
    for i in range(oo):
        params[i] = W /oo
    
    haba = {}
    for i in range(len(f)):
        haba[i] = Hz / len(f)
        
    freq = dict()
    for _ in f:
        freq[_] = list()
    for i in range(oo):
        for _ in range(n):
            if i % n == _:
                freq[f[_]].append(i)
                break
            else:
                pass
    
    # 勾配の初期値を指定（リスト作成）
    grads = {}
    for i in range(oo):
        grads[i] = 0
    
    grahb = {}
    for i in range(len(f)):
        grahb[i] = 0


    ol = 0
    for key in params.keys():
        ol += params[key]
        print(key,":",params[key],"[W]")
    print("総電力：",ol,"[W]")
    print("")
    oi = 0
    for key in haba.keys():
        oi += haba[key]
        print(f[key],":",haba[key],"[Hz]")
    print("総周波数帯域：",oi,"[Hz]")
    print("")
    b.Sowa(params,freq,haba,1)
    print("総和：",b.Wa(0,2))