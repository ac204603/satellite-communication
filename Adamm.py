# -*- coding: utf-8 -*-
import xlrd
import numpy as np
from beam import Beam
from beam import Adam
import matplotlib.pyplot as plt
from geopy.distance import geodesic


if __name__ == '__main__':
    
    chu = [37.5,137]    #中心
    W = 12  #総電力
    Hz = 35*(10**6)     #中心周波数
    f = []              #周波数帯域を格納する配列
    n = 1               #繰り返しビーム数
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
    
    #緯度経度を距離に変換
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

    
    # 各ビームに割り当てる電力の初期値を指定（リスト作成）
    params = {}
    for i in range(oo):
        params[i] = W /oo
    
    # 各周波数帯域の帯域幅を指定
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
    
    # 学習率を指定
    lr = 0.1
    
    # 減衰率を指定
    beta1 = 0.999
    beta2 = 0.9
    
    # 試行回数を指定
    iter_num = 3
    Sowa,Hensa,Z = [],[],[]
    xziku = np.array([x for x in range(iter_num)])
    xziku += 1
    # インスタンスを作成
    optimizer = Adam(lr=lr, beta1=beta1, beta2=beta2)


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
    print("--")
    jinko = b.Per(freq)
    b.Sowa(params,freq,haba,jinko,1)
    print("総和：",b.Wa(0,2),"[bps]")
    print("--")
    print("")
    print("")
    z = b.Sowa(params,freq,haba,jinko,2)
    for _ in range(iter_num):
        # パラメータを更新
        grads = b.df(params,freq,haba,1e-10,jinko,3)
        optimizer.updateW(params,grads,W)
        b.Sowa(params,freq,haba,jinko,2)
        Sowa.append(b.Wa(0,2))
        if z < b.Sowa(params,freq,haba,jinko,2):
            z = b.Sowa(params,freq,haba,jinko,2)
            print("試行回数:",_ + 1)
            ol = 0
            for key in params.keys():
                ol += params[key]
                print(key,":",params[key],"[W]")
            print("総電力",ol,"[W]")
            print("")
            oi = 0
            for key in haba.keys():
                oi += haba[key]
                print(f[key],":",haba[key],"[Hz]")
            print("総周波数帯域：",oi,"[Hz]")
            print("")
            b.Sowa(params,freq,haba,jinko,1)
            print("総和：",b.Wa(0,2),"[bps]")
            print("--")
            print("")
        else:
            print("試行回数:",_ + 1)

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(1,1,1)
    ax2.plot(xziku,Sowa)
    ax2.set_xlabel("iterration")
    ax2.set_ylabel("Sum(bps)")
    print("fin")