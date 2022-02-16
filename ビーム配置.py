import xlrd
import shapely
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.ticker as tick # 目盛り操作に必要なライブラリを読み込みます
import matplotlib.patches as patches


fig,ax = plt.subplots(figsize = (13,13))


df = gpd.read_file('japan.geojson')
df.plot(ax = ax,figsize = (8,8),edgecolor='#444', facecolor='white', linewidth = 0.5,aspect="equal");

xmin = 125
xmax = 150
ymin = 30
ymax = 46

plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.gca().xaxis.set_minor_locator(tick.MultipleLocator(1))
plt.gca().yaxis.set_minor_locator(tick.MultipleLocator(1))

plt.grid(which = "major",alpha = 2.6)
plt.grid(which = "minor",alpha = 0.3)


wb=xlrd.open_workbook('curcle.xlsx')
hsheet = wb.sheet_by_name('180')  #左のシートを選択
re = hsheet.nrows -1           #座標の数を取得
gyo = hsheet.ncols
x = []
y = []
za = []

for i in range(re):
    x.append(hsheet.cell(i + 1,0).value)
    y.append(hsheet.cell(i + 1,1).value)
    tmp = 44 - y[i]
    ra = 2.25 - (tmp*3/80)
    za.append(patches.Circle(xy=(x[i],y[i]), radius=ra,color = (0,0,0,0.2)))
    ax.add_patch(za[i])

    
plt.show()
#120,000m
#fig.savefig("120_19.png")