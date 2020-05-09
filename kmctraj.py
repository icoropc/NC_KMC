# -*- coding: utf-8 -*-
"""
Created on Wed May 24 13:44:28 2017

@author: Igor
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 17:13:15 2013

@author: Igor
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

from mpl_toolkits.mplot3d import Axes3D

#sns.set()
#sns.set_style("white")
sns.set_style("ticks")
#sns.set_palette("BuGn_r")

#data0=pd.read_csv('out_2',index_col=0,skiprows=0)
data0=np.genfromtxt('out_2')


timevals=[]
xvals=[]
yvals=[]
zvals=[]
speed=0
disttrav=0

for i in data0:
    if i[0]==2:
        timevals.append(i[1])
        xvals.append(i[4])
        yvals.append(i[5])
        zvals.append(i[6])
        if (i[1] != 0.0):
            disttrav+=i[3]
            speed=disttrav/i[1]

#plt.plot(xvals,yvals,'k-')
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(xvals,yvals,zvals,"k")        
P=ax.scatter3D(xvals,yvals,zvals,c=timevals, cmap='RdBu_r')        

fig.colorbar(P,shrink=0.7, aspect=20,label='Time (ps)')

ax.set_xlabel("X (nm)")
ax.set_ylabel("Y (nm)")
ax.set_zlabel("Z (nm)")


#plt.plot(datanorm2.iloc[:,0])


#
#ax = plt.gca()
#ax.set_yscale('log')
#ax.set_xscale('log')
#
#
#plt.xlabel(r'q ($\AA^{-1}$)',fontsize=12)
#plt.ylabel('Intensity (AU)',fontsize=12)
##plt.legend(["1 min","3min","2nd add", "later", "redisp"], fontsize=15)
##plt.legend(["all processed in air", "gb proc 0.25M lig sol", "good"], fontsize=15)
##plt.legend(["0s", "3min", "10min","120min"], fontsize=15)
#
#
#plt.xlim(0.02,1)
#plt.ylim(0.006,15)
#
plt.tight_layout()