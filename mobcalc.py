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
import os

#sns.set()
#sns.set_style("white")
sns.set_style("ticks")
#sns.set_palette("BuGn_r")


temp_list=[5,25, 50, 100, 150, 200, 250, 300,500,700,1000,2000,5000,10000]

timevals=[]
xvals=[]
yvals=[]
zvals=[]
speed=0
disttrav=0

disttrav_list = np.zeros(3000)
time_list = np.zeros(3000)
wait_c=0

mobil_list=[]

def run_calcmob(temp):
    cmnd='./ma_hop_bcc_3 < input_ma 2 '  + str(temp) + ' > out_2'
#    cmnd='./ma_hop_bcc_3_MH < input_ma 2 '  + str(temp) + ' > out_2'
    wait_c=os.system(cmnd)
    data0=np.genfromtxt('out_2')
    for i in data0:
        inx=int(i[0])
        if (i[1] != 0.0):
            disttrav_list[inx]+=i[3]
            time_list[inx]=i[1]
    speedlist=disttrav_list/time_list
    return(np.average(speedlist))

#tempmob=run_calcmob(300)
for temp in temp_list:
    mobil_list.append(run_calcmob(temp))

print(mobil_list)

plt.plot(temp_list,mobil_list)