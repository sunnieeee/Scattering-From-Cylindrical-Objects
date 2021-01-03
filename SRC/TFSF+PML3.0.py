# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 15:40:12 2020

@author: David Lyu
"""
# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 15:08:02 2020
TFSF test
@author: lenovo
"""

import numpy as np
from matplotlib import pyplot as plt
import time
start_time = time.time()
lt=380 #时间长度
lx=70 
ly=70  #空间大小
lyy=36 #PML层
lxx=36
lx+=lxx
ly+=lyy

dx=0.1
dy=0.1  #空间步长 μm
dt=0.1  #时间步长 fs
wl=35  #m
#free space的介电常数介磁常数
mu0=1.25
ep0=8.85
sigma=2        #PML内σ大小
sigmam=sigma*mu0/ep0
#入射场的参数
mu1=np.zeros((lx,ly))
ep1=np.zeros((lx,ly))
sigma1=np.zeros((lx,ly))
sigmam1=np.zeros((lx,ly))
#散射场的参数
mu2=np.zeros((lx,ly))
ep2=np.zeros((lx,ly))
sigma2=np.zeros((lx,ly))
sigmam2=np.zeros((lx,ly))
#散射物体的位置与大小
x0=int(lx*3/5)
y0=int(ly*3/5)
r=6 #散射圆柱体半径
#定义散射体内的介电常数
mu1[:,:]=mu0
ep1[:,:]=ep0
mu2[:,:]=mu0
ep2[:,:]=ep0
for i in range(0,lx):
    for j in range(0,ly):
        if (i-x0)**2+(j-y0)**2<=r*r:    #散射物体的参数（只考虑散射场的）在这里改变物体的光学参数！！！！
            mu2[i,j]=2
            ep2[i,j]=10
            sigma2[i,j]=20
            sigmam2[i,j]=20
            
        elif (lxx/2<i<lx-lxx/2)and(lyy/2<j<ly-lyy/2):    #free space的参数（散射场入射场都一样）
            pass
        else:          #PML的参数（散射场入射场都一样）
            sigma1[i,j]=sigma
            sigmam1[i,j]=sigmam
            
            sigma2[i,j]=sigma
            sigmam2[i,j]=sigmam

#定义E与H scattered field & incident field
#入射场 follow原来的update equation
Hx1=np.zeros((lx,ly))   
Hx2=np.zeros((lx,ly))

Hy1=np.zeros((lx,ly))
Hy2=np.zeros((lx,ly))

Ez1=np.zeros((lx,ly))
Ez2=np.zeros((lx,ly))

#散射场 follow新的update equation
Hxs1=np.zeros((lx,ly))   
Hxs2=np.zeros((lx,ly))

Hys1=np.zeros((lx,ly))
Hys2=np.zeros((lx,ly))

Ezs1=np.zeros((lx,ly))
Ezs2=np.zeros((lx,ly))

#迭代需要的系数

gm1=(mu1/dt-sigmam1/2)/(mu1/dt+sigmam1/2)
gm2=1/(mu1/dt+sigmam1/2)
g1=(ep1/dt-sigma1/2)/(ep1/dt+sigma1/2)
g2=1/(ep1/dt+sigma1/2)

ggm1=mu2/dt+sigmam2/2
ggm2=mu2/dt-sigmam2/2
gg1=ep2/dt+sigma2/2
gg2=ep2/dt-sigma2/2
#入射场原点的位置
xx=int(lx*2/5)        
yy=int(ly*2/5)

#计算
for i in range(0,lt-1):     #时间loop

    for j in range(0,lx-1):
        for k in range(0,ly-1):         #入射场Hx与Hy的更新
            Hx2[j,k]=gm1[j,k]*Hx1[j,k]-gm2[j,k]/(dy)*(Ez1[j,k]-Ez1[j,k-1])
            Hy2[j,k]=gm1[j,k]*Hy1[j,k]+gm2[j,k]/(dx)*(Ez1[j,k]-Ez1[j-1,k])

    for j in range(0,lx-1):
        for k in range(0,ly-1):    #入射场Ez的更新
            Ez2[j,k]=g1[j,k]*Ez1[j,k]+g2[j,k]*((Hy2[j+1,k]-Hy2[j,k])/dy-(Hx2[j,k+1]-Hx2[j,k])/dx)
                
    #计算完入射场之后开始计算散射场
    
    for j in range(0,lx-1):
        for k in range(0,ly-1):         #散射场H
            if (lxx/2<j<lx-lxx/2)and(lyy/2<k<ly-lyy/2):     #这个alpha的引入是为了在pml内去掉下方👇迭代式中入射场Hy,Hx的影响
                alpha=1
            else:
                alpha=0
            Hxs2[j,k]=1/ggm1[j,k]*(ggm2[j,k]*Hxs1[j,k]-(Ezs1[j,k]-Ezs1[j,k-1])/dy+alpha*(-(ggm1[j,k]-mu0/dt)*Hx2[j,k]+(ggm2[j,k]-mu0/dt)*Hx1[j,k]))
            Hys2[j,k]=1/ggm1[j,k]*(ggm2[j,k]*Hys1[j,k]+(Ezs1[j,k]-Ezs1[j-1,k])/dx+alpha*(-(ggm1[j,k]-mu0/dt)*Hy2[j,k]+(ggm2[j,k]-mu0/dt)*Hy1[j,k]))

    for j in range(0,lx-1):
        for k in range(0,ly-1):    #Ez的更新（散射场）
            if (lxx/2<j<lx-lxx/2)and(lyy/2<k<ly-lyy/2):     #这个alpha同理
                alpha=1
            else:
                alpha=0

            Ezs2[j,k]=1/gg1[j,k]*(gg2[j,k]*Ezs1[j,k]+((Hys2[j+1,k]-Hys2[j,k])/dy-(Hxs2[j,k+1]-Hxs2[j,k])/dx)+alpha*(-(gg1[j,k]-ep0/dt)*Ez2[j,k]+(gg2[j,k]-ep0/dt)*Ez1[j,k]))
    
    Ez2[xx,yy]=np.sin(np.pi*2*i/wl)  #源点不参与更新               
    Ez1[:,:]=Ez2[:,:]
    Hx1[:,:]=Hx2[:,:]
    Hy1[:,:]=Hy2[:,:]
    
    Ezs1[:,:]=Ezs2[:,:]
    Hxs1[:,:]=Hxs2[:,:]
    Hys1[:,:]=Hys2[:,:]
    

    
        
        
plt.figure("TFSF+")
plt.subplot(2,2,1)
plt.imshow(np.log((Ezs2[int(lxx/2):lx-int(lxx/2),int(lyy/2):ly-int(lyy/2)]+Ez2[int(lxx/2):lx-int(lxx/2),int(lyy/2):ly-int(lyy/2)])**2),cmap='gray_r',vmin=-13, vmax=0)
plt.title('total field')
plt.subplot(2,2,2)
plt.imshow(np.log((Ezs2[int(lxx/2):lx-int(lxx/2),int(lyy/2):ly-int(lyy/2)])**2),cmap='gray_r',vmin=-13, vmax=0)        
plt.title('scattered field')
plt.subplot(2,2,3)
plt.imshow(np.log((Ez2[int(lxx/2):lx-int(lxx/2),int(lyy/2):ly-int(lyy/2)])**2),cmap='gray_r',vmin=-13, vmax=0)
plt.title('incident field')        #plt.imshow(np.log((Ez2)**2),cmap='gray_r',vmin=-6, vmax=0) #限定cbar的范围
plt.subplot(2,2,4)
plt.imshow((ep2[int(lxx/2):lx-int(lxx/2),int(lyy/2):ly-int(lyy/2)])**2,cmap='gray_r')
plt.title('scatter object')          
plt.colorbar()

end_time = time.time()
print('running time: ',end_time-start_time)

