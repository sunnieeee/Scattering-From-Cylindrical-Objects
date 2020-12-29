# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 20:34:19 2020

@author: David Lyu
"""
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 15:08:02 2020
TFSF test
@author: lenovo
"""

import numpy as np
from matplotlib import pyplot as plt
lt=130 #时间长度
lx=100 
ly=100  #空间大小
lyy=36 #PML层
lxx=36
lx+=lxx
ly+=lyy

dx=0.1
dy=0.1  #空间步长
dt=0.1  #时间步长
wl=25
#free space的介电常数介磁常数
mu0=1   
ep0=2
sigma1=2        #PML内σ
sigmam1=sigma1*mu0/ep0
mu=np.zeros((lx,ly))
ep=np.zeros((lx,ly))
sigma=np.zeros((lx,ly))
sigmam=np.zeros((lx,ly))

#散射物体的位置与大小
x0=int(lx*3/5)
y0=int(ly*3/5)
r=8 #6格
#定义散射体内的介电常数
for i in range(0,lx):
    for j in range(0,ly):
        if (i-x0)**2+(j-y0)**2<=r*r:    #散射物体的参数
            mu[i,j]=1.5
            ep[i,j]=2.5
            sigma[i,j]=10
            sigmam[i,j]=10
            
        elif (lxx/2<i<lx-lxx/2)and(lyy/2<j<ly-lyy/2):    #free space的参数
            mu[i,j]=mu0
            ep[i,j]=ep0
            sigma[i,j]=0
            sigmam[i,j]=0
        else:          #PML的参数
            mu[i,j]=mu0
            ep[i,j]=ep0
            sigma[i,j]=sigma1
            sigmam[i,j]=sigmam1
            


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

#画图需要的mesh
X=np.arange(0,lx*dx,dx)
Y=np.arange(0,ly*dy,dy)
X, Y = np.meshgrid(X, Y)
#入射场原点的位置
xx=int(lx/2)        
yy=int(ly/2)


for i in range(0,lt-1):     #时间loop


    for j in range(0,lx-1):
        for k in range(0,ly-1):         #Hx与Hy的更新
            gm1=(mu[j,k]/dt-sigmam[j,k]/2)/(mu[j,k]/dt+sigmam[j,k]/2)
            gm2=1/(mu[j,k]/dt+sigmam[j,k]/2)
            Hx2[j,k]=gm1*Hx1[j,k]-gm2/(dy)*(Ez1[j,k]-Ez1[j,k-1])
            Hy2[j,k]=gm1*Hy1[j,k]+gm2/(dx)*(Ez1[j,k]-Ez1[j-1,k])

    for j in range(0,lx-1):
        for k in range(0,ly-1):    #Ez的更新
            g1=(ep[j,k]/dt-sigma[j,k]/2)/(ep[j,k]/dt+sigma[j,k]/2)
            g2=1/(ep[j,k]/dt+sigma[j,k]/2)

            Ez2[j,k]=g1*Ez1[j,k]+g2*((Hy2[j+1,k]-Hy2[j,k])/dy-(Hx2[j,k+1]-Hx2[j,k])/dx)
                
    #计算完入射场之后开始计算散射场
    
    for j in range(0,lx-1):
        for k in range(0,ly-1):         #散射场H
            if (lxx/2<j<lx-lxx/2)and(lyy/2<k<ly-lyy/2):
                alpha=1
            else:
                alpha=0
            gm1=mu[j,k]/dt+sigmam[j,k]/2
            gm11=(mu[j,k]-mu0)/dt+sigmam[j,k]/2
            gm2=mu[j,k]/dt-sigmam[j,k]/2
            gm22=(mu[j,k]-mu0)/dt-sigmam[j,k]/2
            Hxs2[j,k]=1/gm1*(gm2*Hxs1[j,k]-(Ezs1[j,k]-Ezs1[j,k-1])/dy+alpha*(-gm11*Hx2[j,k]+gm22*Hx1[j,k]))
            Hys2[j,k]=1/gm1*(gm2*Hys1[j,k]+(Ezs1[j,k]-Ezs1[j-1,k])/dx+alpha*(-gm11*Hy2[j,k]+gm22*Hy1[j,k]))

    for j in range(0,lx-1):
        for k in range(0,ly-1):    #Ez的更新（散射场）
            if (lxx/2<j<lx-lxx/2)and(lyy/2<k<ly-lyy/2):
                alpha=1
            else:
                alpha=0
            g1=ep[j,k]/dt+sigma[j,k]/2
            g11=(ep[j,k]-ep0)/dt+sigma[j,k]/2
            g2=ep[j,k]/dt-sigma[j,k]/2
            g22=(ep[j,k]-ep0)/dt-sigma[j,k]/2
            Ezs2[j,k]=1/g1*(g2*Ezs1[j,k]+((Hys2[j+1,k]-Hys2[j,k])/dy-(Hxs2[j,k+1]-Hxs2[j,k])/dx)+alpha*(-g11*Ez2[j,k]+g22*Ez1[j,k]))
    
    Ez2[xx,yy]=np.sin(np.pi*2*i/wl)  #源点不参与更新               
    Ez1[:,:]=Ez2[:,:]
    Hx1[:,:]=Hx2[:,:]
    Hy1[:,:]=Hy2[:,:]
    
    Ezs1[:,:]=Ezs2[:,:]
    Hxs1[:,:]=Hxs2[:,:]
    Hys1[:,:]=Hys2[:,:]
    

    if i==lt-2:
        
        """
        fig = plt.figure()  #定义新的三维坐标轴
        ax = plt.axes(projection='3d')
        ax.plot_surface(X,Y,Ezs2,cmap='rainbow')
        """
        plt.figure("TFSF")
        plt.imshow(np.log((Ezs2[int(lxx/2):lx-int(lxx/2),int(lyy/2):ly-int(lyy/2)])**2),cmap='gray_r',vmin=-6, vmax=0)
        
        #plt.imshow(np.log((Ez2)**2),cmap='gray_r',vmin=-6, vmax=0) #限定cbar的范围
        
        plt.colorbar()
        
