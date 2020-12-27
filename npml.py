import numpy as np
import matplotlib.pyplot as plt

npml = 10
cst1 = 0.33
cst2 = 0.5
cst3 = 0.5 # dt/epsz
IE = 69 
JE = 69 # grid 70*70
c0 = 299792458.0


# TFSF
# set boundaries
iA = 11
iB = IE-iA
jA = 11
jB = JE-jA

ez_inc_low1 = 0
ez_inc_low2 = 0
ez_inc_high1 = 0
ez_inc_high2 = 0

#input source
inc = 20.0
spread = 6.0
lamb = 20e-6
freq = c0/lamb*10**(-6) #frequency in MHz

ez_inc = np.zeors(JE+1)
hx_inc = np.zeors(JE+1)

#ABC for incident buffer
ez_inc[0] = ez_inc_low2
ez_inc_low2 = ez_inc_low1
ez_inc_low1 = ez_inc[1]

ez_inc[JE] = ez_inc_high2
ez_inc_high2 = ez_inc_high1
ez_inc_low1 = ez_inc[JE-1]

#insert source
pulse = sin(2*np.pi)

# Ezinc & Hxinc
for j in range(1,JE+1):
    ez_inc[j] = ez_inc[j] + 0.5*(hx_inc[j-1]-hx_inc[j])
end
for j in range(JE):
    hx_inc[j] = hx_inc[j] + 0.5*(ez_inc[j]-ez_inc[j+1])
end



# incident Hx
for i in range(iA,iB+1):
    Hx[i][jA-1] = Hx[i][jA-1] + 0.5*ez_inc[jA] 
    Hx[i][jB] = Hx[i][jB] - 0.5*ez_inc[jB]
end

#CYLINDER PARAM
rad = 5.0
icenter = IE-1/2
jcenter = JE-1/2
sigma = 2.0
epsilon = 600.0
mur = np.ones(IE,JE)
gb = np.zeros(IE,JE)
iz = np.zeors(IE,JE)

for i in range(IE):
    for j in range(JE):
        xdist = i-icenter
        ydist = j-jcenter
        dist = np.sqrt(xdist**2+ydist**2)
        if dist<=rad :
            ga[i][j] = 1.0/(epsilon+sigma*cst3)
            gb[i][j] = sigma*cst3
        end
    end
end

# x direction
gi2 = np.zeros(IE+1)
gi3 = np.zeros(IE+1)
fi1 = np.zeros(IE+1)
fi2 = np.zeros(IE+1)
fi3 = np.zeros(IE+1)

for i in range(npml):
    xnum = npml - i #格点i到pml内边界距离
    xd = npml #x方向pml区域的厚度
    xxn = xnum / xd
    xn = cst1 * xxn**3
    gi2[i] = 1.0/(1.0+xn) #ｘ方向左pml
    gi3[i] = (1.0-xn)/(1.0+xn)
    gi2[IE-i] = 1.0/(1.0+xn) #x方向右pml
    gi3[IE-i] = (1.0-xn)/(1.0+xn)

    xxn = (xnum-0.5)/xd #i+1/2
    xn = cst1 * xxn**3
    fi1[i] = xn
    fi2[i] = 1.0/(1.0+xn)
    fi3[i] = (1.0-xn)/(1.0+xn)
    fi1[IE-i] = xn
    fi2[IE-i] = 1.0/(1.0+xn)
    fi3[IE-i] = (1.0-xn)/(1.0+xn)
end

# y direction
gj2 = np.zeros(JE+1)
gj3 = np.zeros(JE+1)
fj1 = np.zeros(JE+1)
fj2 = np.zeros(JE+1)
fj3 = np.zeros(JE+1)

for i in range(npml):
    ynum = npml - i #格点j到pml内边界距离
    yd = npml #y方向pml区域的厚度
    yyn = ynum / yd
    yn = cst1 * yyn**3
    gj2[i] = 1.0/(1.0+yn) #y方向上pml
    gj3[i] = (1.0-yn)/(1.0+yn)
    gj2[JE-i] = 1.0/(1.0+yn) #y方向下pml
    gj3[JE-i] = (1.0-yn)/(1.0+yn)

    yyn = (ynum-0.5)/yd #j+1/2
    yn = cst1 * yyn**3
    fj1[i] = yn
    fj2[i] = 1.0/(1.0+yn)
    fj3[i] = (1.0-yn)/(1.0+yn)
    fj1[JE-i] = yn
    fj2[JE-i] = 1.0/(1.0+yn)
    fj3[JE-i] = (1.0-yn)/(1.0+yn)
end

Dz = np.zeros(IE+1.JE+1)
Ez = np.zeros(IE+1.JE+1)
Hx = np.zeros(IE+1.JE+1)
Hy = np.zeros(IE+1.JE+1)
ihx = np.zeros(IE+1,JE+1)
ihy = np.zeros(IE+1,JE+1)
ga = np.zeros(IE+1,JE+1)
mu = np.zeros(IE+1,JE+1)

# Calculate Dz
for j in range(1,JE-1):
    for i in range(1,IE-1):
        Dz[i][j] = gi3[i]*gj3[j]*Dz[i][j] + gi2[i]*gj2[j]*cst2*(Hy[i][j]-Hy[i-1][j]-Hx[i][j]+Hx[i][j-1])
    end
end

## Calculate Ez
#for j in range(1,JE-1):
#    for i in range(1,IE-1):
#        ga[i][j] = 1.0
#        Ez[i][j] = ga[i][j]*Dz[i][j]
#    end
#end

# Calculate Hx
for j in range(1,JE-2):
    for i in range(1,IE-1):
        mu[i][j] = 1.0
        dele = Ez[i][j]-Ez[i][j+1]
        ihx[i][j] = ihx[i][j] + fi1[i]*dele
        Hx[i][j] = fj3[j]*Hx[i][j] + fj2[j]*cst2/mu[i][j]*(dele+ihx[i][j])
    end
end

# Calculate Hy
for j in range(1,JE-1):
    for i in range(1,IE-2):
        mu[i][j] = 1.0
       dele =  Ez[i+1][j]-Ez[i][j]
       ihy[i][j] = ihy[i][j] + fj1[i]*dele
       Hy[i][j] = fi3[i]*Hy[i][j] + fi2[i]*cst2/mu[i][j]*(dele+ihy[i][j])
    end
end

