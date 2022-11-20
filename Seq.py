import numpy as np
from scipy.integrate import *
import matplotlib.pyplot as plt
import math
from math import exp

K = 1
r0 = 0.05
betar = 1 #1,2
v0 = 1
betadelta = 1
delta0 = .02
zeta0 = .03*4 #*4
zeta02 = .01*4 #*4
betazeta = 1
betazeta2 = 5
amax = 2
betaa = 0.5 #0.5,0
vd = 2
ft = 0
k = .01
chemo_start = 1344 #8 weeks
chemo_end = chemo_start+2688 #24 weeks
inter_start = chemo_start+2688
inter_end = inter_start+2688
time=inter_end
min_summ = []

total_dose = 160 #160. Total dose and num doses doubled for combo bc giving for twice the length.
num_doses = 16 #16
txdose = 1

#arr = [0.43865851701724917,0.4315968168640311,0.4279078962932176,0.42565877251212814,0.42417010454285514,0.42305997772104065,0.42222228394801475,0.4215813398682476,0.42105995742555813]
#times = np.linspace(1,5,9)

IC = [.6,1,1]

#Measure total drug (1000) and look at hyper/hypo fractionation, keeping dose constant
#Tx interval vs dose given
#5000 time steps, 1 dose, 10 doses

def grow(v):
    return r0*(v/v0)**betar

def death(v):
    return delta0*math.exp(betadelta*(v-v0))

def chemo(v):
    return zeta0*math.exp(-betazeta*(v-v0))

def target(v):
    return zeta02*((v/v0)**betazeta2)

def fc(t):
    if t>=chemo_start and t<=chemo_end:
        txtime = total_dose/num_doses/txdose
        period = (chemo_end-chemo_start)/num_doses
        if (t-chemo_start)%period < txtime:
            return txdose
        else:
            return 0
    else:
        return 0
    
def ft(t):
    if t>=inter_start and t<=inter_end:
        txtime = total_dose/num_doses/txdose
        period = (inter_end-inter_start)/num_doses
        if (t-inter_start)%period < txtime:
            return txdose
        else:
            return 0
    else:
        return 0
    

def fcrand(t):
    if t>=chemo_start and t<=chemo_end:
        txtime = 160/num_doses/1
        period = (chemo_end-chemo_start)/num_doses
        if (t-chemo_start)%period < txtime:
            return 1
        else:
            return 0
    else:
        return 0
    
def ftrand(t):
    if t>=inter_start and t<=inter_end:
        txtime = 160/num_doses/1
        period = (inter_end-inter_start)/num_doses
        if (t-inter_start)%period < txtime:
            return 1
        else:
            return 0
    else:
        return 0

def evoLV(X, t):    
    
    C = X[0]
    v = X[1]
    v2 = X[2]
    
    if t>1300 and t<1400: #dummy again
        v2+=1
        
    def G(v):
        return grow(v)*(1-C/K)-death(v)-fc(t)*chemo(v)-ft(t)*target(v)
    
    dCdt = C*G(v)
        
    dvdt = k*((betar*r0/v0*(v/v0)**(betar-1))*(1-C/K)-grow(v)*C/K*(-betaa*(amax-1)/amax)-(betadelta*death(v))-fc(t)*(-betazeta*chemo(v))-ft(t)*(betazeta2*zeta02/v0*(v/v0)**(betazeta2-1)))
    
    #dvranddt = .5*((betar*r0/v0*(v/v0)**(betar-1))*(1-C/K)-grow(v)*C/K*(-betaa*(amax-1)/amax)-(betadelta*death(v))-fcrand(t)*(-betazeta*chemo(v))-ftrand(t)*(betazeta2*zeta02/v0*(v/v0)**(betazeta2-1))) #dummy function - txdose
    
    dvranddt = .1*((betar*r0/v0*(v2/v0)**(betar-1))*(1-C/K)-grow(v2)*C/K*(-1*(amax-1)/amax)-(betadelta*death(v2))-fc(t)*(-betazeta*chemo(v2))-ft(t)*(betazeta2*zeta02/v0*(v2/v0)**(betazeta2-1))) #dummy function - betaa
    
    #dvranddt = .3*((betar*r0/v0*(v2/v0)**(betar-1))*(1-C/K)-grow(v2)*C/K*(-1*(amax-1)/amax)-(betadelta*death(v2))-fc(t)*(-betazeta*chemo(v2))-ft(t)*(betazeta2*zeta02/v0*(v2/v0)**(betazeta2-1))) #dummy function - betar. For k just change evolv to .03
    
    #dvranddt = 1*((betar*r0/v0*(v2/v0)**(betar-1))*(1-C/K)-grow(v2)*C/K*(-1*(amax-1)/amax)-(betadelta*death(v2))-fc(t)*(-betazeta*chemo(v2))-ft(t)*(betazeta2*zeta02/v0*(v2/v0)**(betazeta2-1))) #dummy function - zeta0,2

        
    dxvdt = np.array([dCdt, dvdt,dvranddt])

    return dxvdt

# C = 0.2
# v = 1.3

# print(grow(v)*(1-C/K)-death(v)-target(v))
# print((betar*r0/v0*(v/v0)**(betar-1))*(1-C/K)-grow(v)*C/K*(-betaa*(amax-1)/amax)-(betadelta*death(v))-1*(betazeta2*zeta02/v0*(v/v0)**(betazeta2-1)))

# print(grow(v)*(1-C/K)-death(v)-chemo(v))
# print((betar*r0/v0*(v/v0)**(betar-1))*(1-C/K)-grow(v)*C/K*(-betaa*(amax-1)/amax)-(betadelta*death(v))-(-betazeta*chemo(v)))

intxv = np.array(IC)
time_sp = np.linspace(0,time,num=time*100-1)
pop = odeint(evoLV, intxv,time_sp,atol=1e-12,rtol=1e-12)
plt.figure()
plt.subplot(211)
plt.title('Sequential Therapy: High Tx Intensity, Low Evolv')
plt.plot(time_sp,pop[:,0],lw=3,color='b')
plt.xlim(0,inter_end)
plt.xticks([])
ax = plt.gca()
ax.axvspan(chemo_start, chemo_end, facecolor=(0.5,0.5,0.5))
ax.axvspan(inter_start, inter_end, facecolor='r') #firebrick
plt.ylabel('Pop Size, x')
plt.subplot(212)
ax = plt.gca()
ax.axvspan(chemo_start, chemo_end, facecolor=(0.5,0.5,0.5))
ax.axvspan(inter_start, inter_end, facecolor='r')
plt.plot(time_sp,pop[:,1],lw=3,color='b')
plt.ylabel('Biomarker Level, v')
plt.xlabel('Weeks')
plt.xticks([0,1344,4032,6720],[0,8,24,40])
plt.xlim(0,inter_end)
plt.show()

print(min(pop[:,0]))

# for i in np.arange(10,300,50):
#     total_dose = i
#     pop = odeint(evoLV, intxv,time_sp,atol = 1e-12, rtol=1e-12)
#     a = min(pop[:,0])
#     if a < 0:
#         a = 0
#     min_summ.append([i,a])

# min_summ = np.array(min_summ)
# plt.figure()
# #plt.xticks(np.arange(0,20, 2))
# plt.title('Effect of Dose per Tx on Sequential Therapy Protocol')
# plt.ylabel('Nadir of Cancer Population')
# plt.xlabel('Dose per Tx')
# plt.plot(min_summ[:,0],min_summ[:,1],lw=3,color='k')