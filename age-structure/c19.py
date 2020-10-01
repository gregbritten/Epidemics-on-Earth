from numpy import *

N0=1e4

S=N0-5
R=0
D=0
c=zeros(30)
c[0]=5
mu=0.01
beta0=0.2
beta=beta0

quar=90 # quarantine period
socialdist=0.05 # factor reducing beta during quarantine

tmax=1000

# c[0:7] exposed c[7:14] infectious, not ill c[14:30] ill, recovering/ deceased
infectious=hstack((zeros(6),[0.2,0.8],ones(30-8)))
recov=hstack((zeros(14),0.1*ones(14),[0,0]))
dec=hstack((zeros(14),0.02*ones(7),0.1*ones(7),[1,1]))
# for stochastic model, these would be probabilities, recalculated
# each day with the values above being the expected values [or using better
# numbers]

t=0
S_r=[];
R_r=[];
EI_r=[];
I_r=[]
D_r=[]
while t<=tmax:
    N=sum(c)+S+R
    I=sum(c*infectious)

    S_r.append(S)
    R_r.append(R)
    EI_r.append(sum(c))
    I_r.append(I)
    D_r.append(D)
    
    newc=beta*S/N*I   # number of new exposed  [make stochastic]
    R += sum(c*recov)-mu*R
    D += sum(c*dec) + c[29]
    c += -c*recov-c*dec-mu*c  # apply recovery, covid death, natural mortality
    c[1:30] = c[0:29] # number in day 0 _. day 1, day 1-> day 2, etc
    c[0]=newc
    S += -newc+mu*(N0-N)
    t += 1
    # this is the zsimple soicial distancing
    if t>=100 and t<100+quar:
        beta=socialdist*beta0
    else:
        beta=beta0
        
from matplotlib import pyplot as plt
t=arange(0,tmax+1)
#plt.plot(t,EI_r,t,S_r,t,R_r,t,D_r)
#plt.legend(["EI","S","R","D"])
plt.plot(t,EI_r)
plt.title("Exposed+infectious")
#plt.semilogy(t,EI_r)
#plt.show(block=False)
plt.show()
