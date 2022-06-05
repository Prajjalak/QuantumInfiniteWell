import numpy as npy
import matplotlib.pyplot as mplt

v0=1
x0=1
gamma=0.25

def V(x):        # Put your potential here
    if x<x0:
        return 0.0
    else:
        return v0+gamma*((x-x0)**2)

x = npy.linspace(0,10,1000)
h = x[1]-x[0]

#kinetic energy
Tee = npy.zeros((1000-2)**2).reshape(1000-2,1000-2) #kinetic energy matrix
for i in range(1000-2):
    for j in range(1000-2):
        if i==j:
            Tee[i,j]= -2
        elif npy.abs(i-j)==1:
            Tee[i,j]=1
        else:
            Tee[i,j]=0
            
# potential energy
Vee = npy.zeros((1000-2)**2).reshape(1000-2,1000-2) #potential energy matrix

for i in range(1000-2):
    for j in range(1000-2):
        if i==j:
            Vee[i,j]= V(x[i+1])
        else:
            Vee[i,j]=0

# hamiltonian
H = -Tee/(h**2) +Vee
val,vec=npy.linalg.eig(H) #eigenvalues

z = npy.argsort(val)
z = z[0:4] #first 4 elements
engy=val[z]
print('Energy eigenvalues:')
for j in range(4):
   print(engy[j],'in hbar^2/2m_e units')

for i in range(len(z)):
    y = []
    y=npy.append(y,vec[:,z[i]])
    y=npy.append(y,0)
    y=npy.insert(y,0,0)
    mplt.plot(x,y, label=i)
    
mplt.xlabel('x (Amstrong)')
mplt.ylabel(' Wavefunction $\psi(x)$')
mplt.legend()
mplt.title('Normalized Wavefn-s for an Anharmonic Oscillator')
mplt.grid(True)
mplt.axis([0,10,-0.08,0.08])
mplt.show()
