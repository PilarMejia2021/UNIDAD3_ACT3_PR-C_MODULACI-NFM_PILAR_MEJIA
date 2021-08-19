
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import pylab as pl
from math import  pi
from math import log10
plt.close('all')

Vm= float(input('Ingrese el valor de Amplitud de la Moduladora --> '))
Fm= int(input('Ingrese el valor de Frecuencia de la Moduladora --> '))
Vc= float(input('Ingrese el valor de Amplitud de la Portadora --> '))
Fc= int(input('Ingrese el valor de Frecuencia de la Portadora --> '))
kf= float(input('Ingrese el factor de sensibilidad de frecuencia --> '))
n= float(input('Ingrese el numero de periodos --> '))
print()

#variables que usaremos a lo largo del programa
z= 50
∆f=kf*Vm
β=∆f/Fm

Fs=50000
x=0
n0=[]
bessel=[]
f=np.arange(0,10,1)

#ecuaciones para hallar bessel
for i in range (0,len(f)):
    x = round (sp.jv(i,β),2)
    bessel.append(x)

n_positivos=bessel[1:11];
n_negativos=np.flip(n_positivos);
n0.append(bessel[0]);

jn=np.concatenate((n_negativos,n0,n_positivos))

nB=4
BWb=2*Fm*nB
BWc=2*(∆f*Vm)

#valores para las frecuencias

f_ns=[]
f_ps=[]
F0=[]
F0.append(Fc)

for f_inicial in range(0,len(f)):

    if f_inicial==0:
        f_1=Fc-Fm;
        f_inicial=f_1
else:
        f_1=f_1-Fm;
        f_inicial = f_1 ;

f_ns.append(f_inicial);

finv_ns=np.flip(f_ns);

for f_inicial in range (0,len(f)):

    if f_final==0:
        f_1=Fc+Fm;
        f_final = f_1;
    else:
        f_1=f_1+Fm;
        f_final =f_1;

    f_ps.append(f_final);
finv_ps=np.flip(f_ps);
Fn = np.concatenate((finv_ns,F0,f_ps))

t=np.arange(0,n*1/Fm,1/Fs)

#hallar Vc * Jn
f_VcJn=[]
VcJn = 0
VcJn = np.round (abs(jn*vc)/(np.sqrt(2)),2)
f_VcJn.append(VcJn)

#hallar valores en dB de Jn*Vc
f_VndB=[]
VndB=0
VndB=np.round(abs(20*np.log10(VcJn)),2)
f_VndB.append(VndB)

#hallar  potencia en watts (w)
f_PnW=[]
PnW=0
PnW=abs(((jn*Vc)**2)/100)
f_PnW.append(PnW)

#hallar  potencia en dBm
f_PndBm=[]
PndBm=0
PndBm=np.round(abs((10*np.log10(PnW*1000))),2)
f_PndBm.append(PndBm)

#calculo de ecuaciones presentes en el programa de modulación
Vportadora=Vc*np.cos(2*pi*Fc*t);
Vmoduladora=Vm*np.sin(2*pi*Fm*t);
Vfm=Vc*np.cos(2*pi*Fc*t+β*np.sin(2*pi*Fm*t));

#formulas resultantes
print('RESULTADOS MODULACION FM')
print()
print('{:^10} {:^10} {:^10} {:^10}'.format('∆ƒ','β','BWb','BWc'))
print('{:^10} {:^10} {:^10} {:^10}'.format(∆ƒ,β,BWb,BWc))
print()
print('{:^10} {:^9} {:^9} {:^9} {:^9}'.format('Jn','Fn','Vc*Jn','Vn(dB)','Vn(dBm)', 'Vn(dBm)'))
for formatted in map ('{:^10} {:^10} {:^9} {:^10} {:^10}'.format,jn,Fn,VcJn,VndB,PndBm):
    print (formatted);
print();
print("LA ECUACIÓN PORTADORA ES:");
print("Vc(t)=",Vc,"cos(2π",Fc,"t)");
print("LA ECUACIÓN MODULADORA ES:");
print("Vm(t)=",Vm,"sen(2π",Fm,"t)");
print();
print("LA ECUACIÓN GENERAL PARA FM ES:");
print("Vfm(t)=",Vc,"cos[2π",Fc,"t + ",β,"sen(2π",Fm,"t)]");
print()

#graficas resultantes
fig=plt.figure();
fig,plt.subplot (1,1,1);
plt.plot(t,Vportadora,color="blue",linewidth=0.8);
plt.title('Señal Portadora');
plt.xlabel('tiempo');
plt.ylabel ('amplitud');
plt.grid(True);

fig1=plt.figure();
fig1,plt.subplot (1,1,1);
plt.plot(t,Vmoduladora,color="black",linewidth=0.8);
plt.title('Señal Moduladora');
plt.xlabel('tiempo');
plt.ylabel ('amplitud');
plt.grid(True);

fig2=plt.figure();
fig2,plt.subplot (1,1,1);
plt.plot(t,Vfm,color="red",linewidth=0.8);
plt.title('Modulacion de Frecuencia');
plt.xlabel('tiempo');
plt.ylabel ('amplitud');
plt.grid(True)












