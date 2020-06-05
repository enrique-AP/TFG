import matplotlib.pyplot as plt
import numpy as np
import numpy as np
import math
from numpy import sin
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import minimize
import pandas as pd
from scipy.fftpack import fft, fftfreq
data = np.loadtxt('AVIONICA.txt',delimiter='\t',skiprows=1)
#Le o ficheiro no que se atopan os datos do acelerometro
###########################################################
dt = 0.001  #Paso temporal
#E igual a 0.001s xa que a frecuencia de mostreo do 
# acelerometro e 1000HZ
t = data[:,0]
n=len(t)
ac = data[:,1]
Y = fft(ac) / n
frq = fftfreq(n, dt)

###################################################
#Mosta por pantalla as aceleracions medidas
plt.figure()
plt.plot(t, ac,color='blue')
plt.title(r"Aceleracions medidas polo acelerometro ")
plt.xlabel(r'Tempo (s)')
plt.ylabel(r'Aceleracion ')
plt.grid(True)
#Grafica do espectro de frecuencias
plt.figure()
plt.plot(abs(frq), abs(Y),color='blue')
plt.title(r"Espectro de frecuencias ")
plt.xlabel(r'Frecuencia  (Hz)')
plt.ylabel(r' Magnitude ')
plt.grid(True)
plt.show()


#Medimos na grafica dous maximos consecutivos para determinar
#o amortiguamento
max1=2.16409
max2=0.474804
m=2.193 #masa do sistema
w=15.0754*2*math.pi #Frecuencia de resonancia medida na grafica
n=1
#Despexamos a ecuacion do decaemento exponencial para determinar
#o valor de gamma
gamma = (-2*math.pi+math.sqrt(4*math.pi*math.pi+
4*math.log(max1/max2)*math.log(max1/max2)))/(2*math.log(max1/max2))
print(gamma)

#A partir de gamma, calculamos as constantes c e k
k=m*w*w/(1-gamma**2)
c=2*k*gamma/(2*np.sqrt(k*m))
#############################
#Restamoslle o amortiguamento e rixidez de
#Cada muelle para obter as constantes non
#tidas en conta da estrutura
print(c-4*4.5)
print(k-3000*4)