#Modelado caixa tridimensional
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
import control as cnt
import math
import time


#############################
#############################
#Constantes fisicas da caixa

#Masa
m=1
#Inercia
ixx=1/45
ixy=0.0002
ixz=0.00004
iyy=1/60
iyz=0.0005
izz=1/60
I=np.array([[ixx, ixy, ixz],
            [ixy, iyy, iyz],
            [ixz, iyz, izz]])
#Distancias entre os amortiguadores
a=0.2
b=0.2
##################################
##################################
#Constantes muelles e amortiguadores

# k=muelle c=amortiguador (damping)
#Muelle 1
k1=4000
c1=5
#Muelle 2
k2=4000
c2=5
#Muelle 3
k3=4000
c3=5
#Muelle 4
k4=4000
c4=5
###################################
###################################
#Calculo matrices

#Parte de inercia
M=([[ 0.5*a*(-k1-k2+k3+k4), -0.25*a*a*(k1+k2+k3+k4), 0.25*a*b*(-k1+k2-k3+k4),
      0.5*a*(-c1-c2+c3+c4), -0.25*a*a*(c1+c2+c3+c4), 0.25*a*b*(-c1+c2-c3+c4)],
    
    [0.5*b*(-k1+k2+k3-k4),   0.25*a*b*(-k1+k2-k3+k4),-0.25*b*b*(k1+k2+k3+k4),
     0.5*b*(-c1+c2+c3-c4),   0.25*a*b*(-c1+c2-c3+c4),-0.25*b*b*(c1+c2+c3+c4)],
    
    [0,                      0,                       0,
     0,                      0,                       0]])


I_inv = np.linalg.inv(I)
print(I_inv)

w_dot2=I_inv @ M  #Quedanos unha matriz coa ultima fila de 0

#Esta fila representa a derivada
#segunda do xiro arredor do eixo z
#Eliminaremos esta fila xa que
#non non se produce ningun xiro arredor do eixo z
# Oscilacions no eixo z
z_dot2=([[(-k1-k2-k3-k4)/m,
          0.5*a*(-k1-k2+k3+k4)/m,
          0.5*b*(-k1+k2+k3-k4)/m,
          (-c1-c2-c3-c4)/m,
          0.5*a*(-c1-c2+c3+c4)/m,
          0.5*b*(-c1+c2+c3-c4)/m]])


#Matriz A
#Parte derivadas primeras
A=([[0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 1]])
A=np.concatenate((A, z_dot2,w_dot2), axis=0)

A = np.delete(A, (6), axis=0)
#Elimino a ultima fila de 0, xa que a
#rotacion arredor do eixo z non nos interesa


B = np.array([[0,0,0,0,0,0],
              [0,0,0,0,0,0],
              [0,0,0,0,0,0],
              [-(k1+k2+k3+k4)/m,0,0,-(c1+c2+c3+c4)/m,0,0]])
##Engadimos unha columna de o para que as dimensions sexan compatibles
B=np.concatenate((B,w_dot2), axis=0)
B = np.delete(B, (6), axis=0)
print(w_dot2)

#Matriz C
C = np.array([[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0],
              [0, 0, 1, 0, 0, 0],[0, 0, 0, 1, 0, 0],
              [0, 0, 0, 0, 1, 0],[0, 0, 0, 0, 0, 1]])

#Matriz D
D = np.array([[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],[0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],[0, 0, 0, 0, 0, 0]])


################################################
################################################
#Calculo das variables e funcions de transferencia

#Sistema e funcions de transferencia

sys         = signal.lti(A, B, C, D)
num1, den1  = signal.ss2tf(A, B, C, D, 0)
num2, den2  = signal.ss2tf(A, B, C, D, 1)
num3, den3  = signal.ss2tf(A, B, C, D, 2)

transf_11   = cnt.tf(num1[0,:], den1)
transf_12   = cnt.tf(num1[1,:], den1)
transf_13   = cnt.tf(num1[2,:], den1)

transf_21   = cnt.tf(num2[0,:], den2)
transf_22   = cnt.tf(num2[1,:], den2)
transf_23   = cnt.tf(num2[2,:], den2)

transf_31   = cnt.tf(num3[0,:], den3)
transf_32   = cnt.tf(num3[1,:], den3)
transf_33   = cnt.tf(num3[2,:], den3)

#Prints funcions de transferencia do sistema
print(transf_11)
print(transf_12)
print(transf_13)
print(transf_21)
print(transf_22)
print(transf_23)
print(transf_31)
print(transf_32)
print(transf_33)


#Calculo diagramas de bode
w, mag, phase = signal.bode([num1[0,:],den1])
w2, mag2, phase2 = signal.bode([num2[1,:],den2])
w3, mag2, phase3 = signal.bode([num3[2,:],den3])
#Graficas
#Resposta en z ao aplicar un estimulo en z
plt.figure()
plt.subplot(211)
plt.plot(w/(2*math.pi), mag,color='blue')
plt.xlabel(r'Frecuencia (Hz)')
plt.ylabel(r'Magnitude (dB)')
plt.title("Variacion da amplitude coa frecuencia")
plt.grid(True)
plt.subplot(212)
plt.plot(w/(2*math.pi), phase,color='blue')
plt.title("Variacion da fase coa frecuencia")
plt.xlabel(r'Frecuencia (Hz)')
plt.ylabel(r'Fase (graos)')
plt.grid(True)

#Resposta en phi ao aplicar un estimulo en phi
plt.figure()
plt.subplot(211)
plt.plot(w2/(2*math.pi), mag2,color='blue')
plt.xlabel(r'Frecuencia (Hz)')
plt.ylabel(r'Magnitude (dB)')
plt.title("Variacion da amplitude coa frecuencia")
plt.grid(True)
plt.subplot(212)
plt.plot(w2/(2*math.pi), phase2,color='blue')
plt.title("Variacion da fase coa frecuencia")
plt.xlabel(r'Frecuencia (Hz)')
plt.ylabel(r'Fase (graos)')
plt.grid(True)

#Resposta en ksi ao aplicar un estimulo en ksi
plt.figure()
plt.subplot(211)
plt.plot(w3/(2*math.pi), mag,color='blue')
plt.xlabel(r'Frecuencia (Hz)')
plt.ylabel(r'Magnitude (dB)')
plt.title("Variacion da amplitude coa frecuencia")
plt.grid(True)
plt.subplot(212)
plt.plot(w3/(2*math.pi), phase3,color='blue')
plt.title("Variacion da fase coa frecuencia")
plt.xlabel(r'Frecuencia (Hz)')
plt.ylabel(r'Fase (graos)')
plt.grid(True)
plt.show()