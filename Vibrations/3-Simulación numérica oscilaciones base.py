#########################################
#########################################
#Simulacion numerica caixa tridimensional

import matplotlib.pyplot as plt
import numpy as np
import math
from numpy import sin, cos
from numpy import linalg as LA
import scipy.integrate as integrate
from scipy.integrate import solve_ivp
import time


#############################
#############################
#Constantes fisicas da caixa

#Masa
m=1
#Inercia
ixx=1/60
ixy=0
ixz=0
iyy=1/60
iyz=0
izz=1/60
I=np.array([[ixx, ixy, ixz],
            [ixy, iyy, iyz],
            [ixz, iyz, izz]])
#Distancias entre os amortiguadores
a=0.2
b=0.2
##################################
##################################
#Constantes amortiguadores

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
M=([[ 0.5*a*(-k1-k2+k3+k4), -0.25*a*a*(k1+k2+k3+k4),
      0.25*a*b*(-k1+k2-k3+k4),      0.5*a*(-c1-c2+c3+c4),
      -0.25*a*a*(c1+c2+c3+c4),           0.25*a*b*(-c1+c2-c3+c4)],
    [0.5*b*(-k1+k2+k3-k4),   0.25*a*b*(-k1+k2-k3+k4),
     -0.25*b*b*(k1+k2+k3+k4),       0.5*b*(-c1+c2+c3-c4),
     0.25*a*b*(-c1+c2-c3+c4),        -0.25*b*b*(c1+c2+c3+c4)],
    [0,                      0,                            0,
     0,                         0,                               0]])


I_inv = np.linalg.inv(I)


w_dot2=I_inv @ M 
#Eliminaremos a rotacion ao redor do eixo z xa que e despreciable
# Oscilacions no eixo z
z_dot2=([[(-k1-k2-k3-k4)/m,0.5*a*(-k1-k2+k3+k4)/m, 0.5*b*(-k1+k2+k3-k4)/m,
          (-c1-c2-c3-c4)/m, 0.5*a*(-c1-c2+c3+c4)/m, 0.5*b*(-c1+c2+c3-c4)/m]])


#Matriz A
#Parte derivadas primeiras
A=([[0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 1]])
A=np.concatenate((A, z_dot2,w_dot2), axis=0)

A = np.delete(A, (6), axis=0)
#Elimino a ultima fila xa que non consideramos a
#rotacion arredor do eixo z


################################
#Efecto da vibracion da base
#Tolerancia metodo de integracion numerico
tol=2.3/10**13
#Array temporal
nstep = 1000000
t0 = 0.
tf = 1
dt = (tf-t0)/nstep
t = np.linspace(t0,tf,nstep+1)
#Amplitudes e constantes
A1=0.001
A2=0.001
A3=0.001
A4=0.001
print(-k1*A1-k2*A2+k3*A3+k4*A4)
frecuencia=1000
w=2*math.pi * frecuencia

#Condicions iniciais
#Suponhemos que o sistema esta en reposo
#no instante inicial
z_0=0
phi_0=0
ksi_0=0
z_dot_0=0
phi_dot_0=0
ksi_dot_0=0
Yin = ([z_0,phi_0,ksi_0,z_dot_0,phi_dot_0,ksi_dot_0])

t1 = time.time()       

def F(t, Y):
    FF = np.zeros_like(Y)       
    #Forzas e momentos xerados pola base
    F=(-k1*A1-k2*A2-k3*A3-k4*A4)*cos(w*t)+(c1*A1+c2*A2+c3*A3+c4*A4)*w*sin(w*t)
    Mx=0.5*a*(-k1*A1-k2*A2+k3*A3+k4*A4)*cos(w*t)-0.5*a*(-c1*A1-c2*A2+c3*A3+c4*A4)*w*sin(w*t)
    My=0.5*b*(-k1*A1+k2*A2+k3*A3-k4*A4)*cos(w*t)-0.5*b*(-c1*A1+c2*A2+c3*A3-c4*A4)*w*sin(w*t)
    Mz=0
    M_suelo = np.array([[Mx],[My],[Mz]])
    ac_angular=I_inv @ M_suelo

    #Aceleracions angulares
    FF[0] = Y[3]
    FF[1] = Y[4]
    FF[2] = Y[5]
    FF[3] =F/m +Y[0]*A[3][0]+Y[1]*A[3][1]+Y[2]*A[3][2]+Y[3]*A[3][3]+Y[4]*A[3][4]+Y[5]*A[3][5]
    FF[4]= ac_angular[0]+Y[0]*A[4][0]+Y[1]*A[4][1]+Y[2]*A[4][2]+Y[3]*A[4][3]+Y[4]*A[4][4]+Y[5]*A[4][5]
    FF[5]= ac_angular[1]+Y[0]*A[5][0]+Y[1]*A[5][1]+Y[2]*A[5][2]+Y[3]*A[5][3]+Y[4]*A[5][4]+Y[5]*A[5][5]

    return FF

#Solucionamos numericamente
sol = solve_ivp(F, (t0,tf), Yin, t_eval=t, method='RK45', atol=tol, rtol=tol)

###############################################################################
###############################################################################
#Graficas
Y = sol.y

plt.figure()
plt.plot(t, Y[0, :],color='blue')# Bode magnitude plot
plt.xlabel(r'Tempo (s)')
plt.ylabel(r'Desprazamento (m)')
plt.title("Variacion de z debido a oscilacion da base")
plt.grid(True)


plt.figure()
plt.plot(t,Y[1, :]*180/math.pi,color='blue')# Bode phase plot
plt.title(r"Variacion de $\varphi$ debido a oscilacion da base")
plt.xlabel(r'Tempo (s)')
plt.ylabel(r'Angulo(graos)')
plt.grid(True)


plt.figure()
plt.plot(t,Y[2, :]*180/math.pi,color='blue')
plt.title(r"Variacion de $\Psi$ debido a oscilacion da base")
plt.xlabel(r'Tempo (s)')
plt.ylabel(r'Angulo (graos)')
plt.grid(True)

t2 = time.time()
print('Tempo de execucion. =',t2-t1)
#Mide o tempo transcorrido na execucion do programa

plt.show()

