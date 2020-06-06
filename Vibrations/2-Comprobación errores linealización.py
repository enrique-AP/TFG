#Modulos necesarios
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos, exp
import scipy.integrate as integrate
from scipy.integrate import solve_ivp
import time
import math


#############################
#############################
#Masa do hexaedro
m   =  1
#Tensor de inercia do hexaedro
ixx = 1/60
ixy = 0
ixz = 0
iyy = 1/60
iyz = 0
izz = 1/60
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
c1=10
#Muelle 2
k2=4000
c2=10
#Muelle 3
k3=4000
c3=10
#Muelle 4
k4=4000
c4=10
###################################
###################################
#Calculo matrices
I_inv = np.linalg.inv(I) 
#Condicions iniciais do problema
z_0=0  
phi_0=0
ksi_0=0.1
z_dot_0=0 
phi_dot_0=0 
ksi_dot_0=0

# Cos valores iniciais, substituimos t=0 na expresion 
#da nosa solucion para obter as constantes
Yin = ([z_0,phi_0,ksi_0,z_dot_0,phi_dot_0,ksi_dot_0])

#Array temporal,
#representa os puntos de discretizacion temporal
nstep = 1000
t0 = 0.
tf = 0.5 
dt = (tf-t0)/nstep
t = np.linspace(t0,tf,nstep+1)
#Funcion derivada da solucion
#Actualizase co tempo e coa evolucion das variables
def F(t, Y):
     # Crea un array de ceros do tamanho da solucion Y
    FF = np.zeros_like(Y)       
    F1=(-k1*Y[0]-k1*0.5*a*sin(Y[1])-k1*0.5*b*sin(Y[2])
        -c1*Y[3]-c1*0.5*a*cos(Y[1])*Y[4]-c1*0.5*b*cos(Y[2])*Y[5])
    
    F2=(-k2*Y[0]-k2*0.5*a*sin(Y[1])+k2*0.5*b*sin(Y[2])
        -c2*Y[3]-c2*0.5*a*cos(Y[1])*Y[4]+c2*0.5*b*cos(Y[2])*Y[5])
    
    F3=(-k3*Y[0]+k3*0.5*a*sin(Y[1])+k3*0.5*b*sin(Y[2])
        -c3*Y[3]+c3*0.5*a*cos(Y[1])*Y[4]+c3*0.5*b*cos(Y[2])*Y[5])
    
    F4=(-k4*Y[0]+k4*0.5*a*sin(Y[1])-k4*0.5*b*sin(Y[2])
        -c4*Y[3]+c4*0.5*a*cos(Y[1])*Y[4]-c4*0.5*b*cos(Y[2])*Y[5])
    
    
    Mx=0.5*a*cos(Y[1])*(F1+F2-F3-F4)
    My=0.5*b*cos(Y[2])*(+F1-F2-F3+F4)
    Mz=0
    #Vector de derivadas
    FF[0] = Y[3]
    FF[1] = Y[4]
    FF[2] = Y[5]
    FF[3] = F1+F2+F3+F4 #Aceleracion eixo z
    #Aceleracions angulares
    #Multiplicamos a inversa do tensor de inercia
    #Polos momentos
    FF[4]=I_inv[0][0]*Mx+I_inv[0][1]*My
    FF[5]= I_inv[1][0]*Mx+I_inv[1][1]*My

    return FF
################################################################
################################################################
#Solucionamos numericamente
#Tolerancia do metodo numerico
tol=2.3/10**10
sol = solve_ivp(F, (t0,tf), Yin, t_eval=t, method='RK45',
                atol=tol, rtol=tol)

Y = sol.y #Garda os valores da solucion nunha matriz


#Aumentamos a precision do metodo de integracion
tol=2.3/10**14
sol2 = solve_ivp(F, (t0,tf), Yin, t_eval=t, method='RK45',
                 atol=tol, rtol=tol)
Y2 = sol2.y
################################################################
################################################################
#Graficas

#Erros numericos
plt.figure()
plt.plot(t,abs(Y2[2, :]-Y[2, :])*180/math.pi,color='blue')
plt.xlabel(r'Tempo (s)')
plt.ylabel(r'Xiro (graos)')
plt.title(r"Erros numericos", fontsize=16, color='gray')
plt.legend()

#Solucion analitica do problema
ksi=(0.1)*exp(t*(-12.0))*cos(-97.2419662491458*t)

#Diferenza entre os dous metodos
plt.figure()
plt.plot(t,abs(Y2[2, :]-ksi)*180/math.pi,color='blue')
plt.xlabel(r'Tempo (s)')
plt.ylabel(r'Xiro (graos)')
plt.title(r"Erros debido a linealizacion do problema ",
          fontsize=16, color='gray')
plt.legend()

#Representamos ambas solucions para ver o desfase
plt.figure()
plt.plot(t,ksi*180/math.pi,color='blue', label="Linealizado")
plt.plot(t,Y2[2,:]*180/math.pi,color='red', label="Non linealizado")
plt.xlabel(r'Tempo (s)')
plt.ylabel(r'Xiro (graos)')
plt.title(r"Erros debido a linealizacion do problema ",
          fontsize=16, color='gray')
plt.legend()
 
plt.show()