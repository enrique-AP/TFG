#Modelado hexaedro tridimensional
#Importamos modulos necesarios
#para solucionar o problema
import numpy as np
from numpy import sin, cos, exp
import matplotlib.pyplot as plt
from numpy import linalg as LA
import sympy
import math
from sympy.matrices import *
from sympy.plotting import plot
#############################
#############################
#Constantes fisicas do hexaedro
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

# Vector de momentos xerados polas forzas
M=([[ 0.5*a*(-k1-k2+k3+k4),       #Mx
      -0.25*a*a*(k1+k2+k3+k4),
      0.25*a*b*(-k1+k2-k3+k4),
      0.5*a*(-c1-c2+c3+c4),
      -0.25*a*a*(c1+c2+c3+c4),
      0.25*a*b*(-c1+c2-c3+c4)],
    
    [0.5*b*(-k1+k2+k3-k4),       #My
     0.25*a*b*(-k1+k2-k3+k4),
     -0.25*b*b*(k1+k2+k3+k4),
     0.5*b*(-c1+c2+c3-c4),
     0.25*a*b*(-c1+c2-c3+c4),
     -0.25*b*b*(c1+c2+c3+c4)],
    
    [0,   #Mz
     0,
     0,
     0,
     0,
     0]])


I_inv = np.linalg.inv(I) #Inversa do tensor de inercia


w_dot2=I_inv @ M
#Quedanos unha matriz cunha fila de ceros
#Esta fila representa a derivada segunda do xiro
#arredor do eixo z
#Eliminaremola da  nosa matriz ya que dito angulo


# Aceleracions no eixo z
z_dot2=([[(-k1-k2-k3-k4)/m,
          0.5*a*(-k1-k2+k3+k4)/m,
          0.5*b*(-k1+k2+k3-k4)/m,
          (-c1-c2-c3-c4)/m,
          0.5*a*(-c1-c2+c3+c4)/m,
          0.5*b*(-c1+c2+c3-c4)/m]])


#Matriz A
#Parte derivadas primeiras
A=([[0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 1]])
A=np.concatenate((A, z_dot2,w_dot2), axis=0)
#Engadimoslle as filas correspondentes a
#z_dot2 phi_dot2 e ksi_dot2

A     = np.delete(A, (6), axis=0)
#Elimino a ultima fila de 0, xa que a rotacion ao redor
#do eixo z permanece constante

#Declaramos a variable simbolica t para obter unha
#funcion analitica en funcion do tempo
t = sympy.symbols('t')


#Calculamos os autovalores e autovectores da matriz A
w, v  = LA.eig(A)


#Condicions iniciais do problema
z_0=0  
phi_0=0 
ksi_0=3*math.pi/180
z_dot_0=0 
phi_dot_0=0 
ksi_dot_0=0

# Cos valores iniciais, substituimos t=0 na expresion 
#da nosa solucion para obter as constantes
Yin = ([z_0,phi_0,ksi_0,z_dot_0,phi_dot_0,ksi_dot_0])
C_inv = np.linalg.inv(v)
c=C_inv @ Yin   #Vector de constantes C1,C2,C3,C4,C5,C6


#Produto de cada constante polo seu autovector asociado

#Cada fila da matriz v correspondese a
#unha variable que queremos calcular
#Multiplicamos cada autovector pola sua
#constante c e por exp(lambda*t)
#Autovector 0
c0v0=c[0]*v[:,0]*sympy.exp(w[0]*t)  #w=lambda
#Autovector 1
c1v1=c[1]*v[:,1]*sympy.exp(w[1]*t)
#Autovector 2
c2v2=c[2]*v[:,2]*sympy.exp(w[2]*t)
#Autovector 3
c3v3=c[3]*v[:,3]*sympy.exp(w[3]*t)
#Autovector 4
c4v4=c[4]*v[:,4]*sympy.exp(w[4]*t)
#Autovector 5
c5v5=c[5]*v[:,5]*sympy.exp(w[5]*t)


#O vector de solucions x(t) e igual a combinacion
#lineal dos autovectores multiplicados por exp(lambda*t)
#e por unhas constantes c dadas polas condicions iniciais
x_t=c0v0+c1v1+c2v2+c3v3+c4v4+c5v5

z=x_t[0]
z_dot=x_t[3]
phi=x_t[1]
phi_dot=x_t[4]
ksi=x_t[2]
ksi_dot=x_t[5]

#Amosamos por pantalla as expresions
#temporais das nosas variables
print("z(t)=",z)
print("z_dot(t)=",z_dot)
print("phi(t)=",phi)
print("phi_dot(t)=",phi_dot)
print("ksi(t)=",ksi)
print("ksi_dot(t)=",ksi_dot)