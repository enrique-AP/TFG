clear all
syms z1
syms z1dot
syms z2
syms z2dot
syms theta1
syms theta1dot
syms theta2
syms theta2dot
syms phi1
syms phi1dot
syms phi2
syms phi2dot
syms V
syms d

k1=3*2361;
c1=1.4843*3;
k2=3*1988;
c2=2.0066*3;
x=[0.0725 0.0725 -0.0725 -0.0725];
y=[0.13555 -0.13555 0.13555 -0.13555];
x2=[0.04078 0.04078 -0.04078 -0.04078];
y2=[0.04398 -0.04398 0.04398 -0.04398];
V=0;
d=0;
m1=0.44816;
m2=0.72456;
I1=[3 0.2;0.1 1];
I2=[3 0.1;0.2 2];
I1inv=I1^-1;
for i=1:4
    V=V+0.5*3*k1*(z1-x(i)*theta1-y(i)*phi1)^2++0.5*3*k2*(z2-z1+x2(i)*(theta1-theta2)+y2(i)*(phi1-phi2))^2;
    d=V+0.5*3*c1*(z1dot-x(i)*theta1dot-y(i)*phi1dot)^2++0.5*3*c2*(z2dot-z1dot+x2(i)*(theta1dot-theta2dot)+y2(i)*(phi1dot-phi2dot))^2;
end
eq1=-(diff(V,z1)+diff(d,z1dot));
eq2=-(diff(V,z2)+diff(d,z2dot));
eq3=-(diff(V,theta1)+diff(d,theta1dot));
eq4=-(diff(V,theta2)+diff(d,theta2dot));
eq5=-(diff(V,phi1)+diff(d,phi1dot));
eq6=-(diff(V,phi2)+diff(d,phi2dot));

eq1=eq1/m1;
eq2=eq2/m2;
eqvec1=I1^-1*[eq3; eq5];
eqvec2=I2^-1*[eq4; eq6];
I1inv=I1^-1;

eqns = [z1dot==0,
        z2dot==0,
        theta1dot==0,
        theta2dot==0,
        phi1dot==0,
        phi2dot==0,
        eq1 == 0,
        eq2  == 0,
        eqvec1(1)   == 0,
        eqvec2(1)   == 0,
        eqvec1(2)   == 0,
        eqvec2(2)   == 0;];
vars = [z1 z2 theta1 theta2 phi1 phi2 z1dot z2dot theta1dot theta2dot phi1dot phi2dot];
A = equationsToMatrix(eqns,vars);

B=[0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    1 0 0 0 0 0;
    
    0 0 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 0;];
C=[0 1 0 0 0 0 0 0 0 0 0 0];
C=eye(12);
D=0;
A=double(A);
sys = ss(A,B,C,D);
w = linspace(0,200*2*pi,500);
[mag,ph,wout,sdmag,sdphase]  = bode(sys,w);

z1vec=zeros(500,1);
z2vec=zeros(500,1);
theta1vec=zeros(500,1);
theta2vec=zeros(500,1);
phi1vec=zeros(500,1);
phi2vec=zeros(500,1);
for i=1:length(z1vec)
    z1vec(i,1)=mag(2,1,i);
    z2vec(i,1)=mag(2,3,i);
    theta1vec(i,1)=mag(3,3,i);
    theta2vec(i,1)=mag(4,3,i);
    phi1vec(i,1)=mag(5,5,i);
    phi2vec(i,1)=mag(6,5,i);
    
end


pixelz2vec=z1vec*50*6000/(1500*23.5);

 
figure
plot(w/(2*pi),z2vec)
