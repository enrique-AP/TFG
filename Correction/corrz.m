function y=corrz(Vertices)
Vertices=Vertices;

%Parámetros da aproximación polinómica
p00 =      0.0121;  
p10 =    -0.04407;  
p01 =    -0.01099;  
p20 =     0.05713;  
p11 =     0.02823;  
p02 =    0.008713;  
p30 =    -0.02325;  
p21 =    -0.04802;  
p12 =     0.01995;  
p03 =    -0.02042;  
p31 =     0.02181;  
p22 =   -0.006337;  
p13 =   -0.005811;  
p04 =     0.01136;  
p00a =   -0.001331;
p10a =   -0.000368;
p01a =   0.0006739;
p20a =  -0.0001234;
p11a =   0.0003032;
p02a =   0.0003053;
p30a =  -0.0003463;
p21a =  -0.0004219;
p12a =   3.103e-05;
p03a =  -0.0002208;
p40a =   0.0008775;
p31a =    0.000263;
p22a =   -8.15e-05;
p13a =  -0.0001974;
p04a =  -7.956e-05;
p50a =  -0.0003419;
p41a =  -7.643e-05;
p32a =   1.978e-05;
p23a =   0.0001153;
p14a =   1.659e-05;

for i=1:length(Vertices(:,1))
    x=(Vertices(i,4)+Vertices(i,5)+Vertices(i,6))/(3*255);
    y=abs(Vertices(i,3));
    delta = (p00a + p10a*x + p01a*y + p20a*x^2 + p11a*x*y + p02a*y^2 + p30a*x^3 + p21a*x^2*y+ p12a*x*y^2 + p03a*y^3 + p40a*x^4 + p31a*x^3*y + p22a*x^2*y^2 + p13a*x*y^3 + p04a*y^4 + p50a*x^5 + p41a*x^4*y + p32a*x^3*y^2 + p23a*x^2*y^3 + p14a*x*y^4)*4;
    Vertices(i,3)=Vertices(i,3)+delta;
    
            

 
end
writematrix(Vertices,'Correxido.txt','Delimiter','tab')


end