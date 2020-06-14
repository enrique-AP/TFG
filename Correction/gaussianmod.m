function y=gaussianmod(Vertices,k,rmax,intensity)

%Parámetros do axuste da superficies
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


%Axustes do método de entrada
intensidade=intensity;
Vertices=Vertices;
k=k;
rmax=rmax;

%Tesselación da nube de puntos
[idx,C]= kmeans(Vertices(:,1:3),k);
Vertices=[Vertices idx];


%Cálculo dos vectores normais
normal=zeros(k,3);
for i=1:k
    vv2=Vertices(Vertices(:,7)==i,:);
    DM = [vv2(:,1:2) ones(size(vv2(:,1)))];                             % Design Matrix
    B = DM\vv2(:,3); 
    normal(i,:)=B;
end
%Triangulación e interpolación do vector normal
Nx = scatteredInterpolant(C(:,1), C(:,2), C(:,3), normal(:,1));
Ny = scatteredInterpolant(C(:,1), C(:,2), C(:,3), normal(:,2));
Nz = scatteredInterpolant(C(:,1), C(:,2), C(:,3), normal(:,3));
r=zeros(k,k);
vv4=[];
for i=1:k
    inter=i;
    vv2=Vertices(Vertices(:,7)==i,:);
    for j=1:k
        r=(C(i,1)-C(j,1))^2+(C(i,2)-C(j,2))^2+(C(i,3)-C(j,3))^2;
        vv2=Vertices(Vertices(:,7)==i,:);
        vv3=vv2;
        if r<rmax
            vv3=[vv3; Vertices(Vertices(:,7)==j,:)];
       
        end
        
    end
    
    for z=1:length(vv2)
        
            %Cálculo da intensidade do filtro
            x=(vv2(i,4)+vv2(i,5)+vv2(i,6))/(3*255);
            y=abs(Vertices(i,3));
            suma=0;
            ponder=0;
            sigma=p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y + p12*x*y^2 + p03*y^3 + p31*x^3*y + p22*x^2*y^2 + p13*x*y^3 + p04*y^4;
            sigma=sigma*intensidade;
            vec2=[0 0 0];
            ponder=0;
            
            %Cálculo do plano interpolado en cada punto
            normalvec=[Nx(vv2(z,1),vv2(z,2),vv2(z,3)) Ny(vv2(z,1),vv2(z,2),vv2(z,3)) Nz(vv2(z,1),vv2(z,2),vv2(z,3))];

        for w=1:length(vv3)
            diff=[(vv3(w,1)-vv2(z,1)) (vv3(w,2)-vv2(z,2)) (vv3(w,3)-vv2(z,3))] ;
            vec=normalvec*diff';
            ponder=ponder+exp(-diff*diff'/(2*sigma*sigma))/(sigma*sqrt(2*pi));
            vec2=vec2+normalvec*vec*exp(-diff*diff'/(2*sigma*sigma))/(sigma*sqrt(2*pi));
            
        end
        vec3=[vec2 0 0 0 0];
        vv2(z,:)=vv2(z,:)+vec3/ponder;
        
        
    end
    vv4=[vv4;vv2];
    
    
end

writematrix(vv4,'Correxido.txt','Delimiter','tab')


end