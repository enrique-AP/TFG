clear all
ISO=800;
d=1500;
dfar=2000;
f=50;
pixel=linspace(0,10,1000);
anchosensor=23.5;
resolucion=6000;
iluminancia=5378.1512
EV100=log(iluminancia/2.5)/log(2)
%%%%%%%%%%%%%%%%%%%%
list=[3.5 4 4.5 5 5.6 6.3 7.1 8 9 10 11 13 14 16 18 20 22];
fnumber=sqrt(375*(d*f/(d-f)-dfar*f/(dfar-f)));
for i=1:length(list)
    listmin(i)=(list(i)-fnumber);
    if listmin(i)>0
        listmin(i)=100;
    end
    listmin(i)=abs(listmin(i));
end
minimum = min(min(listmin));
[x,y]=find(listmin==minimum)
fnumber=list(y);


    

    
for i=1:length(pixel)
    EV=EV100+log(ISO/100)/log(2);
    t(i)=fnumber/(2^EV);
    deltax(i)=d*pixel(i)*anchosensor/(fnumber*resolucion*1000);
    v(i)=deltax(i)/t(i);

end
figure
plot(pixel,v);

clear all
ISO=linspace(0,25600,1000);
d=1500;
dfar=2000;
f=50;
pixel=1;
anchosensor=23.5;
resolucion=6000;
iluminancia=5378.1512
EV100=log(iluminancia/2.5)/log(2)
%%%%%%%%%%%%%%%%%%%%
list=[3.5 4 4.5 5 5.6 6.3 7.1 8 9 10 11 13 14 16 18 20 22];
fnumber=sqrt(375*(d*f/(d-f)-dfar*f/(dfar-f)));
for i=1:length(list)
    listmin(i)=(list(i)-fnumber);
    if listmin(i)>0
        listmin(i)=100;
    end
    listmin(i)=abs(listmin(i));
end
minimum = min(min(listmin));
[x,y]=find(listmin==minimum)
fnumber=list(y);


    

    
for i=1:length(ISO)
    EV(i)=EV100+log(ISO(i)/100)/log(2);
    t(i)=fnumber/(2^EV(i));
    deltax(i)=d*pixel*anchosensor/(fnumber*resolucion*1000);
    v(i)=deltax(i)/t(i);

end
figure
plot(ISO,v);

clear all
ISO=800;
d=1500;
dfar=2000;
f=50;
pixel=1;
anchosensor=23.5;
resolucion=6000;
EV100=linspace(1,13,1000);
%%%%%%%%%%%%%%%%%%%%
list=[3.5 4 4.5 5 5.6 6.3 7.1 8 9 10 11 13 14 16 18 20 22];
fnumber=sqrt(375*(d*f/(d-f)-dfar*f/(dfar-f)));
for i=1:length(list)
    listmin(i)=(list(i)-fnumber);
    if listmin(i)>0
        listmin(i)=100;
    end
    listmin(i)=abs(listmin(i));
end
minimum = min(min(listmin));
[x,y]=find(listmin==minimum)
fnumber=list(y);


    

    
for i=1:length(EV100)
    EV(i)=EV100(i)+log(ISO/100)/log(2);
    t(i)=fnumber/(2^EV(i));
    deltax(i)=d*pixel*anchosensor/(fnumber*resolucion*1000);
    v(i)=deltax(i)/t(i);

end
figure
plot(EV100,v);


