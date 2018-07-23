%Simulation of the frequency response of SAW delay line using the
%transmission matrix approach.
close all;
Vs=3488;  %SAW velocity on the free sections of YZ-LiNb03 substrate
Vm=3353.0924;  %SAW velocity on the metallized sections of YZ-LiNb03 substrate
%with  metallization height 200nm and wavelength 34.88um
Co=4.6e-10; %Capacitance per finger pair per unit length on YZ-LiNb03 substrate
K2=0.045;  %Electromechanical coupling coefficient for YZ-LiNb03 substrate
lambdao=34.88e-6;  %SAW wavelength
fo=Vs/lambdao;
W=80*lambdao; %Acoustic aperture
d=127.75*lambdao;  %Length of the delay line
fm=Vm/lambdao;
df=0.125*lambdao;
dm=0.25*lambdao;
Zo=1/(fo*Co*W*K2);  %acoustic impedance for the free sections
Zm=1/(fm*Co*W*K2); %acoustic impedance for the metallized sections
%Zm=169.66e3;
Nt1 =80;  %Number of electrodes in the input IDT
Nt3=30;  %Number of electrodes in the output IDT
h=200e-9;
k11=(0.018+(0.3*h/lambdao))*(2*pi/lambdao); %Self coupling coefficient for YZ-LiNb03 substrate
i=1;
for f=93.5e6:10e6/799:103.5e6
   lambda=Vs/f;
   omega=2*pi*f;


   %computing the ABCD matrix for a single finger
   thetaf=((2*pi*f*df)/(Vs));
   thetam=((2*pi*f*dm)/(Vm));
   Af=cos(thetaf);
   Bf=sqrt(-1)*Zo*sin(thetaf);
   Cf=sqrt(-1)*sin(thetaf)/Zo;
   Df=cos(thetaf);
   Am=cos(thetam);
   Bm=sqrt(-1)*Zm*sin(thetam);
   Cm=sqrt(-1)*sin(thetam)/Zm;
   Dm=cos(thetam);
   Afinger=[Af Bf;Cf Df]*[Am  Bm;Cm  Dm]*[Af Bf;Cf Df];
   Ase=Afinger(1,1);
   Bse=Afinger(1,2);
   Cse=Afinger(2,1);
   Dse=Afinger(2,2);
   thetae=acos(Ase);
   Ze=Bse/(sqrt(-1)*sin(thetae));
   %computing the 2x2 matrix for a single finger in the  IDT
   t11 =0.5*(2*Ase+(Bse/Zo)+Zo*Cse);
   t12=0.5*(Zo*Cse-(Bse/Zo));
   t13=((sqrt(-1 )*tan(thetae/2)*(Zo^0.5))/(2*Ze))*(-Ase-1 -(Bse/Zo));
   t21=-t12;
   t22=conj(t11);
   t23=sqrt(-1)*tan(thetae/2)*(Zo^0.5)*(1 +Ase-(Bse/Zo))/(2*Ze);
   t31=2*t13;
   t32=-2*t23;
   t33=sqrt(-1)*omega*Co*W*0.5+sqrt(-1)*2*(tan(thetae/2)/Ze)-sqrt(-1)*(sin(thetae)*(tan(thetae/2)^2))/Ze;
   %computing the 2x2  IDT matrix
   t1=[t11  t12;t21  t22]^Nt1;
   t3=[t11  t12;t21  t22]^Nt3;
   t111=t1(1,1);
   t113=t3(1,1);
   t121=t1(1,2);
   t123=t3(1,2);
   t211=t1(2,1);
   t213=t3(2,1);
   t221=t1(2,2);
   t223=t3(2,2);
   Bp=[t13;t23]+[t11  t12;t21  t22]*[-t13;-t23];
   Cp=[t31  t32]*[t11  t12;t21  t22]+[-t31  -t32];
   t33p=2*t33+[t31  t32]*[-t13;-t23];
   Tp=[t11 t12;t21 t22]^2;
   tauin=[0;0];
   tauprimein=[0 0];
   tauout=[0;0];
   tauprimeout=[0 0];
   t333=(Nt3/2)*t33p;
   t331=(Nt1/2)*t33p;
   %computing t13, t23, t31, t32, t33 values for the overall  IDT
   for i1=1:(Nt1/2)
       tauin=tauin+(Tp^(i1-1))*Bp;
       tauprimein=tauprimein+Cp*Tp^(i1 -1);
       t331 =t331 +((Nt1/2)-i1 )*Cp*Tp^(i1 -1 )*Bp;
   end
       for i2=1 :(Nt3/2)
           tauout=tauout+(Tp^(i2-1 ))*Bp;
           tauprimeout=tauprimeout+Cp*Tp^(i2-1);
           t333=t333+((Nt3/2)-i2)*Cp*Tp^(i2-1)*Bp;
       end
       t131=tauin(1,1);
       t133=tauout(1,1);
       t231=tauin(2,1);
       t233=tauout(2,1);
       t311 =tauprimein(1,1);
       t313=tauprimeout(1,1); 
       t321=tauprimein(1,2);
       t323=tauprimeout(1,2);
       %computing the matrix elements for the delay path
       thetad=2*pi*f*d/Vs;
       Ad=cos(thetad);
       Bd=sqrt(-1)*Zo*sin(thetad);
       Cd=sqrt(-1)*sin(thetad)/Zo;
       Dd=cos(thetad);
       thetaed=acos(Ad);
       Ze=Bd/(sqrt(-1)*sin(thetaed));
       d11 =0.5*(2*Ad+(Bd/Zo)+Zo*Cd);
       d12=0.5*(Zo*Cd-(Bd/Zo));
       d21=-d12;
       d22=0.5*(2*Ad-(Bd/Zo)-Zo*Cd);
       d2=[d11 d12;d21 d22];
       
       m=t1 *d2*t3;
       k=t1 *d2*[t133;t233];
       p=[t311  t321]*d2*t3;
       l=[t311  t321]*d2*[t133;t233];
       %computing y-parameters for the SAW delay line
       y11(i)=t331 -(p(1,1)*t131/m(1,1));
       y21(i)=-t313*t131/m(1,1);
       y12(i)=l(1,1)-(p(1,1)*k(1,1)/m(1,1));
       y22(i)=t333-(t313*k(1,1)/m(1,1));
       
       %computing s21 using y-parameters
       s21(i)=2*y21(i)/((1+y11(i))*(1+y22(i))-y12(i)*y21(i));
       s12(i)=2*y12(i)/((1+y11(i))*(1+y22(i))-y12(i)*y21(i));
       s11(i)=((1-y11(i))*(1+y22(i))+y12(i)*y21(i))/((1+y11(i))*(1+y22(i)-y12(i)*y21(i)));
       s22(i)=((1+y11(i))*(-+y22(i))+y12(i)*y21(i))/((1+y11(i))*(1+y22(i)-y12(i)*y21(i)));
       
       z11=((1+s11(i)*(1-s22(i))+s12(i)*s21(i))/((1-s11(i))*(1-s22(i))-s12(i)*s21(i)));
       z12=2*s12(i)/((1-s11(i))*(1-s22(i)-s12(i)*s21(i)));
       z21=2*s21(i)/((1-s11(i))*(1-s22(i)-s12(i)*s21(i)));
       z22=((1-s11(i)*(1+s22(i))+s12(i)*s21(i))/((1-s11(i))*(1-s22(i))-s12(i)*s21(i)));


       
       
       
       
       
       
       
       
       %modified s-parameter including source and load  impedance
       s21new(i)=-50*z12/(z12^2-(z11+50)*(z22+50));
       i=i+1;
end
f=93.5e6:10e6/799:103.5e6;
figure(1);
[ax,h1,h2]=plotyy(f,20*log10(abs(s21new)),f,angle(s21new)*180/pi);
title('Frequency Response of delay line with  Np1=50');
xlabel('Frequency (MHz)');
set(h1,'LineWidth',2)
set(h2,'LineWidth',0.5)
set(get(ax(1),'Ylabel'),'String','dB')
set(get(ax(2),'Ylabel'),'String','Phase angle(degrees)')

       


       