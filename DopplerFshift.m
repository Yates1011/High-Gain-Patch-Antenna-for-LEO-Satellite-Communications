U = 6.67408*10^-11;
R = 6371*10^3;
H = 640*10^3;
f = 11.7*10^9;
c = 300*10^6;

a = 0:0.1:180;
alpha = deg2rad(a);
Beta = (asin((R*cos(alpha))/(R+H))); %0.497656838 radians ~28.5136 deg
Vc = sqrt((U/R+H));
Vr = Vc * sin(Beta);
T1 = (sin((1.5708+alpha))/(R+H));
T2 = (sin(Beta)/R);

fdopp = (((Vc*f)/c)*(R*cos(alpha)/(R+H)));


figure(1)
plot(a,fdopp/1000)    
title('Doppler shift with respect to satellite elevation');
xlabel('Elevation in degrees')
ylabel('Frequecy shift in MHz')



