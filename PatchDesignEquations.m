er = 2.2;
Fr = 11.95*10^9;
h = 0.158;
c = 300*10^6;
t = 35*10^-6;

W = (c/(2*Fr))*(sqrt(2/(er+1))); %Answer in metres *100 for cm
Wcm = W*100;
ereff = ((er+1)/2)+((er-1)/2)*(1+12*(h/Wcm))^-0.5;
Lambda = c/Fr
%Lambda = c/(Fr*sqrt(ereff));
LambdaDiv2 = (Lambda/2);
%LambdaDiv4 = Lambda/4;
WoverH = Wcm/0.158;


EffL = (0.412 * h)*((ereff+0.3)*(((Wcm/h)+0.264))/((ereff-0.258)*((Wcm/h)+0.8))) %Answer in metres *100 for cm
EffLmm = EffL*10;
L = (c/(2*Fr*(sqrt(ereff)))) - (2*EffL) %Answer in metres *100 for cm
Le = L + (2*EffL)

Rin = 90*(((er)^2)/(er-1))*(L/W)
Wo = (7.48*L)/(exp(50*((sqrt(er + 1.41))/87))) - 1.25*t;
Zc = (120*pi)/((sqrt(ereff))*((Wo/L)+1.393+0.667*log(Wo/L+1.444)))

Z1 = sqrt(Zc*Rin)

Z0 = (120*pi*ereff)/(WoverH + 1.393 + 0.667*log(WoverH + 1.444));
