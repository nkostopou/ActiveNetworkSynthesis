%Kostopoulou Natalia 
%AEM 9146

%BE Chebyshev
clear;
clc;
% initialization variables
a1=9;
a2=1;
a3=4;
a4=6;

f0=1.05*1000; 
f1=750 +(25*a4); 
f2 = (f0^2)/f1;
D = (1/1.8)*((f0^2-f1^2)/f1);
f3 = (-D + sqrt(D^2 +4*f0^2))/2;
f4 = (f0^2)/f3;
amin= 25 +((a3-5)/10); 
amax= 0.5 + (a4/10) ;

%metatropi suxnotitwn
w0 = 2*pi*f0;
w11 = 2*pi*f1;
w21 = 2*pi*f2;
w31 = 2*pi*f3;
w41 = 2*pi*f4;
bw=w21-w11;

Wp=1;
Ws=(w21-w11)/(w41-w31);

%taksi filtrou
n = acosh(sqrt((10^(amin/10) -1)/(10^(amax/10) -1)))/acosh(Ws);
n=ceil(n);

%suntelestis diakimanshs
e=sqrt(10^(amax/10)-1);
a=(1/n)*asinh(1/e);

%suxnothta hmisews isxuos 
whp = cosh((1/n)*(acosh(1/e)));

%gwnies butterworth gia n=4
y1=22.5;
y2=-y1;
y3=67.5;
y4=-y3;
%real
s1=-sinh(a)*cosd(y1);
s3=-sinh(a)*cosd(y3);
%imaginary
w1=cosh(a)*sind(y1);
w3=cosh(a)*sind(y3);

%poloi
p1 = s1 + 1i*w1;
p2 = s1 - 1i*w1;
p3 = s3 + 1i*w3;
p4 = s3 - 1i*w3;

%metro polwn
W01 = norm(p1);
W03 = norm(p3);

%antistrofi polwn
wo1 = (1/W01);
wo3 = (1/W03);

Q1 =W01/(2*abs(s1));
Q3 =W03/(2*abs(s3));

%poloi inverse chebyshev
p11 = wo1*exp(1i*angle(p1));
p22 = wo1*exp(1i*angle(p2));
p33 = wo3*exp(1i*angle(p3));
p44 = wo3*exp(1i*angle(p4));

%Geffe polou 1,2

Sigma12 = abs(real(p11));
W12 = abs(imag(p11));
qc= w0/bw;

C1 = Sigma12^2 + W12^2;
D1 = (2*Sigma12)/qc;
E1 = 4 + (C1/(qc^2));
G1 = sqrt(E1^2 - (4*(D1^2)));
Q1g=(1/D1)*sqrt(0.5*(E1+G1));
k1 = (Sigma12*Q1g)/qc;
W1g = k1 + sqrt(k1^2 -1);
w2g = W1g*w0;      %w02
w1g = (1/W1g)*w0;  %w01


%Geffe polou 3,4
Sigma34 = abs(real(p33));
W34 = abs(imag(p33));

C3 = Sigma34^2 + W34^2;
D3 = (2*Sigma34)/qc;
E3 = 4 + (C3/(qc^2));
G3 = sqrt(E3^2 - (4*(D3^2)));
Q3g=(1/D3)*sqrt(0.5*(E3+G3));
k3 = (Sigma34*Q3g)/qc;
W3g = k3 + sqrt(k3^2 -1);
w4g = W3g*w0;     %w04
w3g = (1/W3g)*w0; %w03

%Monada I LPN (w0,w1g)
Wz1=w0/w1g;
R11=1;
R14=1;
R12=4*(Q1g^2);
R13=Wz1^2/(2*(Q1g^2));
R15 = 4*Q1g^2/(Wz1^2 -1);
C11 = 1/(2*Q1g);
k1H = 1/(1+R13); 
k1L = k1H*((w0/w1g)^2);

%klimakopoihsh
kf1 = w1g;
km1 = C11/(kf1*10^(-6));

R11new=R11*km1;
R12new=R12*km1;
R13new=R13*km1;
R14new=R14*km1;
R15new=R15*km1;

%UNIT II HPN
Wz2=w0/w2g;
R21=1 ;
R23=1;
k21 = (1/ Wz2)^2 -1 ;
k22 = (2 + k21) * Q1g^2 / ((2+k21) * Q1g^2 +1) ;
R22= (2 + k21) ^2 * Q1g^2;
R24= (2 + k21) * Q1g ^2;
C22 = 1/ ( (2+k21) * Q1g );
C21= k21 * C22;
k2H = k22 * (1/ Wz2^2);
k2L = k2H*(Wz2^2);
%klimakopoihsh
kf2 = w2g;
km2= C22/ ( kf2 * 10^(-6));
R21new= R21 * km2;
R22new= R22 * km2;
R23new= R23 * km2;
R24new= R24 * km2;
C21new= C21/ (km2 * kf2);
C22new= C22/ (km2 * kf2);

%UNIT III LPN
Wz3= w0/w3g;
R31=1;
R34=1;
C31=1/(2*Q3g);
R32=4*(Q3g^2);
R35=R32/(Wz3^2-1);                    
R33=Wz3^2/(2*Q3g^2);
k3H=1/(R33+1);
k3L= k3H*(Wz3^2);
%klimakopoihsh UNIT III
kf3=w3g;
km3= C31/ ( kf3 * 10^(-6));
R31new=R31*km3;
R32new=R32*km3;
R33new=R33*km3;
R34new=R34*km3;
R35new=R35*km3;

%UNIT IV HPN
Wz4= w0/w4g;
R41=1 ;
R43=1;
k41 = (1/Wz4^2) -1 ;
k42 = (2 + k41) * Q3g^2 / ((2+k41) * Q3g^2 +1) ;
R42= (2 + k41) ^2 * Q3g^2;
R44= (2 + k41) * Q3g^2;
C42 = 1/ ( (2+k41) * Q3g );
C41 = k41 * C42;
k4H = k42 * (1/ Wz4^2);
k4L = k4H*(Wz4^2);

%klimakopoihsh UNIT IV
kf4=w4g;
km4= C42/ ( kf4 * 10^(-6));
R41new= R41 * km4;
R42new= R42 * km4;
R43new= R43 * km4;
R44new= R44 * km4;
C41new= C41 / (km4 * kf4);
C42new=C42/(km4 * kf4);

%synartiseis metaforas
T1=k1H*tf([1,0,w0^2], [1, w1g/Q1g, w1g^2]);
T2=k2H*tf([1,0,w0^2], [1, w2g/Q1g, w2g^2]);
T3=k3H*tf([1,0,w0^2], [1, w3g/Q3g, w3g^2]);
T4=k4H*tf([1,0,w0^2], [1, w4g/Q3g, w4g^2]);

%Overall TF
TBE = series( T1, T2 );
TBE = series (TBE, T3);
TBE = series (TBE, T4);

%sunoliko kerdos stis ypshles
ktot=k1H*k2H*k3H*k4H;

%rithmisi kerdous
kfinal = db2mag(10)/ktot;

%teliki synartisi metaforas
TBEf=kfinal*TBE;

%  plot_transfer_function( T1, [f1 f2 f3 f4] );
%  plot_transfer_function( T2, [f1 f2 f3 f4] );
%  plot_transfer_function( T3, [f1 f2 f3 f4] );
%  plot_transfer_function( T4, [f1 f2 f3 f4] );
%  plot_transfer_function( TBE, [f1 f2 f3 f4] );
 plot_transfer_function( TBEf, [f1 f2 f3 f4] );
%  plot_transfer_function( inv(TBEf), [f1 f2 f3 f4] );
%  ltiview ({'bodemag'}, T1,T2,T3,T4,TBEf);

% eisagw periodiko shma
f11=(w0- (w0-w31)/2);
f12=(w0+(w0+w31)/2);
f13= 0.5*w11;
f14=(2.4*w21);
f15= (3.5*w21);

syms t;
input_signal(t) = 0.8*cos(f11*t) + 0.1*cos(f12*t) + cos(f13*t)+ 0.8*cos(f14*t) + 0.4*cos(f15*t);
figure;

time_sym = 0:0.001:0.15;
ezplot(input_signal(t),time_sym);
Fs = 20000;                              
T = 1/Fs;                  
L = 1500;           
time = (0:L-1)*T;

x = 0.8*cos(f11*time ) + 0.1*cos(f12*time ) + cos(f13*time )+ 0.8*cos(f14*time ) + 0.4*cos(f15*time );

Y = fft(x);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

sim = lsim(TBEf,x,time);
Y = fft(sim);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

function plot_transfer_function( tf, frequency_markers )
%PLOT_TRANSFER_FUNCTION Plots bode of a transfer function with markers
%
%   tf                - The transfer function (created using tf)
%   frequency_markers - A matrix of frequencies in Hz
%
%   Example:
%       plot_transfer_function( tf([1000], [1 1000]), [10 1000 10000] );

figure;
x_space = logspace(1,5,5000); % 5000 points between 10^1 and 10^5
x_space = 2 * pi * x_space; % to rad / sec
[mag,~,wout] = bode(tf,x_space);
mag = squeeze(mag);
wout = squeeze(wout);
mag = 20*log10(mag);
wout = wout/2/pi;
semilogx(wout,mag,'-b');
axis([min(wout) max(wout) min(mag)-10 max(mag)+10]);
[num,den] = tfdata(tf);
syms s;
d1 = digits(5);
ltx = latex(vpa(poly2sym(cell2mat(num),s)/poly2sym(cell2mat(den),s)));
digits(d1);
title(strcat('$',ltx,'$'), 'Interpreter','latex', 'FontSize', 24);
xlabel('Frequency (Hz)', 'FontSize', 18);
ylabel('Magnitude (dB)', 'FontSize', 18);
grid on;
hold on;
[dbMarks,~,frequency_markers] = bode(tf,2 * pi * frequency_markers);
dbMarks = squeeze(dbMarks);
frequency_markers = squeeze(frequency_markers);
dbMarks = 20*log10(dbMarks);
frequency_markers = frequency_markers/2/pi;
Aw = cell(size(frequency_markers, 1) + 1, 1);
Aw{1} = 'Transfer function';
for i = 1 : size(frequency_markers, 1)
    semilogx(frequency_markers(i),dbMarks(i),'o');
    Aw{i+1} = sprintf('Attenuation at %.2f Hz is %.2f dB', ...
        frequency_markers(i), dbMarks(i));
end
legend(Aw,'Location','best','FontSize',12);
set(gca,'FontSize',14);
end



