%Kostopoulou Natalia 
%AEM 9146

%HP Chebyshev
clear;
clc;
%initialization variables
a1=9;
a2=1;
a3=4;
a4=6;

mu = 0;
fp = (3 + mu)*(10^3);
wp = 2*pi*fp;
fs = fp/1.8;
ws = 2*pi*fs;

amin = 25 + a3*(4/9);
amax = 0.5 + a4*(0.25/9);

% transforming specifications
Wp = 1;
Ws = wp/ws;
%taksi filtrou
n = acosh(((10^(amin/10) - 1)/(10^(amax/10)-1))^(1/2))/acosh(Ws);
n = ceil(n);
%suntelestis diakimanshs
e = sqrt(10^(amax/10) - 1);
a = asinh(1/e)/n;


%Sixnotita imiseias isxios
Whp = cosh((1/n)*((acosh(10^(amax/10) - 1))^(-1/2)));

%Gwnies Butterworth gia n=4
ps1 = 22.5;
ps2 = -22.5;
ps3 = 67.5;
ps4 = -67.5;

%Poloi tou Chebyshev prokiptoun ws exis
real1 = -sinh(a)*cosd(ps1);
imag1 = cosh(a)*sind(ps1);

real2 = -(sinh(a)*cosd(ps2));
imag2 = cosh(a)*sind(ps2);

real3 = -sinh(a)*cosd(ps3);
imag3 = -(cosh(a)*sind(ps3));

real4 = -sinh(a)*cosd(ps4);
imag4 = -(cosh(a)*sind(ps4));


W12 = sqrt(real2^2 + imag2^2);
Q12 = W12/(-2*real1);
W34 = sqrt(real4^2 + imag4^2);
Q34 = W34/(-2*real3);

%anwdiavati sinartisi 
whp = wp/Whp;
%poloi 
wmega12 = wp/W12;
wmega34 = wp/W34;

%UNIT 1
C11 = 1;
C12 = 1;
R12 = 1;
R11 = 1;
k1 = 3 - 1/Q12;
r11=1;
r12= 2-1/Q12;

%Klimakopoiisi
kf1 = wmega12;
km1 = C11/(kf1*10^(-6));
R11 = km1*R11;
R12 = km1*R12;
C11 = 10^(-6);
C12 = C11;
r12=km1*r12;
%Monada 2
C21 = 1;
C22 = 1;
R22 = 1;
R21 = 1;
k2 = 3 - 1/Q34;
r21=1;
r22= 2-1/Q34;

%Klimakopoiisi
kf2 = wmega34;
km2 = C21/(kf2*10^(-6));
R21 = km2*R21;
R22 = km2*R22;
C21 = 10^(-6);
C22 = C21;
r22=km2*r22;

%synartiseis metaforas
TH1 = tf([k1 0 0],[1 wmega12/Q12 wmega12^2]);
TH2 = tf([k2 0 0],[1 wmega34/Q34 wmega34^2]);
%Overall TF
THP = series( TH1, TH2 );

%Sinoliko kerdos 
ktot = k1*k2; %ktot = 4.39005523
%rithmisi kerdous sta 10dB
kfinal = db2mag(10)/ktot;
Ra = 10^3;
Rb = kfinal*Ra;

%teliki synartisi metaforas
THPf=kfinal*THP;

%  plot_transfer_function( TH1, [fs fp] );
%  plot_transfer_function( TH2, [fs fp] );
% plot_transfer_function( THP, [fs fp] );
% plot_transfer_function( THPf, [fs fp] );
%plot_transfer_function( inv(THPf), [fs fp] );
% ltiview ({'bodemag'}, TH1,TH2,THPf);

% eisagw periodiko shma
f11=(0.2*ws);
f12=(0.7*ws);
f13= (1.6*wp);
f14=(2.4*wp);
f15= (5*wp);

syms t;
input_signal(t) = cos(f11*t) + 0.6*cos(f12*t) + 1.5*cos(f13*t)+ 0.7*cos(f14*t) + 0.4*cos(f15*t);
figure;

time_sym = 0:0.001:0.15;
ezplot(input_signal(t),time_sym);
Fs = 20000;                              
T = 1/Fs;                  
L = 1500;           
time = (0:L-1)*T;

x = cos(f11*time ) + 0.6*cos(f12*time ) + 1.5*cos(f13*time )+ 0.7*cos(f14*time ) + 0.4*cos(f15*time );

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

sim = lsim(THPf,x,time);
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
