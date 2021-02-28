%Kostopoulou Natalia 
%AEM 9146

%LP Butterworth
clear;
clc;
% initialization variables
a1 = 9;
a2 = 1;
a3 = 4;
a4 = 6;
m = 2;

fp = 0.6 * (3 + m) * 1000;
fs = 2.0 * fp;
ws = 2 * pi * fs;
wp = 2 * pi * fp;
amin = 17.5 + (max(1,a4) - 5) * 0.5;
amax = 0.6 + (max(1,a3) - 3) / 10;

nt = (log((10^(amin/10)-1) / ((10^(amax/10)-1)))) / (2*log(ws/wp));
n = ceil(nt);
w0 = wp / (((10^(amax/10)-1)) ^ (1/(2*n)));

%For n = 5
y1 = 0;
y2 = 0.628319;
y3 = 1.25664;

Q1 = 0.5;
Q2 = 0.62;
Q3 = 1.62;

p1 = -cos(y1) + 1i * sin(y1);
p2 = -cos(y2) + 1i * sin(y2);
p3 = -cos(y3) + 1i * sin(y3);

kf = w0;
Cnew = 0.1 * 10^(-6);
%first unit
C1 = 1;
R1 = 1;
km1 = C1 / (kf * Cnew);
C1new = Cnew;
R1new = km1 * R1;

% second unit
k = 1;
R21 = 1;
R22 = 1;
C21 = 2 * Q2;
C22 = 1/(2 * Q2);

km2 = C21 / (kf * Cnew);
C21new = Cnew;
C22new = C22 / (kf * km2);
R21new = R21 * km2;
R22new = R22 * km2;  

% third unit
k = 1;
R31 = 1;
R32 = 1;
C31 = 2 * Q3;
C32 = 1/(2 * Q3);

km3 = C31 / (kf * Cnew);
C31new = Cnew;
C32new = C32 / (kf * km3);
R31new = R31 * km3;
R32new = R32 * km3;

% transfer function
T1 = tf( w0 , [1 w0] );
T2 = tf( w0^2 , [1 w0/Q2 w0^2] );
T3 = tf( w0^2 , [1 w0/Q3 w0^2] );

Tlp = T1 * T2 * T3;

% ltiview({'bodemag'}, T1)
% ltiview({'bodemag'}, T2)
% ltiview({'bodemag'}, T3)
% ltiview ({'bodemag'}, Tlp)
% ltiview({'bodemag'}, T1, T2, T3, Tlp)
% plot_transfer_function(Tlp , [fs fp] );

% attenuation transfer function
 invTlp = inv(Tlp);
 plot_transfer_function(invTlp, [fs fp] );

%fourier 
T = 10*(1/2000);
t = 0:1/1000000:T-1/1000000;
x = sawtooth(2*pi*2000*t,0.5);
Fs=40000;

figure;
plot(t,x);
grid
% fasma eisodou
X=fft(x);
n=1500;
f = Fs*(0:(n/2))/n;
Px = abs(X/n);
figure
plot(f,Px(1:n/2+1));
grid
title('Frequency Domain of input')
xlabel('Frequency (f)')
ylabel('|Yx|')

% fasma eksodou
y = lsim(Tlp,'r',x,t);
Y=fft(y);
Py = abs(Y/n);
figure
plot(f,Py(1:n/2+1))
title('Frequency Domain of output')
xlabel('Frequency (f)')
ylabel('|Yx|')
grid on

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

