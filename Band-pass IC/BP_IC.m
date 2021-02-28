%Kostopoulou Natalia 
%AEM 9146

%BP Inverse Chebyshev
clear;
clc;
% initialization variables
a1 = 9;
a2 = 1;
a3 = 4;
a4 = 6;

f0=1000;
f1=650 + 25*a4;
f2= f0^2/f1;
D= 2.1*(f0^2-f1^2)/f1;
f3= (-D + sqrt(D^2+4*f0^2))/2;
f4=f0^2/f3;

amin = 35 - a3; 
amax = 0.4 + a4/36;

% prodiagrafes.
w1 = 2 * pi * f1;
w2= 2 * pi * f2;
w3= 2* pi * f3;
w4 = 2 * pi * f4;
wo = 2 * pi * f0;

Wp = 1;
Ws = ( w4 - w3 ) / ( w2 - w1 );
bw = w2 - w1 ;
qc = wo / bw;

%taksi filtrou
n = acosh( sqrt( ( 10^( amin / 10 ) -1 ) / ( 10^( amax / 10 ) - 1 ) ) ) / acosh( Ws );
n = ceil(n);
%suntelestis diakimanshs
e=1/sqrt(10^(amin/10)-1);

a = ( 1 / n ) * ( asinh( 1 / e ) );
whp = 1 / cosh( acosh( 1 / e ) / n);

%gwnies butterworth gia n=4
y1=22.5;
y2=-y1;
y3=67.5;
y4=-y3;
%real
s1=-sinh(a)*cosd(y1);
s3=-sinh(a)*cosd(y3);
%imaginary
u1=cosh(a)*sind(y1);
u3=cosh(a)*sind(y3);

%poloi
p1 = s1 + 1i*u1;
p2 = s1 - 1i*u1;
p3 = s3 + 1i*u3;
p4 = s3 - 1i*u3;

W12 = norm(p1);
W34 = norm(p3);
Q12 = W12/(-2*real(p1));
Q34 = W34/(-2*real(p3));

% Norma twn polwn Inverse Chebyshev
InvW12 = 1 / W12;
InvW34 = 1 / W34;

% Inverse Chebyshev poles.
InvW1= InvW12 * Ws ;
InvW2= InvW34 * Ws;

%Thesi polwn
S12 = -InvW1 / ( 2 * Q12 );
S34= -InvW2 / ( 2 * Q34);
W12 = sqrt( InvW1^2 - S12^2 );
W34 = sqrt( InvW2^2 - S34^2 );

% Zeros.
z1 = sec(pi/(2*n));
z2 = sec(3*pi/(2*n));
%klimakopoihsh mhdenikwn
Z1 = z1 * Ws;
Z2= z2 * Ws;

%metasximatismos polou 1,2
p1Inv = S12 + ( W12 * 1i );
C1= S12^2 + W12^2;
D1= -2* S12 / qc;
E1= 4 + C1/ qc^2;
G1= sqrt (E1^2 - 4* D1^2);
Q1_2= 1/D1 * sqrt ( 1/2* ( E1+ G1) );
k1= (-S12 * Q1_2 ) /qc;
W1= k1 + sqrt( k1^2 -1);
wo1 = 1/W1 * wo;
wo2 = W1* wo;

%metasximatismos polou 3,4
p3Inv = S34 + ( W34 * 1i );
C2= S34^2 + W34^2;
D2= -2* S34 / qc;
E2= 4 + C2/ qc^2;
G2= sqrt (E2^2 - 4* D2^2);
Q3_4= 1/D2 * sqrt ( 1/2* ( E2+ G2) );
k2= -S34 * Q3_4 /qc;
W2= k2 + sqrt( k2^2 -1);
wo3 = 1/W2 * wo;
wo4 = W2 * wo;

%metasximatismos zero1
Kzero1 = 2 + (Z1^2) / (qc^2);
x1 = ( Kzero1 + sqrt( Kzero1^2 - 4 ) ) / 2;
wz1 = wo * ( sqrt(x1) );
wz2 = wo / ( sqrt(x1) );

%metasximatismos zero2
Kzero2 = 2 + (Z2^2) / (qc^2);
x2 = ( Kzero2+ sqrt( Kzero2^2 - 4 ) ) / 2;
wz3 = wo * ( sqrt(x2) );
wz4 = wo / ( sqrt(x2) );

%UNIT I LPN
wzo1= wz1/wo1;
R11=1;
R14=1;
C11=1/(2*Q1_2);
R12=4*(Q1_2^2);
R15=R12/(wzo1^2-1);
R13=wzo1^2/(2*Q1_2^2);
k1=1/(R13+1);
%klimakopoihsh UNIT I
kf1=wo1;
km1= C11 / ( kf1 * 10^(-6));
R11new=R11*km1;
R12new=R12*km1;
R13new=R13*km1;
R14new=R14*km1;
R15new=R15*km1;

%UNIT II HPN
wzo2= wz2/ wo2;
R21=1 ;
R23=1;
k21 = 1/ wzo2^2 -1 ;
k22 = (2 + k21) * Q1_2^2 / ((2+k21) * Q1_2^2 +1) ;
R22= (2 + k21) ^2 * Q1_2^2;
R24= (2 + k21) * Q1_2 ^2;
C22 = 1/ ( (2+k21) * Q1_2 );
C21= k21 * C22;
k2 = k22 * 1/ wzo2^2;
%klimakopoihsh UNIT II
kf2=wo2;
km2= C22/ ( kf2 * 10^(-6));
R21new= R21 * km2;
R22new= R22 * km2;
R23new= R23 * km2;
R24new= R24 * km2;
C21new= C21/ (km2 * kf2);

%UNIT III LPN
wzo3= wz3/wo3;
R31=1;
R34=1;
C31=1/(2*Q3_4);
R32=4*(Q3_4^2);
R35=R32/(wzo3^2-1);                    
R33=wzo3^2/(2*Q3_4^2);
k3=1/(R33+1);
%klimakopoihsh UNIT III
kf3=wo3;
km3= C31/ ( kf3 * 10^(-6));
R31new=R31*km3;
R32new=R32*km3;
R33new=R33*km3;
R34new=R34*km3;
R35new=R35*km3;

%UNIT IV HPN
wzo4= wz4/ wo4;
R41=1 ;
R43=1;
k41 = (1/wzo4^2) -1 ;
k42 = (2 + k41) * Q3_4^2 / ((2+k41) * Q3_4^2 +1) ;
R42= (2 + k41) ^2 * Q3_4^2;
R44= (2 + k41) * Q3_4^2;
C42 = 1/ ( (2+k41) * Q3_4 );
C41 = k41 * C42;
k4 = k42 * 1/ wzo4^2;
%klimakopoihsh UNIT IV
kf4=wo4;
km4= C42/ ( kf4 * 10^(-6));
R41new= R41 * km4;
R42new= R42 * km4;
R43new= R43 * km4;
R44new= R44 * km4;
C41new= C41 / (km4 * kf4);

%synartiseis metaforas
T1 = tf( [k1 0 ( k1 * wz1^2 ) ], [ 1 ( wo1 / Q1_2 ) wo1^2 ] );
T2 = tf( [k2 0 ( k2 * wz2^2 ) ], [ 1 ( wo2 / Q1_2 ) wo2^2 ] );
T3 = tf( [k3 0 ( k3 * wz3^2 ) ], [ 1 ( wo3 / Q3_4 ) wo3^2 ] );
T4 = tf( [k4 0 ( k4 * wz4^2 ) ], [ 1 ( wo4 / Q3_4 ) wo4^2 ] );

%Overall TF
TBP = series( T1, T2 );
TBP = series (TBP, T3);
TBP = series (TBP, T4);

%rithmisi sta 10 db
totalGain = abs(evalfr(TBP, wo * 1i));
totalGain1 = abs(evalfr(T1, wo * 1i));
totalGain2 = abs(evalfr(T2, wo * 1i));
totalGain3 = abs(evalfr(T3, wo * 1i));
totalGain4 = abs(evalfr(T4, wo * 1i));

again = ( 10^ ( 0.5 ) )/totalGain;
TBPf = again * TBP;

InverseTBP = inv(TBPf);
% plot_transfer_function(TBP, [f1 f2 f3 f4])
% plot_transfer_function(TBPf, [f1 f2 f3 f4])
% plot_transfer_function(T1, [f1 f2 f3 f4])
% plot_transfer_function(T2, [f1 f2 f3 f4])
% plot_transfer_function(T3, [f1 f2 f3 f4])
% plot_transfer_function(T4, [f1 f2 f3 f4])
 plot_transfer_function(InverseTBP, [f1 f2 f3 f4])
 ltiview ({'bodemag'}, T1,T2,T3,T4,TBPf)

% eisagw periodiko shma
f11= (wo- (wo-w1)/3);
f12=(wo+(wo+w1)/4);
f13= 0.5*w3;
f14=(2.4*w4);
f15= (3*w4);

syms t;
input_signal(t) = cos(f11*t) + 0.6*cos(f12*t) + 0.7*cos(f13*t)+ 0.8*cos(f14*t) + 0.6*cos(f15*t);
figure;

time_sym = 0:0.001:0.15;
ezplot(input_signal(t),time_sym);
Fs = 20000;                              
T = 1/Fs;                  
L = 1500;           
time = (0:L-1)*T;

x = cos(f11*time ) + 0.6*cos(f12*time ) + 0.7*cos(f13*time )+ 0.8*cos(f14*time ) + 0.6*cos(f15*time );

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

sim = lsim(TBPf,x,time);
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

