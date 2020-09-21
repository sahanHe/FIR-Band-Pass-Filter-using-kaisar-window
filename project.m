close all;
clear all;
index_no=[1 6 0 1 9 8];
%constructing Filter parameters using index number
A_p=0.05+(0.01*index_no(4));                                    %maximum Passband attenuation
A_a=40+index_no(5);                                             %minimum stopband attenuation
wp1=index_no(6)*100 +400;                                       %Lower passband edge
wp2=index_no(6)*100 +700;                                       %Upper passband edge
wa1=index_no(6)*100 +250;                                       %Lower stopband edge
wa2=index_no(6)*100 +800;                                       %Upper stopband edge
ws=2*(index_no(6)*100 +1200);                                   %sampling frequancy
L=4048;                                                         %number of samples for input signals
Bt=min((wp1-wa1),(wa2-wp2));                                    %Critical transition width
wc1=wp1-(Bt/2);                                                 %Lower cutoff frequency
wc2=wp2+(Bt/2);                                                 %Higher cuttoff drequency
delta_p=(10^(0.05*A_p) -1)/(10^(0.05*A_p) +1);                  %passband ripple 
delta_a=10^(-0.05*A_a);                                         %stopband ripple
delta=min(delta_p,delta_a);                                     %actual ripple

Aa=-20*log10(delta);                                            %Actual Attenuation
alpha=.5842*(Aa-21)^(.4)+0.07886*(Aa-21);                       
D=(Aa-7.95)/14.36;
N=(ws*D/Bt)+1;                                                  %calculating N
if mod(fix(N),2)==0
    N=N+1;
else
    N=N+2;
end
N=fix(N);
M=(N-1)/2;                                                      
n_vec=linspace(-M,M,2*M+1);                                     %vector of N samples between -N-1/2 and N+1/2  
B=alpha*sqrt(1-(n_vec./M).^2);                                  %calculatng beta
w_knt=Bessel1(B)/getBessel(alpha);                              %Kaisar window generation 
T=2*pi/ws;                                                      %Sampling period
h_n=sin(wc2*n_vec*T)./(pi*n_vec)-sin(wc1*n_vec*T)./(pi*n_vec);  %impulse response of the Band pass filter
h_n(M+1)=(wc2-wc1)*2/ws;                                        
hw=w_knt.*h_n;                                                  %FIR Filter

figure;                                                         %Plotting kaisar window
stem(n_vec,w_knt);
xlabel('n');
ylabel('amplitude');
title('kaisar window')
hold on;

figure;                                                        %plotting the impulse response of the FIR Filter
stem(n_vec,hw);
xlabel('n');
ylabel('amplitude');
title('impulse response of the windowd FIR Filter')
hold on;

[Hw,w]=freqz(hw);                                               %Frequency response of the filter
w=w*(ws/(2*pi));                                                %freq domain correction
H_db=20*log10(abs(Hw));                                         %get the magnitude in db scale
start_freq=round((2*wc1/ws)*512);                               %starting frequancy index of for the passband
end_freq=round((2*wc2/ws)*512);                                 %ending frequancy index for the passband

PB_rippleU=0.06.*ones(length(w(start_freq:end_freq+10)));       %maximum passband ripple
PB_rippleL=-0.06.*ones(length(w(start_freq:end_freq+10)));

figure;                                                         %plotting the passband frequancy response
plot(w(start_freq:end_freq+10),H_db(start_freq:end_freq+10));
xlabel('frequancy(rad/s)');
ylabel('gain(dB)');
title('Frequancy response of the Filter in passband')

hold on;
plot(w(start_freq:end_freq+10),PB_rippleU);
hold on;
plot(w(start_freq:end_freq+10),PB_rippleL);
u_step= zeros(1,117,'double');
hw_c=hw(M+1:2*M+1);
figure;                                                         %plotting the causal impulse response of the FIR filter
stem(hw_c);
xlabel('n ');
ylabel('Amplitude');
title('causal impulse response')
hold on;
figure;
plot(w,abs(Hw));                                                %plotting FIR filter frequancy response in normal range 
hold on;

%Testing for Sinusoids input

n_vec2=linspace(-L/2,L/2 -1,L);                                 %sample vector sinusoids
Sins=sin(wc1*(pi/ws)*n_vec2)+sin((wc1+wc2)*(pi/ws)*n_vec2)+sin(((ws/2)+wc2)*(pi/ws)*n_vec2);    %input sinusoids
figure;                                                         %plotting input signal
stem(n_vec2,Sins);                                              
xlabel('n');
ylabel('Amplitude');
title('input signal')
hold on;
f=[-L/2:(L-2)/2]*(ws/L);                                        %frequancy range for frequancy domain
Sinsh2=fftshift(fft(Sins,L));                                   %Frequancy domain representation of input

figure;                                                         %ploting input in frequancy domain
plot(f,abs(Sinsh2)/(L));                
xlabel('rad/s');
ylabel('amplitude');
title('input signal in frequancy domain')
hold on;
H4=fftshift(fft(hw,L));                                         %FFT of FIR filter
figure;                                                         %ploting frequancy response of FIR filter
plot(f,abs(H4));
xlabel('rad/s');
ylabel('amplitude');
title('FIR filter in frequancy domain')
hold on;
filt_Sins2=Sinsh2.*H4;                                          %filtering the input signal
figure;                                                         %plotting filtered Inut in frequncy domain
plot(f,(1/L)*abs(filt_Sins2));
xlabel('rad/s');
ylabel('amplitude');
title('Output signal from FIR filter')
hold on;

figure;                                                         %plotting Ideal Filter OUTPUT
plot(f,(1/L)*abs(fftshift(fft(sin((wc1+wc2)*(pi/ws)*n_vec2),L))));
xlabel('rad/s');
ylabel('amplitude');
title('Output signal from Ideal filter')
hold on;













    
    
   
    
    



