function [FFTS] = ffts(FE,X,sig)

dt = FE.time(3);   sf = 1/dt;   lsig=length(sig(1,:))

t = ( 0 : lsig-1 ) * dt;


NFFT = 2^nextpow2(lsig); % Next power of 2 from length of signal

ampw = fft(sig,NFFT)/lsig;

freq = sf / 2 * linspace(0,1,NFFT/2+1);

length(ampw)
length(freq)

plot(freq,2*abs(ampw(1:NFFT/2+1)),'.-') 
