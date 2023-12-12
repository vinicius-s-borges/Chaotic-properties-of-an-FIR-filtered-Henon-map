clc
clear all

c = poly(-ones(1,1));
c = poly(-ones(1,1))/sum(c);

figure
subplot(1,2,1)
zplane(c)
ax.XLim = [-1 1];
ax.YLim = [-1 1];
subplot(1,2,2)
[h,w] = freqz(c,1,'whole',2001);
plot(w/pi,20*log10(abs(h)))
ax.XLim = [0 1];
ax.YLim = [0 -40];
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')