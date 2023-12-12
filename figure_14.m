clc
clear all

Ns=19;
wpassante = 0.5;

c=fir1(Ns-1,wpassante,'low',hamming(Ns));
figure
subplot(1,2,1)
[hz1, hp1, ht1] = zplane(c); % Plot zeros and poles when coefficients are unquantized
set(findobj(ht1, 'Type', 'line'),'Color','k','linestyle','-','linewidth',2);
ax.XLim = [-1 1];
ax.YLim = [-1 1];
set(gca,'layer','top','FontSize',24)
subplot(1,2,2)
[h,w] = freqz(c,1,'whole',2001);
plot(w/pi,20*log10(abs(h)))
ax.XLim = [0 1];
ax.YLim = [0 -40];
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')