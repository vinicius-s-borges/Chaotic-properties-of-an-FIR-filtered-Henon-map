clc
clear all

c = poly([exp(1j*pi/4) exp(-1j*pi/4)]);
c = poly([exp(1j*pi/4) exp(-1j*pi/4)])/sum(c);

figure
subplot(1,2,1)
zplane(c)
hold on
th = linspace( pi/4, 0, 100);
R = 0.5;  
x = R*cos(th) + 0;
y = R*sin(th) + 0;
plot(x,y);
x = [0 0.7071];
y = [0 0.7071];
line(x,y,'Color','red','LineStyle','--')
x = [0 0.7071];
y = [0 -0.7071];
line(x,y,'Color','red','LineStyle','--')
hold off
ax.XLim = [-1 1];
ax.YLim = [-1 1];
subplot(1,2,2)
[h,w] = freqz(c,1,'whole',4001);
plot(w/pi,0.1*log10(abs(h)))
ax.XLim = [0 1];
ax.YLim = [0 -40];
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')