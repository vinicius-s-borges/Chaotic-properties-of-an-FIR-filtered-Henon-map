clc
clear all

%% Temporal series 
a=1.4;b=0.3;N=3000; %Map parameters 
x=zeros(1,N);y=zeros(1,N); 
x(1)=0;y(1)=0; %initial conditions
for n=1:N
x(n+1)=a-(x(n))^2+b*y(n);
y(n+1)=x(n);
end

x1=zeros(1,N);y1=zeros(1,N);
x1(1)=0.0001;y1(1)=0.0001; %initial conditions
for n=1:N
x1(n+1)=a-(x1(n))^2+b*y1(n);
y1(n+1)=x1(n);
end

%% phase space

N = 1e6; % number of points

xh=zeros(1,N);yh=zeros(1,N);
xh(1)=0;yh(1)=0;
for n=1:N
xh(n+1)=a-(xh(n))^2+b*yh(n);
yh(n+1)=xh(n);
end

%% DEP
Nfft = 1e4; % number of FFT points
sinal = xh(1,1:end);
Num_Sinais = 100; %number of signals  
[pxx1] = pwelch(sinal,64,32);
w = linspace(0,1,2*64-1);

%%
pontos = 0:1:50;

figure
subplot(2,2,[1 2])
plot(pontos,x(1:51),'Color','black','LineWidth',2);
hold on
plot(pontos,x1(1:51),'--','Color','black','LineWidth',2);
ylabel('$x_{1}$','Interpreter','latex','FontSize', 40)
xlabel('$n$','Interpreter','latex','FontSize', 40)
grid on
set(gca,'FontSize',20,'LineWidth',2)
subplot(2,2,3)
plot(x(500:3000),y(500:3000),'.','Color','black','LineWidth',2);
ylabel('$x_{2}$','Interpreter','latex','FontSize', 40)
xlabel('$x_{1}$','Interpreter','latex','FontSize', 40)
grid on
set(gca,'FontSize',20,'LineWidth',2)
subplot(2,2,4)
plot(w,pxx1(2:end-1)/max(pxx1(2:end-1)),'Color','black','LineWidth',2); %ylim([-0.1 3]);grid;
ylabel('$$S_{XX}(\omega)$$','Interpreter','Latex','FontSize',20);
xlabel('$$\omega/\pi$$','Interpreter','Latex','FontSize',20);
grid on;
ylim([0 1.1])
set(gca,'FontSize',20,'LineWidth',2)

