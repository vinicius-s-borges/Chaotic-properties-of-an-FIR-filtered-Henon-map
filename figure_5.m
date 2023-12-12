clc
clear all

%Simulation parameters
Nitera=3000;
len_pontos=150;
pontos=0:1:len_pontos;

%Map parameters
a=1.4;
b=0.3;

%%
Ns=0;
ganho = 1.52433;
c = [1 zeros(1,Ns) 1];
c = ([1 zeros(1,Ns) 1]/sum(abs(c)))*ganho;
Ns=length(c);

x1 = 1*ones(Ns,1);
x2 = -1*ones(Ns,1);
for n=1:1:Nitera
    x1(:,n+1) = Henon_N_2(x1(:,n),a,b,c);
end

for n=1:1:Nitera
    x2(:,n+1) = Henon_N_2(x2(:,n),a,b,c);
end

%%
Ns=14;
ganho = 0.5;
c = [1 zeros(1,Ns) 1];
c = ([1 zeros(1,Ns) 1]/sum(abs(c)))*ganho;
Ns=length(c);

xpf = 1*ones(Ns,1);
xpf1 = -1*ones(Ns,1);
for n=1:1:Nitera
    xpf(:,n+1) = Henon_N_4(xpf(:,n),a,b,c);
end

for n=1:1:Nitera
    xpf1(:,n+1) = Henon_N_4(xpf1(:,n),a,b,c);
end

%% 
Ns=17;
ganho = 0.65;
c = [1 zeros(1,Ns) 1];
c = ([1 zeros(1,Ns) 1]/sum(abs(c)))*ganho;
Ns=length(c);

xp = 0.5*ones(Ns,1);
xp1 = -0.5*ones(Ns,1);
for n=1:1:Nitera
    xp(:,n+1) = Henon_N_4(xp(:,n),a,b,c);
end

for n=1:1:Nitera
    xp1(:,n+1) = Henon_N_4(xp1(:,n),a,b,c);
end

%% Phase space
Ns=14;
ganho = 1;
c = [1 zeros(1,Ns) 1];
c = ([1 zeros(1,Ns) 1]/sum(abs(c)))*ganho;
Ns=length(c);

xc = 0.5*ones(Ns,1);
xc1 = 0.52*ones(Ns,1);
for n=1:1:Nitera
    xc(:,n+1) = Henon_N_4(xc(:,n),a,b,c);
end

for n=1:1:Nitera
    xc1(:,n+1) = Henon_N_4(xc1(:,n),a,b,c);
end

%% DEP

N = 1e6; % number of points
Nfft = 1e4; % number of FFT points

for n=1:1:N
    xc(:,n+1) = Henon_N_4(xc(:,n),a,b,c);
end

sinal = xc(1,1:end);
Num_Sinais = 100; %number of signals   
[pxx1] = pwelch(sinal,64,32);
w = linspace(0,1,2*64-1);

%%
figure
subplot(4,2,[1 2])
plot(pontos,xpf(1,1:length(pontos)),'Color','blue','LineWidth',1)
hold on
plot(pontos,xpf1(1,1:length(pontos)),'Color','black','LineWidth',1)
hold off
ylabel('$x_{1}$','Interpreter','latex','FontSize', 40)
grid on
ylim([-2.1 2.1])
yticks([-2 0 2])
xticks([0 25 50 75 100 125 150])
set(gca,'XTickLabel',[],'FontSize',20,'LineWidth',2)
subplot(4,2,[3 4])
plot(pontos,xp(1,1:length(pontos)),'Color','blue','LineWidth',1)
hold on
plot(pontos,xp1(1,1:length(pontos)),'Color','black','LineWidth',1)
hold off
ylabel('$x_{1}$','Interpreter','latex','FontSize', 40)
grid on
ylim([-2.1 2.1])
yticks([-2 0 2])
xticks([0 25 50 75 100 125 150])
set(gca,'XTickLabel',[],'FontSize',20,'LineWidth',2)
subplot(4,2,[5 6])
plot(pontos,xc(1,1:length(pontos)),'Color','blue','LineWidth',1)
hold on
plot(pontos,xc1(1,1:length(pontos)),'Color','black','LineWidth',1)
hold off
ylabel('$x_{1}$','Interpreter','latex','FontSize', 40)
xlabel('$n$','Interpreter','latex','FontSize', 40)
grid on
ylim([-2.1 2.1])
yticks([-2 0 2])
xticks([0 25 50 75 100 125 150])
set(gca,'FontSize',20,'LineWidth',2)
subplot(4,2,7)
plot(xc(1,500:end),xc(2,500:end),'.','Color','black','MarkerSize',1);
grid on
ylim([-2.5 2.2])
xlim([-2.6 2.1])
ylabel('$x_{1}$','Interpreter','latex','FontSize', 40)
xlabel('$x_{2}$','Interpreter','latex','FontSize', 40)
set(gca,'FontSize',20,'LineWidth',2)
subplot(4,2,8)
plot(w,pxx1(2:end-1)/max(pxx1(2:end-1)),'Color','black','LineWidth',2); %ylim([-0.1 3]);grid;
ylabel('$$S_{XX}(\omega)$$','Interpreter','Latex','FontSize',20);
xlabel('$$\omega/\pi$$','Interpreter','Latex','FontSize',20);
grid on;
ylim([0 1.1])
set(gca,'FontSize',20,'LineWidth',2)


function [x] = Henon_N_4(x,a,b,c)

x1=[a-x(3)^2+b*x(2);
    x(1)];   
x2=[c(1)*(a-x(3)^2+b*x(2))+c(2)*x(1)+c(3)*x(2)];
x3=[c(4:end)*x(4:length(c))];
x4=x2+x3;

x5=[x(2)]; 
x6=[x(4:length(c)-1)];

x7=vertcat(x1,x4);
x8=vertcat(x7,x5);
x=vertcat(x8,x6);

end

function [x] = Henon_N_2(x,alpha,beta,c)
x=[alpha-(c(1)*x(1)+c(2)*x(2))^2+beta*x(2);
   x(1);];
end