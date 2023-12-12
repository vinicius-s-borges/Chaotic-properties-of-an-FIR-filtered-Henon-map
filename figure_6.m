clear all
clc

%Parametros simulação
Numero_bifurc = 1000; 
N_itera = 1500; 
descartN_itera = 1000; 
inicio_bifurc = 0.0001; 
fim_bifurc = 1.8; 
bifurc = linspace(inicio_bifurc,fim_bifurc,Numero_bifurc);

%Parametros mapa
ganho=inicio_bifurc;
a=1.4;
b=0.3;

%%

%Vetor x1 e x2
orbitasx_1 = zeros(N_itera,Numero_bifurc);
orbitasx_2 = zeros(N_itera,Numero_bifurc);

%Condição inicial - Ponto fixo
p1A = (-(1-b)+sqrt((1-b)^2+4*a*(ganho^2)))/(2*ganho^2);
p2A = p1A;
p3A = ganho*p1A;

orbitasx_1(1,1) = p1A;
orbitasx_2(1,1) = p2A;

for inda = 1:Numero_bifurc,
    %Calculo dos coeficientes
    ganho = bifurc(inda);
    Ns=0;
    c = [1 1];
    c = ([1 1]/sum(c))*ganho;
    
    p1A = (-(1-b)+sqrt((1-b)^2+4*a*(ganho^2)))/(2*ganho^2);
    p2A = p1A;
    p3A = ganho*p1A;
      
    if inda>1,
        if isnan(orbitasx_1(n+1,inda-1))==1
        orbitasx_1(1,inda)=p1A;
        orbitasx_2(1,inda)=p2A;
        else
        orbitasx_1(1,inda)=orbitasx_1(n+1,inda-1);
        orbitasx_2(1,inda)=orbitasx_2(n+1,inda-1);
        end
    end
       
    %Vetor x
    x=[orbitasx_1(1,inda);orbitasx_2(1,inda)];
    for n = 1:N_itera-1,
    x(:,n+1) = Henon_N_2(x(:,n),a,b,c);
    end
    orbitasx_1(:,inda)= x(1,:);
    orbitasx_2(:,inda)= x(2,:);
end

bifurc1 = bifurc;
orbitasx_11 = orbitasx_1;

%%

%Parametros simulação
Numero_bifurc = 1000; 
N_itera = 1500; 
descartN_itera = 1000; 
inicio_bifurc = 1.35; 
fim_bifurc = 1.55; 
bifurc = linspace(inicio_bifurc,fim_bifurc,Numero_bifurc);

%Parametros mapa
ganho=inicio_bifurc;
a=1.4;
b=0.3;

%%

%Vetor x1 e x2
orbitasx_1 = zeros(N_itera,Numero_bifurc);
orbitasx_2 = zeros(N_itera,Numero_bifurc);

%Condição inicial - Ponto fixo
p1A = (-(1-b)+sqrt((1-b)^2+4*a*(ganho^2)))/(2*ganho^2);
p2A = p1A;
p3A = ganho*p1A;

orbitasx_1(1,1) = p1A;
orbitasx_2(1,1) = p2A;

for inda = 1:Numero_bifurc,
    %Calculo dos coeficientes
    ganho = bifurc(inda);
    Ns=0;
    c = [1 1];
    c = ([1 1]/sum(c))*ganho;
    
    p1A = (-(1-b)+sqrt((1-b)^2+4*a*(ganho^2)))/(2*ganho^2);
    p2A = p1A;
    p3A = ganho*p1A;
      
    if inda>1,
        if isnan(orbitasx_1(n+1,inda-1))==1
        orbitasx_1(1,inda)=p1A;
        orbitasx_2(1,inda)=p2A;
        else
        orbitasx_1(1,inda)=orbitasx_1(n+1,inda-1);
        orbitasx_2(1,inda)=orbitasx_2(n+1,inda-1);
        end
    end
    
    %Vetor x
    x=[orbitasx_1(1,inda);orbitasx_2(1,inda)];
    for n = 1:N_itera-1,
    x(:,n+1) = Henon_N_2(x(:,n),a,b,c);
    end
    orbitasx_1(:,inda)= x(1,:);
    orbitasx_2(:,inda)= x(2,:);
end

%%
figure
subplot(2,1,1)
plot(bifurc1,orbitasx_11(descartN_itera+1:end,:)','k.', 'MarkerSize',1)
ylabel('$$x_{1}$$','Interpreter','Latex','FontSize',20)
ylim([-2.5 2])
xlim([0 1.8])
grid on
hold on 
x1=inicio_bifurc;
x2=fim_bifurc;
y1=-2.5;
y2=2;
x = [x1, x2, x2, x1, x1];
y = [y1, y1, y2, y2, y1];
plot(x, y, 'r-', 'LineWidth', 2);
hold off
set(gca,'FontSize',20,'LineWidth',2)
subplot(2,1,2)
plot(bifurc,orbitasx_1(descartN_itera+1:end,:)','k.', 'MarkerSize',1)
xlabel('$$G$$','Interpreter','Latex','FontSize',20)
ylabel('$$x_{1}$$','Interpreter','Latex','FontSize',20)
ylim([-2.5 2])
xlim([inicio_bifurc fim_bifurc])
grid on
set(gca,'FontSize',20,'LineWidth',2)


function [x] = Henon_N_2(x,alpha,beta,c)
x=[alpha-(c(1)*x(1)+c(2)*x(2))^2+beta*x(2);
   x(1);];
end


