clear all
clc
%%

Numeroasbifurc = 2000; 
descartebifurc = 1000; 
Nitera = 1500; 
iniciobifurc = 0; 
fimbifurc = 1.2; 
cond_inicial=0.8;
ganho=0.5;

a=1.4;
b=0.3;
asbifurc = linspace(iniciobifurc,fimbifurc,Numeroasbifurc);
%%

p1A = (-(1-b)+sqrt((1-b)^2+4*a*(ganho^2)))/(2*ganho^2);
p2A = p1A;
p3A = ganho*p1A;

orbitas1_1 = zeros(Nitera,Numeroasbifurc);
orbitas1_1(1,1) = cond_inicial;
orbitas2_1 = zeros(Nitera,Numeroasbifurc);
orbitas2_1(1,1) = cond_inicial;
orbitas3_1 = zeros(Nitera,Numeroasbifurc);
orbitas3_1(1,1) = cond_inicial;

for inda = 1:Numeroasbifurc,
    ganho = asbifurc(inda);
    ws = pi/2;
    c = poly([exp(1j*ws) exp(-1j*ws)]);
    c = (poly([exp(1j*ws) exp(-1j*ws)])/sum(c))*ganho;
    
    p1A = (-(1-b)+sqrt((1-b)^2+4*a*(ganho^2)))/(2*ganho^2);
    p2A = p1A;
    p3A = ganho*p1A;
    
    if inda>1,
        if isnan(orbitas1_1(n+1,inda-1))==1
        orbitas1_1(1,inda)=p1A;
        orbitas2_1(1,inda)=p2A;
        orbitas2_1(1,inda)=p3A;
        else
        orbitas1_1(1,inda)=orbitas1_1(n+1,inda-1);
        orbitas2_1(1,inda)=orbitas2_1(n+1,inda-1);
        orbitas3_1(1,inda)=orbitas3_1(n+1,inda-1);
        end
    end    

    x=[orbitas1_1(1,inda);orbitas2_1(1,inda);orbitas3_1(1,inda)];%;orbitas4_1(1,inda)];
    for n = 1:Nitera-1,
    x(:,n+1) = Henon_N_3(x(:,n),a,b,c);
    end
    orbitas1_1(:,inda)= x(1,:);
    orbitas2_1(:,inda)= x(2,:);
    orbitas3_1(:,inda)= x(3,:);
end
%%
figure
plot(asbifurc,orbitas1_1(descartebifurc+1:end,:)','k.', 'MarkerSize',1);
xlabel('$$G$$','Interpreter','Latex','FontSize',18)
ylabel('$$x_{1}$$','Interpreter','Latex','FontSize',18)
ylim([-3.2 2.2])
xlim([iniciobifurc fimbifurc])
grid on
set(gca,'FontSize',24,'LineWidth',2)

function [x] = Henon_N_3(x,alpha,beta,c)
x=[alpha-x(3)^2+beta*x(2);
   x(1);
   c(1)*(alpha-x(3)^2+beta*x(2))+c(2)*x(1)+c(3)*x(2);];
end
