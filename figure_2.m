clc
clear all

%Map parameters
a=1.4;
b=0.3;

%Filter parameters
Nganhos = 1005
G = linspace(-10,10,Nganhos);

%Fixed point
for iterag=1:Nganhos, 
    p1mais(iterag) = (-(1-b)+sqrt((1-b)^2+4*a*(G(iterag)^2)))/(2*G(iterag)^2);
    p1menos(iterag) = (-(1-b)-sqrt((1-b)^2+4*a*(G(iterag)^2)))/(2*G(iterag)^2);
end

figure
plot(G,p1mais,'--','Color','black','LineWidth',1)
hold on
plot(G,p1menos,'Color','black','LineWidth',1)
hold off
ylabel('$p^{+;-}_{1}$','Interpreter','latex','FontSize', 1)
xlabel('$G$','Interpreter','latex','FontSize', 1)
ylim([-10 10])
ylim([-5 2])
grid on
legend('$p^{+}_{1}$','$p^{-}_{1}$','Interpreter','latex')
set(gca,'FontSize',20,'LineWidth',2)