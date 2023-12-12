clear all
clc

%Simulation parameters
Nitera=500; % number of iterations discarded
N_itera_mapa=3000; % number of iterations 
ninicial=25; % number of random variables 


%Map parameters
a=1.4;
b=0.3;

%Filter parameters
coef = 0:1:40; % zeros
Ncoef = length(coef); % number of zeros
Nganhos = 100; % number of gain
ganhos = linspace(0.000001,2.5,Nganhos); % gain

%% Fixed Point

for iteraNs=1:Ncoef,
for iterag=1:Nganhos,  
    ganho = ganhos(iterag);
    c = [1 zeros(1,coef(iteraNs)) 1];
    c = ([1 zeros(1,coef(iteraNs)) 1]/sum(abs(c)))*ganho;
    Ns=length(c);
    gamma = ganho;
    
    if Ns==1
        
        %Fixed Point +
        p1A = (-(1-b)+sqrt((1-b)^2+4*a*(gamma^2)))/(2*gamma^2);
        p2A = p1A;
        p3A = gamma*p1A;
        pontofixoA = [p1A;p2A;p3A;p1A*ones(Ns-3,1)];

        JacobianoA = dhenonNws_1_pf(pontofixoA,Ns,c);

        autovaloresA = abs(eig(JacobianoA));
        
        autovA(iteraNs,iterag) = max(autovaloresA);
        
        %Fixed Point -
        p1B = (-(1-b)-sqrt((1-b)^2+4*a*(gamma^2)))/(2*gamma^2);
        p2B = p1B;
        p3B = gamma*p1B;
        pontofixoB = [p1B;p2B;p3B;p1B*ones(Ns-3,1)];

        JacobianoB = dhenonNws_1_pf(pontofixoB,Ns,c);

        autovaloresB = abs(eig(JacobianoB));    
        
        autovB(iteraNs,iterag) = max(autovaloresB);
        
    end

    if Ns==2
        
        %Fixed Point +
        p1A = (-(1-b)+sqrt((1-b)^2+4*a*(gamma^2)))/(2*gamma^2);
        p2A = p1A;
        p3A = gamma*p1A;
        pontofixoA = [p1A;p2A;p3A;p1A*ones(Ns-3,1)];

        JacobianoA = dhenonNws_2_pf(pontofixoA,Ns,c);

        autovaloresA = abs(eig(JacobianoA));
        
        autovA(iteraNs,iterag) = max(autovaloresA);
        
        %Fixed Point -
        p1B = (-(1-b)-sqrt((1-b)^2+4*a*(gamma^2)))/(2*gamma^2);
        p2B = p1B;
        p3B = gamma*p1B;
        pontofixoB = [p1B;p2B;p3B;p1B*ones(Ns-3,1)];

        JacobianoB = dhenonNws_2_pf(pontofixoB,Ns,c);

        autovaloresB = abs(eig(JacobianoB));    
        
        autovB(iteraNs,iterag) = max(autovaloresB);
        
    end
    
    if Ns==3
        
        %Fixed Point +
        p1A = (-(1-b)+sqrt((1-b)^2+4*a*(gamma^2)))/(2*gamma^2);
        p2A = p1A;
        p3A = gamma*p1A;
        pontofixoA = [p1A;p2A;p3A;p1A*ones(Ns-3,1)];

        JacobianoA = dhenonNws_3_pf(pontofixoA,Ns,c);

        autovaloresA = abs(eig(JacobianoA));
        
        autovA(iteraNs,iterag) = max(autovaloresA);
        
        %Fixed Point -
        p1B = (-(1-b)-sqrt((1-b)^2+4*a*(gamma^2)))/(2*gamma^2);
        p2B = p1B;
        p3B = gamma*p1B;
        pontofixoB = [p1B;p2B;p3B;p1B*ones(Ns-3,1)];

        JacobianoB = dhenonNws_3_pf(pontofixoB,Ns,c);

        autovaloresB = abs(eig(JacobianoB));    
        
        autovB(iteraNs,iterag) = max(autovaloresB);
        
    end    
    
    if length(c)>3
        %Fixed Point +
        p1A = (-(1-b)+sqrt((1-b)^2+4*a*(gamma^2)))/(2*gamma^2);
        p2A = p1A;
        p3A = gamma*p1A;
        pontofixoA = [p1A;p2A;p3A;p1A*ones(Ns-3,1)];

        JacobianoA = dhenonNws_4_pf(pontofixoA,Ns,c);

        autovaloresA = abs(eig(JacobianoA));
        
        autovA(iteraNs,iterag) = max(autovaloresA);
        
        %Fixed Point -
        p1B = (-(1-b)-sqrt((1-b)^2+4*a*(gamma^2)))/(2*gamma^2);
        p2B = p1B;
        p3B = gamma*p1B;
        pontofixoB = [p1B;p2B;p3B;p1B*ones(Ns-3,1)];

        JacobianoB = dhenonNws_4_pf(pontofixoB,Ns,c);

        autovaloresB = abs(eig(JacobianoB));    
        
        autovB(iteraNs,iterag) = max(autovaloresB);
    end
    
    autovAA(iteraNs,iterag) = autovA(iteraNs,iterag);
    if autovA(iteraNs,iterag)>1
        autovAp(iteraNs,iterag) = NaN;
    else
        autovAp(iteraNs,iterag) = 2;
    end
    autovBB(iteraNs,iterag) = autovB(iteraNs,iterag);
    if autovB(iteraNs,iterag)>1
        autovBp(iteraNs,iterag) = NaN;
    else
        autovBp(iteraNs,iterag) = 2;
    end    
end
end

%%

figure
pcolor(ganhos,coef,autovAp)
title('Ponto fixo p1')
xlabel({'$ganho$'},'Interpreter','latex')
ylabel({'$N_{S}$'},'Interpreter','latex')
shading flat
grid on
set(gca,'layer','top','FontSize',24)

%% Lyapunov
for iteraNs=1:Ncoef,
    iteraNs
    for iterag=1:Nganhos,  
        ganho = ganhos(iterag);
        
        if coef(iteraNs)==0        
        c = [1 1];
        c = ([1 1]/sum(c))*ganho;
        end
        
        if coef(iteraNs) > 0 
        c = [1 zeros(1,coef(iteraNs)) 1];
        c = ([1 zeros(1,coef(iteraNs)) 1]/sum(c))*ganho;
        end
        
        Ns=length(c);
    
        for condinicial=1:ninicial,
            
        p1A = (-(1-b)+sqrt((1-b)^2+4*a*(ganho^2)))/(2*ganho^2);
        p2A = p1A;
        p3A = ganho*p1A;
        pontofixoA = [p1A;p2A;p3A;p1A*ones(Ns-3,1)];
        
            if Ns==1
                x=0.01*rand(2,1)+p1A;
                
                for n=1:1:Nitera
                    [x] = Henon_N_1(x,a,b,c);
                end
            
            if isnan(x)==1
                h_vetor(condinicial)=NaN;
                break;
            end
        
            [lyapunov] = lyapunov_N_1(x,N_itera_mapa,Ns,c,a,b);
    
            h_vetor(condinicial)=lyapunov;
            end
        
            if Ns==2
                x=0.01*rand(Ns,1)+[p1A; p2A];
                
                for n=1:1:Nitera
                    [x] = Henon_N_2(x,a,b,c);
                end
        
            if isnan(x)==1
                h_vetor(condinicial)=NaN;
                break;
            end
        
            [lyapunov] = lyapunov_N_2(x,N_itera_mapa,Ns,c,a,b);
    
            h_vetor(condinicial)=lyapunov;
            end
        
            if Ns==3
                x=0.01*rand(Ns,1)+pontofixoA;
                
                for n=1:1:Nitera
                    [x] = Henon_N_3(x,a,b,c);
                end
        
            if isnan(x)==1
                h_vetor(condinicial)=NaN;
                break;
            end
        
            [lyapunov] = lyapunov_N_3(x,N_itera_mapa,Ns,c,a,b);
    
            h_vetor(condinicial)=lyapunov;
            end
        
            if Ns>3
                x=0.01*rand(Ns,1)+pontofixoA;
                
                for n=1:1:Nitera
                    [x] = Henon_N_4(x,a,b,c);
                end

            if isnan(x)==1
                h_vetor(condinicial)=NaN;
                break;
            end
        
            [lyapunov] = lyapunov_N_4(x,N_itera_mapa,Ns,c,a,b);
    
            h_vetor(condinicial)=lyapunov;
            end        
        
        end
    h(iteraNs,iterag)=mean(h_vetor);
    h_desvio(iteraNs,iterag)=std(h_vetor);
    end

iteraNs
end

save('variables_lyapunov_teste_1_sum')

%%
figure
pcolor(ganhos,coef,h)
xlabel({'$\frac{\omega_{s}}{\pi}$'},'Interpreter','latex')
ylabel({'$N_{Z}$'},'Interpreter','latex')
shading flat
colorbar 
grid on
set(gca,'layer','top','FontSize',24)

%%

for iteraNs=1:Ncoef,
    for iterag=1:Nganhos,
        if h(iteraNs,iterag)<0;
            h_color(iteraNs,iterag)=-1;
        end
        if h(iteraNs,iterag)>0
            h_color(iteraNs,iterag)=0;
        end
        if isnan(h(iteraNs,iterag))==1
            h_color(iteraNs,iterag)=1;
        end             
    end
end

%%

for iteraNs=1:Ncoef,
    for iterag=1:Nganhos,
        if h(iteraNs,iterag)<0 & autovAp(iteraNs,iterag)==2;
            h_color_1(iteraNs,iterag)=-2;
            h_color_2(iteraNs,iterag)=NaN;
            h_color_3(iteraNs,iterag)=NaN;
            h_color_4(iteraNs,iterag)=NaN;
        end
        if h(iteraNs,iterag)<0 & autovAp(iteraNs,iterag)~=2;
            h_color_1(iteraNs,iterag)=NaN;
            h_color_2(iteraNs,iterag)=-1;
            h_color_3(iteraNs,iterag)=NaN;
            h_color_4(iteraNs,iterag)=NaN;
        end        
        if h(iteraNs,iterag)>0
            h_color_1(iteraNs,iterag)=NaN;
            h_color_2(iteraNs,iterag)=NaN;           
            h_color_3(iteraNs,iterag)=0;
            h_color_4(iteraNs,iterag)=NaN;
        end
        if isnan(h(iteraNs,iterag))==1
            h_color_1(iteraNs,iterag)=NaN;
            h_color_2(iteraNs,iterag)=NaN;           
            h_color_3(iteraNs,iterag)=NaN;
            h_color_4(iteraNs,iterag)=1;
        end             
    end
end

%%
figure
pcolor(ganhos,coef,h_color)
xlabel({'$\gamma$'},'Interpreter','latex')
ylabel({'$N_{Z}$'},'Interpreter','latex')
shading flat
grid on
hold on
pcolor(ganhos,coef,autovAp)
shading flat
grid on
hold off
hold on
plot(0.3,10,'x') 
set(gca,'layer','top','FontSize',24)
plot(0.5,20,'x') 
set(gca,'layer','top','FontSize',24)
plot(0.6,30,'x') 
set(gca,'layer','top','FontSize',24)
hold off
set(gca,'layer','top','FontSize',24)

figure
pcolor(coef,ganhos,h_color')
ylabel({'$G$'},'Interpreter','latex')
xlabel({'$N_{Z}$'},'Interpreter','latex')
shading flat
grid on
hold on
pcolor(coef,ganhos,autovAp')
legend;
shading flat
grid on
hold off
hold on
plot(0.5,41.5,'x') 
set(gca,'layer','top','FontSize',24)
plot(0.8,25.5,'x') 
set(gca,'layer','top','FontSize',24)
plot(1.25,25.5,'x') 
set(gca,'layer','top','FontSize',24)
hold off
set(gca,'layer','top','FontSize',24)
hold on
xline(0,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(1,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(2,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(3,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(4,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(5,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(6,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(7,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(8,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(9,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(10,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(11,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(12,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(13,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(14,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(15,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(16,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(17,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(18,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(19,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(20,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(21,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(22,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(23,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(24,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(25,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(26,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(27,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(28,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(29,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(30,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(31,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(32,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(33,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(34,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(35,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(36,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(37,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(38,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
xline(39,'-k','LineWidth',2)
set(gca,'layer','top','FontSize',24)
hold off
grid on
ax = gca;
ax.XTick = 0.5:1:40;
ax.XTickLabels = 1:1:40;
set(gca,'layer','top','FontSize',24)

%%

ncoef_int=1:2*Ncoef;
k=1;
for iteraNs=1:2*Ncoef,
    if mod(ncoef_int(iteraNs),2)==1;
        h_color_int(iteraNs,:)=h_color(k,:);
        autovAp_int(iteraNs,:)=autovAp(k,:);
        k=k+1;
    else
        h_color_int(iteraNs,:)=NaN*ones(1,100);
        autovAp_int(iteraNs,:)=NaN*ones(1,100);
    end        
end

%%
ncoef_int=1:2*Ncoef;
k=1;

for iteraNs=1:2*Ncoef,
        if mod(ncoef_int(iteraNs),2)==1;
            h_color_1_int(iteraNs,:)=h_color_1(k,:);
            h_color_2_int(iteraNs,:)=h_color_2(k,:);
            h_color_3_int(iteraNs,:)=h_color_3(k,:);
            h_color_4_int(iteraNs,:)=h_color_4(k,:);
            k=k+1;
        else
            h_color_1_int(iteraNs,:)=NaN*ones(1,100);
            h_color_2_int(iteraNs,:)=NaN*ones(1,100);
            h_color_3_int(iteraNs,:)=NaN*ones(1,100);
            h_color_4_int(iteraNs,:)=NaN*ones(1,100);
        end
end

%%
figure
pcolor(ncoef_int(1,1:80),ganhos,h_color_1_int(1:80,:)')
ylabel({'$G$'},'Interpreter','latex')
xlabel({'$N_{Z}$'},'Interpreter','latex')
shading flat
hold on
pcolor(ncoef_int(1,1:80),ganhos,h_color_2_int(1:80,:)')
pcolor(ncoef_int(1,1:80),ganhos,h_color_3_int(1:80,:)')
pcolor(ncoef_int(1,1:80),ganhos,h_color_4_int(1:80,:)')
shading flat
plot(41.5,0.5,'x') 
set(gca,'layer','top','FontSize',24)
plot(25.5,0.8,'x') 
set(gca,'layer','top','FontSize',24)
plot(25.5,1.25,'x') 
set(gca,'layer','top','FontSize',24)
xline(1,'-k','LineWidth',2)
xline(2,'-k','LineWidth',2)
xline(3,'-k','LineWidth',2)
xline(4,'-k','LineWidth',2)
xline(5,'-k','LineWidth',2)
xline(6,'-k','LineWidth',2)
xline(7,'-k','LineWidth',2)
xline(8,'-k','LineWidth',2)
xline(9,'-k','LineWidth',2)
xline(10,'-k','LineWidth',2)
xline(11,'-k','LineWidth',2)
xline(12,'-k','LineWidth',2)
xline(13,'-k','LineWidth',2)
xline(14,'-k','LineWidth',2)
xline(15,'-k','LineWidth',2)
xline(16,'-k','LineWidth',2)
xline(17,'-k','LineWidth',2)
xline(18,'-k','LineWidth',2)
xline(19,'-k','LineWidth',2)
xline(20,'-k','LineWidth',2)
xline(21,'-k','LineWidth',2)
xline(22,'-k','LineWidth',2)
xline(23,'-k','LineWidth',2)
xline(24,'-k','LineWidth',2)
xline(25,'-k','LineWidth',2)
xline(26,'-k','LineWidth',2)
xline(27,'-k','LineWidth',2)
xline(28,'-k','LineWidth',2)
xline(29,'-k','LineWidth',2)
xline(30,'-k','LineWidth',2)
xline(31,'-k','LineWidth',2)
xline(32,'-k','LineWidth',2)
xline(33,'-k','LineWidth',2)
xline(34,'-k','LineWidth',2)
xline(35,'-k','LineWidth',2)
xline(36,'-k','LineWidth',2)
xline(37,'-k','LineWidth',2)
xline(38,'-k','LineWidth',2)
xline(39,'-k','LineWidth',2)
xline(40,'-k','LineWidth',2)
xline(41,'-k','LineWidth',2)
xline(42,'-k','LineWidth',2)
xline(43,'-k','LineWidth',2)
xline(44,'-k','LineWidth',2)
xline(45,'-k','LineWidth',2)
xline(46,'-k','LineWidth',2)
xline(47,'-k','LineWidth',2)
xline(48,'-k','LineWidth',2)
xline(49,'-k','LineWidth',2)
xline(50,'-k','LineWidth',2)
xline(51,'-k','LineWidth',2)
xline(52,'-k','LineWidth',2)
xline(53,'-k','LineWidth',2)
xline(54,'-k','LineWidth',2)
xline(55,'-k','LineWidth',2)
xline(56,'-k','LineWidth',2)
xline(57,'-k','LineWidth',2)
xline(58,'-k','LineWidth',2)
xline(59,'-k','LineWidth',2)
xline(60,'-k','LineWidth',2)
xline(61,'-k','LineWidth',2)
xline(62,'-k','LineWidth',2)
xline(63,'-k','LineWidth',2)
xline(64,'-k','LineWidth',2)
xline(65,'-k','LineWidth',2)
xline(66,'-k','LineWidth',2)
xline(67,'-k','LineWidth',2)
xline(68,'-k','LineWidth',2)
xline(69,'-k','LineWidth',2)
xline(70,'-k','LineWidth',2)
xline(71,'-k','LineWidth',2)
xline(72,'-k','LineWidth',2)
xline(73,'-k','LineWidth',2)
xline(74,'-k','LineWidth',2)
xline(75,'-k','LineWidth',2)
xline(76,'-k','LineWidth',2)
xline(77,'-k','LineWidth',2)
xline(78,'-k','LineWidth',2)
xline(79,'-k','LineWidth',2)
xline(80,'-k','LineWidth',2)
xline(81,'-k','LineWidth',2)
xline(82,'-k','LineWidth',2)
hold off
grid on
ax = gca;
ax.XTick = [1.5 9.5 19.5 29.5 39.5 49.5 59.5 69.5 79.5];
ax.XTickLabels = [1 5 10 15 20 25 30 35 40];
legend
set(gca,'layer','top','FontSize',18)

%%
coef_impar = 1:2:39;
coef_par = 0:2:40;
k=1;
n=1;
for iteraNs=1:Ncoef,
    if mod(coef(iteraNs),2)==1;
        h_color_impar(k,:)=h_color(iteraNs,:);
        autovAp_impar(k,:)=autovAp(iteraNs,:);
        k=k+1;
    else
        h_color_par(n,:)=h_color(iteraNs,:);
        autovAp_par(n,:)=autovAp(iteraNs,:);
        n=n+1;
    end        
end

%%
for iteraNs=1:Ncoef,
    if mod(coef(iteraNs),2)==1;
        h_color_impar_2(iteraNs,:)=h_color(iteraNs,:);
        h_color_par_2(iteraNs,:)=NaN*ones(1,100);
        autovAp_impar_2(iteraNs,:)=autovAp(iteraNs,:);
        autovAp_par_2(iteraNs,:)=NaN*ones(1,100);
    else
        h_color_par_2(iteraNs,:)=h_color(iteraNs,:);
        h_color_impar_2(iteraNs,:)=NaN*ones(1,100);
        autovAp_par_2(iteraNs,:)=autovAp(iteraNs,:);
        autovAp_impar_2(iteraNs,:)=NaN*ones(1,100);
    end        
end

%%

figure
pcolor(coef,ganhos,h_color_impar_2')
ylabel({'$G$'},'Interpreter','latex')
xlabel({'$N_{Z}$'},'Interpreter','latex')
shading flat
hold on
pcolor(coef,ganhos,autovAp_impar_2')
shading flat
grid on
hold off
hold on 
plot(5.5,0.45,'x') 
set(gca,'layer','top','FontSize',24)
xline(1,'-k','LineWidth',2)
xline(2,'-k','LineWidth',2)
xline(3,'-k','LineWidth',2)
xline(4,'-k','LineWidth',2)
xline(5,'-k','LineWidth',2)
xline(6,'-k','LineWidth',2)
xline(7,'-k','LineWidth',2)
xline(8,'-k','LineWidth',2)
xline(9,'-k','LineWidth',2)
xline(10,'-k','LineWidth',2)
xline(11,'-k','LineWidth',2)
xline(12,'-k','LineWidth',2)
xline(13,'-k','LineWidth',2)
xline(14,'-k','LineWidth',2)
xline(15,'-k','LineWidth',2)
xline(16,'-k','LineWidth',2)
xline(17,'-k','LineWidth',2)
xline(18,'-k','LineWidth',2)
xline(19,'-k','LineWidth',2)
xline(20,'-k','LineWidth',2)
xline(21,'-k','LineWidth',2)
xline(22,'-k','LineWidth',2)
xline(23,'-k','LineWidth',2)
xline(24,'-k','LineWidth',2)
xline(25,'-k','LineWidth',2)
xline(26,'-k','LineWidth',2)
xline(27,'-k','LineWidth',2)
xline(28,'-k','LineWidth',2)
xline(29,'-k','LineWidth',2)
xline(30,'-k','LineWidth',2)
xline(31,'-k','LineWidth',2)
xline(32,'-k','LineWidth',2)
xline(33,'-k','LineWidth',2)
xline(34,'-k','LineWidth',2)
xline(35,'-k','LineWidth',2)
xline(36,'-k','LineWidth',2)
xline(37,'-k','LineWidth',2)
xline(38,'-k','LineWidth',2)
xline(39,'-k','LineWidth',2)
xline(40,'-k','LineWidth',2)
hold off
grid on
ax = gca;
ax.XTick = 1.5:2:40;
ax.XTickLabels = 1:2:40;
set(gca,'layer','top','FontSize',24)

figure
pcolor(coef,ganhos,h_color_par_2')
ylabel({'$G$'},'Interpreter','latex')
xlabel({'$N_{Z}$'},'Interpreter','latex')
shading flat
hold on
pcolor(coef,ganhos,autovAp_par_2')
shading flat
grid on
hold off
hold on 
plot(10.5,0.2,'x') 
set(gca,'layer','top','FontSize',24)
plot(10.5,0.6,'x') 
set(gca,'layer','top','FontSize',24)
xline(0,'-k','LineWidth',2)
xline(1,'-k','LineWidth',2)
xline(2,'-k','LineWidth',2)
xline(3,'-k','LineWidth',2)
xline(4,'-k','LineWidth',2)
xline(5,'-k','LineWidth',2)
xline(6,'-k','LineWidth',2)
xline(7,'-k','LineWidth',2)
xline(8,'-k','LineWidth',2)
xline(9,'-k','LineWidth',2)
xline(10,'-k','LineWidth',2)
xline(11,'-k','LineWidth',2)
xline(12,'-k','LineWidth',2)
xline(13,'-k','LineWidth',2)
xline(14,'-k','LineWidth',2)
xline(15,'-k','LineWidth',2)
xline(16,'-k','LineWidth',2)
xline(17,'-k','LineWidth',2)
xline(18,'-k','LineWidth',2)
xline(19,'-k','LineWidth',2)
xline(20,'-k','LineWidth',2)
xline(21,'-k','LineWidth',2)
xline(22,'-k','LineWidth',2)
xline(23,'-k','LineWidth',2)
xline(24,'-k','LineWidth',2)
xline(25,'-k','LineWidth',2)
xline(26,'-k','LineWidth',2)
xline(27,'-k','LineWidth',2)
xline(28,'-k','LineWidth',2)
xline(29,'-k','LineWidth',2)
xline(30,'-k','LineWidth',2)
xline(31,'-k','LineWidth',2)
xline(32,'-k','LineWidth',2)
xline(33,'-k','LineWidth',2)
xline(34,'-k','LineWidth',2)
xline(35,'-k','LineWidth',2)
xline(36,'-k','LineWidth',2)
xline(37,'-k','LineWidth',2)
xline(38,'-k','LineWidth',2)
xline(39,'-k','LineWidth',2)
hold off
grid on
legend
ax = gca;
ax.XTick = 0.5:2:40;
ax.XTickLabels = 0:2:40;
set(gca,'layer','top','FontSize',24)

%%

function [x] = dhenonNws_1_pf(x,d,c)
b=0.3;

x=[-2*c(1)^2*x(1)       b             
    1                    0;];
end

function [x] = dhenonNws_2_pf(x,d,c)
b=0.3;

x=[-((2*c(1)^2*x(1))+(2*c(1)*c(2)*x(2)))       -((2*c(1)*x(1))+(2*c(1)*c(2)^2*x(2)))+b             
    1                                           0;];
end

function [x] = dhenonNws_3_pf(x,d,c)
b=0.3;

x=[0        b              -2*x(3)    
   1        0              0
   c(2)     c(1)*b+c(3)    -2*c(1)*x(3);];
end

function [dx] = dhenonNws_4_pf(x,d,c)
b=0.3;

x1=[0       b              -2*x(3);
    1       0              0;
    c(2)    c(1)*b+c(3)    -2*c(1)*x(3);
    0       1              0;];

x2=zeros(2,d-3);

x3=zeros(d-4,3);

x4=c(4:end);

x5=zeros(1,d-3);

x6=eye(d-4,d-4);

x7=zeros(d-4,1);

dx1=vertcat(x1,x3);

dx2=vertcat(x2,x4,x5);

dx3=horzcat(x6,x7);

dx4=vertcat(dx2,dx3);

dx=horzcat(dx1,dx4);
end

function [lyapunov] = lyapunov_N_1(x,Nitera,Ncoeficientes,coeficientes,alfa,beta)

w = eye(Ncoeficientes); %Initial orthonormal base

    for i = 1:Nitera, 
        dx = dhenon_N_1(x,beta,coeficientes); 
        z = dx*w; 
         
        %Gram-Schmidt orthogonalization 
        y = gsog(z); 
        w=[]; 
         
        %Standardization 
        for k=1:Ncoeficientes 
            r(k,i)=norm(y(:,k)); 
            w = [w y(:,k)/norm(y(:,k))]; 
        end 
             
        %Map iteration 
        x=Henon_N_1(x,alfa,beta,coeficientes); 
    end 
    
    %Calculation of Lyapunov exponents
    for k=1:Ncoeficientes
    lyapunov(k) = sum(log(r(k,:)))/Nitera;
    end
    
    lyapunov=max(lyapunov);
    
end

function [x] = Henon_N_1(x,alpha,beta,c)
x=[alpha-(c(1)*x(1))^2+beta*x(2);
   x(1);];
end

function [dx] = dhenon_N_1(x,beta,c)
dx=[-2*c(1)^2*x(1) beta;
    1              0;];
end

function [lyapunov] = lyapunov_N_2(x,Nitera,Ncoeficientes,coeficientes,alfa,beta)

w = eye(Ncoeficientes); %Initial orthonormal basis

    for i = 1:Nitera, 
        dx = dhenon_N_2(x,beta,coeficientes); 
        z = dx*w; 
         
        %Gram-Schmidt orthogonalization 
        y = gsog(z); 
        w=[]; 
         
        %Standardization 
        for k=1:Ncoeficientes 
            r(k,i)=norm(y(:,k)); 
            w = [w y(:,k)/norm(y(:,k))]; 
        end 
             
        %Map iteration 
        x=Henon_N_2(x,alfa,beta,coeficientes); 
    end 
    
    %Calculation of Lyapunov exponents
    for k=1:Ncoeficientes
    lyapunov(k) = sum(log(r(k,:)))/Nitera;
    end
    
    lyapunov=max(lyapunov);
    
end

function [x] = Henon_N_2(x,alpha,beta,c)
x=[alpha-(c(1)*x(1)+c(2)*x(2))^2+beta*x(2);
   x(1);];
end

function [dx] = dhenon_N_2(x,b,c)
dx=[-2*c(1)^2*x(1)-2*c(1)*c(2)*x(2) -2*c(1)*c(2)*x(1)-2*c(2)^2*x(2)+b;
    1                         0;];
end

function [lyapunov] = lyapunov_N_3(x,Nitera,Ncoeficientes,coeficientes,alfa,beta)

w = eye(Ncoeficientes); %Initial orthonormal base

    for i = 1:Nitera, 
        dx = dhenon_N_3(x,beta,coeficientes); 
        z = dx*w; 
         
        %Gram-Schmidt orthogonalization 
        y = gsog(z); 
        w=[]; 
         
        %Normalization 
        for k=1:Ncoeficientes 
            r(k,i)=norm(y(:,k)); 
            w = [w y(:,k)/norm(y(:,k))]; 
        end 
             
        %Map iteration 
        x=Henon_N_3(x,alfa,beta,coeficientes); 
    end 
    
    %Calculation of Lyapunov exponents
    for k=1:Ncoeficientes
    lyapunov(k) = sum(log(r(k,:)))/Nitera;
    end
    
    lyapunov=max(lyapunov);
    
end

function [x] = Henon_N_3(x,alpha,beta,c)
x=[alpha-x(3)^2+beta*x(2);
   x(1);
   c(1)*(alpha-x(3)^2+beta*x(2))+c(2)*x(1)+c(3)*x(2);];
end

function [dx] = dhenon_N_3(x,beta,c)
dx=[0       beta           -2*x(3)     ;
    1       0              0           ;
    c(2)    (c(1)*beta)+c(3) -2*c(1)*x(3);];
end

function [lyapunov] = lyapunov_N_4(x,Nitera,Ncoeficientes,coeficientes,alfa,beta)

w = eye(Ncoeficientes); %Initial orthonormal basis

    for i = 1:Nitera
        dx = dhenon_N_4(x,Ncoeficientes,coeficientes); 
        z = dx*w; 
         
        %Gram-Schmidt orthogonalization 
        y = gsog(z); 
        w=[]; 
         
        %Normalization 
        for k=1:Ncoeficientes 
            r(k,i)=norm(y(:,k)); 
            w = [w y(:,k)/norm(y(:,k))]; 
        end 
             
        %Map iteration 
        x=Henon_N_4(x,alfa,beta,coeficientes); 
    end 
    
    %Calculation of Lyapunov exponents
    for k=1:Ncoeficientes
    lyapunov(k) = sum(log(r(k,:)))/Nitera;
    end
    
    lyapunov=max(lyapunov);
end
    
function [Q, R] = gsog(X)
% Gram-Schmidt orthogonalization
% Written by Mo Chen (sth4nth@gmail.com).
[d,n] = size(X);
m = min(d,n);
R = eye(m,n);
Q = zeros(d,m);
D = zeros(1,m);
for i = 1:m
    R(1:i-1,i) = bsxfun(@times,Q(:,1:i-1),1./D(1:i-1))'*X(:,i);
    Q(:,i) = X(:,i)-Q(:,1:i-1)*R(1:i-1,i);
    D(i) = dot(Q(:,i),Q(:,i));
end
R(:,m+1:n) = bsxfun(@times,Q,1./D)'*X(:,m+1:n);
end

function [x] = henonNws_2(x,a,b,c)

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

function [dx] = dhenonNws_2(x,d,c)
b=0.3;

x1=[0       b              -2*x(3);
    1       0              0;
    c(2)    c(1)*b+c(3)    -2*c(1)*x(3);
    0       1              0;];

x2=zeros(2,d-3);
x3=zeros(d-4,3);
x4=c(4:end);
x5=zeros(1,d-3);
x6=eye(d-4,d-4);
x7=zeros(d-4,1);
dx1=vertcat(x1,x3);
dx2=vertcat(x2,x4,x5);
dx3=horzcat(x6,x7);
dx4=vertcat(dx2,dx3);
dx=horzcat(dx1,dx4);
end

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

function [dx] = dhenon_N_4(x,d,c)
b=0.3;

x1=[0       b              -2*x(3);
    1       0              0;
    c(2)    c(1)*b+c(3)    -2*c(1)*x(3);
    0       1              0;];

x2=zeros(2,d-3);

x3=zeros(d-4,3);

x4=c(4:end);

x5=zeros(1,d-3);

x6=eye(d-4,d-4);

x7=zeros(d-4,1);

dx1=vertcat(x1,x3);

dx2=vertcat(x2,x4,x5);

dx3=horzcat(x6,x7);

dx4=vertcat(dx2,dx3);

dx=horzcat(dx1,dx4);
end