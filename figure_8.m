clear all
clc

%%

%Simulation parameters
Nitera=500; % number of iterations discarded
N_itera_mapa=3000; % number of iterations 
ninicial=25; % number of random variables 


%Map parameters
a=1.4;
b=0.3;

%Filter parameters
Nws = 500; % number of frequencies
w=linspace(0.00001,3.14,Nws); % frequencies
Nganhos = 500; % number of gain
ganhos = linspace(0.0001,3,Nganhos); % gain
%%

for iteraws=1:Nws,
for iterag=1:Nganhos,  
    ganho = ganhos(iterag);
    ws = w(iteraws);
    c = poly([exp(1j*ws) exp(-1j*ws)]);
    c = (poly([exp(1j*ws) exp(-1j*ws)])/sum(c))*ganho;
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
        
        autovA(iteraws,iterag) = max(autovaloresA);
        
        %Fixed Point -
        p1B = (-(1-b)-sqrt((1-b)^2+4*a*(gamma^2)))/(2*gamma^2);
        p2B = p1B;
        p3B = gamma*p1B;
        pontofixoB = [p1B;p2B;p3B;p1B*ones(Ns-3,1)];

        JacobianoB = dhenonNws_1_pf(pontofixoB,Ns,c);

        autovaloresB = abs(eig(JacobianoB));    
        
        autovB(iteraws,iterag) = max(autovaloresB);
        
    end

    if Ns==2
        
        %Fixed Point +
        p1A = (-(1-b)+sqrt((1-b)^2+4*a*(gamma^2)))/(2*gamma^2);
        p2A = p1A;
        p3A = gamma*p1A;
        pontofixoA = [p1A;p2A;p3A;p1A*ones(Ns-3,1)];

        JacobianoA = dhenonNws_2_pf(pontofixoA,Ns,c);

        autovaloresA = abs(eig(JacobianoA));
        
        autovA(iteraws,iterag) = max(autovaloresA);
        
        %Fixed Point -
        p1B = (-(1-b)-sqrt((1-b)^2+4*a*(gamma^2)))/(2*gamma^2);
        p2B = p1B;
        p3B = gamma*p1B;
        pontofixoB = [p1B;p2B;p3B;p1B*ones(Ns-3,1)];

        JacobianoB = dhenonNws_2_pf(pontofixoB,Ns,c);

        autovaloresB = abs(eig(JacobianoB));    
        
        autovB(iteraws,iterag) = max(autovaloresB);
        
    end
    
    if Ns==3
        
        %Fixed Point +
        p1A = (-(1-b)+sqrt((1-b)^2+4*a*(gamma^2)))/(2*gamma^2);
        p2A = p1A;
        p3A = gamma*p1A;
        pontofixoA = [p1A;p2A;p3A;p1A*ones(Ns-3,1)];

        JacobianoA = dhenonNws_3_pf(pontofixoA,Ns,c);

        autovaloresA = abs(eig(JacobianoA));
        
        autovA(iteraws,iterag) = max(autovaloresA);
        
        %Fixed Point -
        p1B = (-(1-b)-sqrt((1-b)^2+4*a*(gamma^2)))/(2*gamma^2);
        p2B = p1B;
        p3B = gamma*p1B;
        pontofixoB = [p1B;p2B;p3B;p1B*ones(Ns-3,1)];

        JacobianoB = dhenonNws_3_pf(pontofixoB,Ns,c);

        autovaloresB = abs(eig(JacobianoB));    
        
        autovB(iteraws,iterag) = max(autovaloresB);
        
    end    
    
    if length(c)>3
        %Fixed Point +
        p1A = (-(1-b)+sqrt((1-b)^2+4*a*(gamma^2)))/(2*gamma^2);
        p2A = p1A;
        p3A = gamma*p1A;
        pontofixoA = [p1A;p2A;p3A;p1A*ones(Ns-3,1)];

        JacobianoA = dhenonNws_4_pf(pontofixoA,Ns,c);

        autovaloresA = abs(eig(JacobianoA));
        
        autovA(iteraws,iterag) = max(autovaloresA);
        
        %Fixed Point -
        p1B = (-(1-b)-sqrt((1-b)^2+4*a*(gamma^2)))/(2*gamma^2);
        p2B = p1B;
        p3B = gamma*p1B;
        pontofixoB = [p1B;p2B;p3B;p1B*ones(Ns-3,1)];

        JacobianoB = dhenonNws_4_pf(pontofixoB,Ns,c);

        autovaloresB = abs(eig(JacobianoB));    
        
        autovB(iteraws,iterag) = max(autovaloresB);
    end
    
    autovAA(iteraws,iterag) = autovA(iteraws,iterag);
    if autovA(iteraws,iterag)>1
        autovAp(iteraws,iterag) = NaN;
    else
        autovAp(iteraws,iterag) = -2;
    end
    autovBB(iteraws,iterag) = autovB(iteraws,iterag);
    if autovB(iteraws,iterag)>1
        autovBp(iteraws,iterag) = NaN;
    else
        autovBp(iteraws,iterag) = -2;
    end    
end
end

%%

figure
pcolor(w/pi,ganhos,autovAp')
title('Ponto fixo p1')
xlabel({'$G$'},'Interpreter','latex')
ylabel({'$w_{S}$'},'Interpreter','latex')
shading flat
grid on
set(gca,'layer','top','FontSize',24)

%%
for iteraws=1:Nws,
for iterag=1:Nganhos,  
    ganho = ganhos(iterag);
    ws = w(iteraws);
    c = poly([exp(1j*ws) exp(-1j*ws)]);
    c = ((poly([exp(1j*ws) exp(-1j*ws)]))/sum(c))*ganho;
    Ns=length(c);
    
        for condinicial=1:ninicial,
            
                    
        p1A = (-(1-b)+sqrt((1-b)^2+4*a*(ganho^2)))/(2*ganho^2);
        p2A = p1A;
        p3A = ganho*p1A;
        pontofixoA = [p1A;p2A;p3A;p1A*ones(Ns-3,1)];

            
            x=0.1*rand(Ns,1)+pontofixoA;
        
            if Ns==1
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
        
h(iteraws,iterag)=mean(h_vetor);
h_desvio(iteraws,iterag)=std(h_vetor);
end

iteraws
end

save('variables_lyapunov_teste_2_szpi')
%%

figure
pcolor(w/pi,ganhos,h')
ylabel({'$G$'},'Interpreter','latex')
xlabel({'$\frac{w_s}{\pi}$'},'Interpreter','latex')
shading flat
colorbar 
grid on
set(gca,'layer','top','FontSize',24)

%%

for iteraws=1:Nws,
    for iterag=1:Nganhos,
        if h(iteraws,iterag)<0 && h(iteraws,iterag)<h_desvio(iteraws,iterag) && h(iteraws,iterag)>-h_desvio(iteraws,iterag);
            h_color(iteraws,iterag)=-1;
        end
        if h(iteraws,iterag)<0 && h(iteraws,iterag)>h_desvio(iteraws,iterag) && h(iteraws,iterag)<-h_desvio(iteraws,iterag);
            h_color(iteraws,iterag)=0;
        end
        if isnan(h(iteraws,iterag))==1
            h_color(iteraws,iterag)=2;
        end             
    end
end


%%

for iteraws=1:Nws,
    for iterag=1:Nganhos,
        if h(iteraws,iterag)<0;
            h_color(iteraws,iterag)=-1;
        end
        if h(iteraws,iterag)>0 & h(iteraws,iterag)<0.002
            h_color(iteraws,iterag)=0;
        end
        if h(iteraws,iterag)>0.002
            h_color(iteraws,iterag)=1;
        end
        if isnan(h(iteraws,iterag))==1
            h_color(iteraws,iterag)=2;
        end             
    end
end

%%
figure
pcolor(w/pi,ganhos,h_color') 
ylabel({'$G$'},'Interpreter','latex')
xlabel({'$\frac{w_s}{\pi}$'},'Interpreter','latex')
shading flat
grid on
hold on
yline(0.4)
yline(0.6)
pcolor(w/pi,ganhos,autovAp')
shading flat
grid on
hold off
set(gca,'layer','top','FontSize',24)
legend

%%
figure
subplot(2,1,1)
pcolor(w/pi,ganhos,h')
ylabel({'$G$'},'Interpreter','latex')
xlabel({'$\frac{w_s}{\pi}$'},'Interpreter','latex')
shading flat
colorbar 
grid on
set(gca,'layer','top','FontSize',24)
subplot(2,1,2)
pcolor(w/pi,ganhos,h_color') 
ylabel({'$G$'},'Interpreter','latex')
xlabel({'$\frac{w_s}{\pi}$'},'Interpreter','latex')
shading flat
grid on
hold on
yline(0.4)
yline(0.6)
pcolor(w/pi,ganhos,autovAp')
shading flat
grid on
legend
hold off
set(gca,'layer','top','FontSize',24)


%%
for iteraws=1:Nws,
    for iterag=1:Nganhos,
        if h(iteraws,iterag)<0 & autovAp(iteraws,iterag)==-2;
            h_color_1(iteraws,iterag)=-2;
            h_color_2(iteraws,iterag)=NaN;
            h_color_3(iteraws,iterag)=NaN;
            h_color_4(iteraws,iterag)=NaN;
            h_color_5(iteraws,iterag)=NaN;
        end
        if h(iteraws,iterag)<0 && h(iteraws,iterag)<h_desvio(iteraws,iterag) && h(iteraws,iterag)>-h_desvio(iteraws,iterag);
            h_color_1(iteraws,iterag)=NaN;
            h_color_2(iteraws,iterag)=-1;
            h_color_3(iteraws,iterag)=NaN;
            h_color_4(iteraws,iterag)=NaN;
            h_color_5(iteraws,iterag)=NaN;
        end
        if h(iteraws,iterag)<0 && h(iteraws,iterag)>h_desvio(iteraws,iterag) && h(iteraws,iterag)<-h_desvio(iteraws,iterag);
            h_color_1(iteraws,iterag)=NaN;
            h_color_2(iteraws,iterag)=NaN;
            h_color_3(iteraws,iterag)=0;
            h_color_4(iteraws,iterag)=NaN;
            h_color_5(iteraws,iterag)=NaN;
        end        
        if h(iteraws,iterag)>0
            h_color_1(iteraws,iterag)=NaN;
            h_color_2(iteraws,iterag)=NaN;           
            h_color_3(iteraws,iterag)=NaN;
            h_color_4(iteraws,iterag)=1;
            h_color_5(iteraws,iterag)=NaN;
        end
        if isnan(h(iteraws,iterag))==1
            h_color_1(iteraws,iterag)=NaN;
            h_color_2(iteraws,iterag)=NaN;           
            h_color_3(iteraws,iterag)=NaN;
            h_color_4(iteraws,iterag)=NaN;
            h_color_5(iteraws,iterag)=2;
        end             
    end
end

%%
figure
subplot(2,1,1)
pcolor(w/pi,ganhos,h')
ylabel({'$G$'},'Interpreter','latex')
xlabel({'$\frac{w_s}{\pi}$'},'Interpreter','latex')
shading flat
colorbar 
grid on
set(gca,'layer','top','FontSize',24)
subplot(2,1,2)
pcolor(w/pi,ganhos,h_color_1') 
hold on
pcolor(w/pi,ganhos,h_color_2')
pcolor(w/pi,ganhos,h_color_3')
pcolor(w/pi,ganhos,h_color_4')
pcolor(w/pi,ganhos,h_color_5')
shading flat
ylabel({'$G$'},'Interpreter','latex')
xlabel({'$\frac{w_s}{\pi}$'},'Interpreter','latex')
shading flat
grid on
hold on
yline(0.4)
yline(0.6)
shading flat
grid on
legend
hold off
set(gca,'layer','top','FontSize',24)

%%
for iteraws=1:Nws,
    for iterag=1:Nganhos,
        if h(iteraws,iterag)<0 && autovAp(iteraws,iterag)==-2
            h_color_1(iteraws,iterag)=-2;
            h_color_2(iteraws,iterag)=NaN;
            h_color_3(iteraws,iterag)=NaN;
            h_color_4(iteraws,iterag)=NaN;
        end
        if h(iteraws,iterag)<0 && autovAp(iteraws,iterag)~=-2
            h_color_1(iteraws,iterag)=NaN;
            h_color_2(iteraws,iterag)=-1;
            h_color_3(iteraws,iterag)=NaN;
            h_color_4(iteraws,iterag)=NaN;
        end        
        if h(iteraws,iterag)>0 && autovAp(iteraws,iterag)~=-2
            h_color_1(iteraws,iterag)=NaN;
            h_color_2(iteraws,iterag)=NaN;           
            h_color_3(iteraws,iterag)=0;
            h_color_4(iteraws,iterag)=NaN;
        end
        if isnan(h(iteraws,iterag))==1
            h_color_1(iteraws,iterag)=NaN;
            h_color_2(iteraws,iterag)=NaN;           
            h_color_3(iteraws,iterag)=NaN;
            h_color_4(iteraws,iterag)=1;
        end             
    end
end


%%
figure
pcolor(w/pi,ganhos,h_color_1')
hold on
pcolor(w/pi,ganhos,h_color_2')
pcolor(w/pi,ganhos,h_color_3')
pcolor(w/pi,ganhos,h_color_4')
hold off
shading flat
ylabel({'$G$'},'Interpreter','latex')
xlabel({'$\frac{w_s}{\pi}$'},'Interpreter','latex')
shading flat
grid on
legend
set(gca,'layer','top','FontSize',24)

%%
figure
subplot(2,1,1)
pcolor(w/pi,ganhos,h')
ylabel({'$G$'},'Interpreter','latex')
xlabel({'$\frac{w_s}{\pi}$'},'Interpreter','latex')
shading flat
colorbar 
grid on
set(gca,'layer','top','FontSize',24)
subplot(2,1,2)
pcolor(w/pi,ganhos,h_color_1') 
hold on
pcolor(w/pi,ganhos,h_color_2')
pcolor(w/pi,ganhos,h_color_3')
pcolor(w/pi,ganhos,h_color_4')
shading flat
ylabel({'$G$'},'Interpreter','latex')
xlabel({'$\frac{w_s}{\pi}$'},'Interpreter','latex')
shading flat
grid on
hold on
yline(0.4)
yline(0.6)
shading flat
grid on
legend
hold off
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

w = eye(Ncoeficientes); %Initial orthonormal base

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
         
        %Standardization 
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

w = eye(Ncoeficientes); %Initial orthonormal base

    for i = 1:Nitera
        dx = dhenon_N_4(x,Ncoeficientes,coeficientes); 
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
