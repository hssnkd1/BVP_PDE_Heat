# BVP_PDE_Heat

clc; clear all;
close all;
rng('default');
profile on;
 
 
T =1;        %  Lenth of the T
n =63;
omgasize=n+1;
 
t = linspace(0, T, n);
x=linspace(0, 1, n);
% Exact Solution
g=sin((pi/2)*x); 

x0=0.5;
 xy=0;  
for j=1:n

    for i=1:n
      xy=0;          
    if(j >= i)
    for ki=0:10
       

        xy=xy + (  ( (exp(-((0.5+ki)^2)*(pi^2)*(t(j)-x(i))) )* ...
            ( cos(pi*(ki+0.5)*x0) ) /( ((0.5+ki))*(pi ))  ));
     end

        xy=1+ (2*xy);  
     end
    
  if(j < i)
      xy=0;  
  end
    
     if(xy < 0)
      xy=0;  
  end
    
    
    kxy(i)=xy;
   
    
    end
       A(j,:)=(1/n)*kxy;
end
An=null(A);
E=eye(n,n);
[U,S,V] = svd(A);
C=V*S*V';
Q=U*V';
f=A*g';
f_d1=f*0;
for i=1:n
f_d1(i) =f(i)+(0.008*rand);
end
delta=norm(f_d1-f);

u=zeros(n,1);
I=eye(n,n);
fin=f_d1;

% % Landweber Iteration methods
uCLI=u*0;
uMLI=u*0;
itMLI_1=u*0;
uMLI_1=u*0;
for i=1:500
   %CLI
   uCLI=uCLI+(0.01*(V*S*U')*(fin-(A*uCLI))); 
   itCLI(i)=norm((A*uCLI)-fin)^2;
   
   
   % MLI
   uMLI=uMLI+(0.01*(U*S*U')*(fin-(A*uMLI)));
   itMLI(i)=norm((A*uMLI)-fin)^2;
   
end

figure(1)
plot(g,'--','LineWidth',1,'DisplayName','Exact')
hold
plot(uCLI,'-','LineWidth',1,'DisplayName','u_\alpha,\alpha=0.01 and delta=0.04')
lgd = legend;
lgd.FontSize = 12;
xlabel('CLI', 'FontSize', 12)
 
figure(2)
plot(g,'--','LineWidth',1,'DisplayName','Exact')
hold
plot(uMLI,'-','LineWidth',1,'DisplayName','u_\alpha,\alpha=0.01 and \delta=0.04')
lgd = legend;
lgd.FontSize = 12;
xlabel('MLI', 'FontSize', 12)

figure(3)
plot(itCLI,'-')
hold
plot(itMLI,'-.') 
legend('CLI','MLI')
xlabel('Convergence rate', 'FontSize', 12)


