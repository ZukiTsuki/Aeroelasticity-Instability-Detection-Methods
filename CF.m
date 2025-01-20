function [] = CF()
clear
clc
%%verification
a = -0.2;
e = -0.1;
mu = 20; %mass ratio
r = sqrt(6/25); %r^2
s = 0.4;
n=101;

% mass=[526.708 446.708];
% m=mass(1);
% Ip=[7 4];
% rho=[1.225 1.1346 1.091 0.9093];
% 
% l=sqrt(7.4*16.6);
% c=16.6/l;
% b=c/2;
% a = 0.6;
% e = 0.525;
% r=sqrt(Ip(1)/m*b^2);
% mu=m/(pi*rho(1)*b^2);
% wth=(pi/2)*sqrt(10^5/(7*b));
% wh=((1.8751)^3)*sqrt((2*10^5)/m*l^3);
% s=wth/wh;

z=zeros(n,2);
o=zeros(n,2);
Ck=zeros(1,n);
xtheta = e-a;
k=linspace(0,1,n);

Ck=(0.01365+0.2808i.*k-(k.^2)/2)./(0.01365+0.3455i.*k-(k.^2));

lh=1- 2*1i*Ck./k;
lt=-a-(1i./k)-(2*Ck./k.^2)-(2*1i*(0.5-a)*Ck./k);
mh=-a+(2*1i*(0.5+a)*Ck./k);
mt=(1/8)+a^2-1i*(0.5-a)./k+2*(0.5+a)*Ck./k.^2+2*1i*(0.25-a^2)*Ck./k;

A=(mu^2)*(r^2)*(s^2);
B=-(mu^2)*(r^2)-mu*(r^2)*(s^2)-mu*mt*(s^2)-mu*(r^2).*lh;
C=(mu^2)*(r^2)+mu*mt+mu*(r^2)*lh+mt.*lh-(mu^2)*(xtheta^2)-mu*xtheta*mh-mu*xtheta*lt-lt.*mh;


for i=2:n
    p=[A B(i) C(i)];
    y=roots(p);
    z(i,1)=real(y(1));
    o(i,1)=imag(y(1));
    z(i,2)=real(y(2));
    o(i,2)=imag(y(2));
end

% figure (1);
% scatter(V,o);
% hold on;
% grid on;
% legend('\Omega/\omega_\theta','\Omega/\omega_\theta','\Omega/\omega_\theta','\Omega/\omega_\theta');
% xlabel('$\frac{U}{b\omega_\theta}$','Interpreter','latex');
% ylabel('$\frac{\Omega}{\omega_\theta}$','Interpreter','latex');
% set(gca,'fontsize', 18);
% hold off;
% figure (2);
% title('Graph of tan\lambdal and 3/7\lambdal');
% scatter(V,z);
% hold on;
% grid on;
% legend('\Gamma/\omega_\theta','\Gamma/\omega_\theta','\Gamma/\omega_\theta','\Gamma/\omega_\theta');
% xlabel('$\frac{U}{b\omega_\theta}$','Interpreter','latex');
% ylabel('$\frac{\Gamma}{\omega_\theta}$','Interpreter','latex');
% set(gca,'fontsize', 18);
% hold off;

figure (3);
scatter(k,o);
hold on;
scatter(k,z);
grid on;
legend('\Omega/\omega_\theta','\Omega/\omega_\theta','\Omega/\omega_\theta','\Omega/\omega_\theta');
xlabel('$\frac{b\omega}{U}$','Interpreter','latex');
ylabel('$\frac{\omega_\theta}{\omega}$','Interpreter','latex');
set(gca,'fontsize', 18);
hold off;

