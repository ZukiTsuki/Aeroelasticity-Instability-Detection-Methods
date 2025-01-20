function [] = pmethod()
%By Kazuki Arai

clear
clc
%verification
% a = -0.2;
% e = -0.1;
% mu = 20; %mass ratio
% r = sqrt(6/25); %r^2
% s = 0.4;
% xtheta = e-a;

mass=[526.708 446.708];
m=mass(1);
rho=[1.225 1.1346 1.091 0.9093];

l=sqrt(7.4*16.6);
c=16.6/l;
b=c/2;
a = -0.2;
e = [-0.3 -0.1];
Ip=[7 4];
r=sqrt(Ip(1)/(m*b^2));
mu=m/(pi*rho(2)*l*b^2);
wth=(pi/2)*sqrt((10^5)/(Ip(1)*l));
wh=((1.8751)^2)*sqrt((2*10^5)/(m*l^3));
s=wh/wth;
xtheta = e(1)-a;

n=1000;
z=zeros(n,2);
o=zeros(n,2);
V=linspace(0,5,n);


A=(r^2-xtheta^2);
B=(r^2./V.^2)-(1/mu)-(2*a/mu)+(s^2*r^2./V.^2)-(2*xtheta/mu);
C=(s^2*r^2./V.^4)-s^2./(mu.*V.^2)-(2*a*s^2)./(mu.*V.^2);

for i=2:n
    p=[A 0 B(i) 0 C(i)];
    y=roots(p);
    z(i,1)=real(y(1))*V(i);
    o(i,1)=imag(y(1))*V(i);
    z(i,2)=real(y(2))*V(i);
    o(i,2)=imag(y(2))*V(i);
    z(i,3)=real(y(3))*V(i);
    o(i,3)=imag(y(3))*V(i);
    z(i,4)=real(y(4))*V(i);
    o(i,4)=imag(y(4))*V(i);
end


figure (1);
scatter(V,o);
hold on;
grid on;
legend('\Omega/\omega_\theta','\Omega/\omega_\theta','\Omega/\omega_\theta','\Omega/\omega_\theta');
xlabel('$\frac{U}{b\omega_\theta}$','Interpreter','latex');
ylabel('$\frac{\Omega}{\omega_\theta}$','Interpreter','latex');
set(gca,'fontsize', 18);
hold off;
figure (2);
title('Graph of tan\lambdal and 3/7\lambdal');
scatter(V,z);
hold on;
grid on;
legend('\Gamma/\omega_\theta','\Gamma/\omega_\theta','\Gamma/\omega_\theta','\Gamma/\omega_\theta');
xlabel('$\frac{U}{b\omega_\theta}$','Interpreter','latex');
ylabel('$\frac{\Gamma}{\omega_\theta}$','Interpreter','latex');
set(gca,'fontsize', 18);
hold off;

% figure (3);
% scatter(V,o);
% hold on;
% scatter(V,z);
% grid on;
% legend('\Omega/\omega_\theta','\Omega/\omega_\theta','\Omega/\omega_\theta','\Omega/\omega_\theta','\Gamma/\omega_\theta','\Gamma/\omega_\theta','\Gamma/\omega_\theta','\Gamma/\omega_\theta');
% xlabel('$\frac{U}{b\omega_\theta}$','Interpreter','latex');
% ylabel('$\frac{\Gamma \:or\: \Omega}{\omega_\theta}$','Interpreter','latex');
% set(gca,'fontsize', 18);
% hold off;
% end
