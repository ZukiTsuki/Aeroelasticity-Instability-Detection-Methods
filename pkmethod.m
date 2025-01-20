function [] = pkmethod()
%By Kazuki Arai
clear
clc
%verification
% a = -0.2;
% e = -0.1;
% mu = 20; %mass ratio
% r = sqrt(6/25); %r^2
% s = 0.4;
% x_theta = e-a;

n=1000;
mass=[526.708 446.708];
m=mass(2);
Ip=[7 4];
rho=[1.225 1.1346 1.091 0.9093];
l=sqrt(7.4*16.6);
c=16.6/l;
b=c/2;
a = -0.2;
e = [-0.3 -0.1];
x_theta = e(2)-a;
r=sqrt(Ip(2)/(m*b^2));
mu=m/(pi*rho(4)*l*b^2);
wth=(pi/2)*sqrt((10^5)/((Ip(2))*l));
wh=((1.8751)^2)*sqrt((2*10^5)/(m*l^3));
s=wh/wth;
z=zeros(n,4);
o=zeros(n,4);
k=0.01;
k_n=zeros(n,4);
gam=zeros(n,4);
V=linspace(0,5,n);

for i=2:n
    i
    Ck=(0.01365+0.2808i*k-(k^2)/2)/(0.01365+0.3455i*k-(k^2));
    f11=(s^2)/(V(i)^2)-((k^2)/mu)+ 2*1i*Ck/mu;
    f12=(k*(1i+a*k)+(2+1i*k*(1-2*a))*Ck)/mu;
    f21=(a*(k^2)-1i*k*(1+2*a)*Ck)/mu;
    f22=((8*mu*r^2)/(V(i)^2)+4*1i*(1+2*a)*(2*1i-k*(1-2*a))*Ck-k*(k-4*1i+8*a*(1+a*k)))/(8*mu);
    A=(r^2)-(x_theta^2);
    B=f22+f11*(r^2)-f21*x_theta-f12*x_theta;
    C=f11.*f22-f12.*f21;
    p=[A 0 B 0 C];
    y=roots(p);
    for j=1:4
        z(i,j)=real(y(j));
        o(i,j)=imag(y(j));
    end
    for j=1:4
        gam(i,j)=z(i,j)/abs(o(i,j));
        k_n(i,j)=o(i,j);       
        if abs((abs(k_n(i,j))-k))<10^-4
            continue;
        else
            t=0;
            while abs((abs(k_n(i,j))-k))>10^-4
                t=t+1;
                k=abs(k_n(i,j));
                Ck=(0.01365+0.2808i*k-(k^2)/2)/(0.01365+0.3455i*k-(k^2));
                f11=(s^2)/(V(i)^2)-((k^2)/mu)+ 2*1i*Ck/mu;
                f12=(k*(1i+a*k)+(2+1i*k*(1-2*a))*Ck)/mu;
                f21=(a*(k^2)-1i*k*(1+2*a)*Ck)/mu;
                f22=((8*mu*r^2)/(V(i)^2)+4*1i*(1+2*a)*(2*1i-k*(1-2*a))*Ck-k*(k-4*1i+8*a*(1+a*k)))/(8*mu);
                A=(r^2)-(x_theta^2);
                B=f22+f11*(r^2)-f21*x_theta-f12*x_theta;
                C=f11*f22-f12*f21;
                p=[A 0 B 0 C];
                if any(isinf(p)==1) | any(isnan(p)==1)
                    k_n(i,j)=rand;
                    break;
                else
                y=roots(p); 
                k_n(i,j)=imag(y(j));
                end
                if t>4
                    k_n(i,j)=rand;
                    break;
                end
            end
        end
    end
end


disp("Plotting, please wait...");
figure (1);
for ii = 2:n
    % for jj=1:4
        plot(V(ii),(o(ii,1)*V(ii)),'*','Color',"red");
        hold on;
        plot(V(ii),(o(ii,2)*V(ii)),'*','Color',"blue");
        hold on;        
        plot(V(ii),(o(ii,3)*V(ii)),'*','Color',"green");
        hold on;
        plot(V(ii),(o(ii,4)*V(ii)),'*','Color',"magenta");
        hold on;
        xlim auto;
        ylim auto;
    % end
end

grid on;
% legend('\omega/\omega_\theta','\omega/\omega_\theta','\omega/\omega_\theta','\omega/\omega_\theta');
xlabel('$\frac{U}{b\omega_\theta}$','Interpreter','latex');
ylabel('$\frac{\omega}{\omega_\theta}$','Interpreter','latex');
set(gca,'fontsize', 15);
hold off;

figure (2);
for ii = 2:n
    % for jj=1:4
        plot(V(ii),(gam(ii,1)),'o','Color',"red");
        hold on;
        plot(V(ii),(gam(ii,2)),'o','Color',"blue");
        hold on;
        plot(V(ii),(gam(ii,3)),'o','Color',"green");
        hold on;
        plot(V(ii),(gam(ii,4)),'o','Color',"magenta");
        hold on;
        xlim auto;
        ylim auto;
    % end
end

grid on;
% legend('\gamma_1','\gamma_2','\gamma_3','\gamma_4');
xlabel('$\frac{U}{b\omega_\theta}$','Interpreter','latex');
ylabel('\gamma');
set(gca,'fontsize', 15);
hold off;
%--------------------------------------------------------------------------
%Approach 2 to find k_F, (omega_F/omega_theta) and V_F
%--------------------------------------------------------------------------
%finding numerically k_F, (omega_F/omega_theta) and V_F
% temp1_1 = sign(o(:,1)); %"sign" function outputs -1 for a negative number;
% %+1 for a positive number and 0 for 0
% temp1_2 = diff(temp1_1); %look for "diff" function in MATLAB's help
% ind1 = find(abs(temp1_2)==2,1,'first'); %look for "find" function in MATLAB's help
% if ~isempty(ind1)
%     k_F_1 = k(ind1);
%     omega_ratio_1 = 1/sqrt(z(ind1,1));
%     V_F_1 = omega_ratio_1/k_F_1;
%     fprintf('First branch');
%     fprintf('The dimensionless flutter speed is %6.4f\n',V_F_1);
% end
% 
% temp2_1 = sign(o(:,2));
% temp2_2 = diff(temp2_1);
% ind2 = find(abs(temp2_2)==2,1,'first');
% if ~isempty(ind2)
%     k_F_2 = k(ind2);
%     omega_ratio_2 = 1/sqrt(z(ind2,2));
%     V_F_2 = omega_ratio_2/k_F_2;
%     fprintf('Second branch\n');
%     fprintf('The dimensionless flutter speed is %6.4f\n',V_F_2);
% end

end
