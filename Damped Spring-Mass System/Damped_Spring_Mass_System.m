clc;clear;
h = 0.1;t0 = 0;x0 = 1;v0 = 0;tn = 15; %initial values
n = (tn-t0)/h;
t(1) = t0; x(1) = x0; v(1) = v0;      %initial values
c=5;                                  %damping coefficient
f = @(t,x,v) v; g = @ (t,x,v) -x-v/c; %equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%                   RK2                       %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
for i=1:n
    t(i+1) = t0 + i*h;
    k1 = h*f(t(i),x(i),v(i));
    m1 = h*g(t(i),x(i),v(i));
    k2 = h*f(t(i+1),x(i)+k1,v(i)+k1);    
    m2 = h*g(t(i+1),x(i)+k1,v(i)+m1);
    x(i+1) = x(i) + (1/2)*(k1+k2);
    v(i+1) = v(i) + (1/2)*(m1+m2);    
end    
plot(t,x,'LineWidth',2)
hold on
xlabel('t(s)')
ylabel('x(m)')
title('RK2 method Solution with different damping coefficients')
set(gca,'FontSize',20)
legend('c=5','c=40','c=200')
title(legend,'damping coeff.')
set(gcf, 'PaperPosition', [-4 4 30 20])
set(groot, 'DefaultFigurePosition', [500 600 1000 700])
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%                  Analytical                 %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%                c=5                  %%%%%%%%%%%%%%%%%%%%

t_a = 0:0.1:15;
x_a = exp(-t_a/8).*(cos((3*sqrt(7).*t_a)/8)+(1/(3*sqrt(7)))*sin((3*sqrt(7).*t_a)/8));
plot(t_a,x_a,'LineWidth',2)
hold on
%%%%%%%%%%%%%%%%%%                c=40                 %%%%%%%%%%%%%%%%%%%%
t_a = 0:0.1:15;
x_a = exp(-1*t_a)+exp(-1*t_a).*t_a;
plot(t_a,x_a,'LineWidth',2)
%%%%%%%%%%%%%%%%%%                c=200                %%%%%%%%%%%%%%%%%%%%
t_a = 0:0.1:15;
x_a = (1+((-2*sqrt(6)+5)/(4*sqrt(6))))*exp((-5+2*sqrt(6))*t_a)...
    -((-2*sqrt(6)+5)/(4*sqrt(6)))*exp((-5-2*sqrt(6))*t_a);
plot(t_a,x_a,'LineWidth',2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%          Analytical Solution Plot            %%%%%%%%%%%%%%%%
xlabel('t(s)')
ylabel('x(m)')
title('Analytical Solution with different damping coefficients')
set(gca,'FontSize',20)
legend('c=5(underdamped)','c=40(critically damped)','c=200(overdamped)')
title(legend,'damping coefficient')
set(gcf, 'PaperPosition', [-4 4 30 20])
set(groot, 'DefaultFigurePosition', [500 600 1000 700])
%}