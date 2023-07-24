clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 1;                    %Length
alpha = 0.1;              %diffusion coefficient
dx = 0.05;                %Δx
dt = 0.01;                %Δt
d = (alpha*dt)/(dx^2);    %diffusion number
%for t_max = 0:0.1:0.5     %time interval loop(0.0,0.1,0.2,0.3,0.4,0.5)
t_max=0.5;
 M = round(t_max/dt + 1); %# of nodes in time
 N = round(L/dx + 1);     %# of nodes in distance
 t = zeros(M,1);
 x = zeros(N,1);          
 for i=1:N                %distance node values
     x(i) = 0 + (i-1)*dx;
 end    
 for i=1:M                %time node values
     t(i) = 0 + (i-1)*dt;
 end  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 T = 100*ones(M,N);   %Initial conditions
 T(:,1) = 300;        %Boundary conditions
 T(:,N) = 300;        %Boundary conditions
 T(1,2:N-1) = 100;    %Initialconditions
 T_initial = T;       %Initial temperatures
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
                        %//FTCS Equation loop//%
      for n=1:M
        for i=2:N-1
            T(n+1,i) = d*T(n,i+1) + (1-2*d)*T(n,i) + d*T(n,i-1);
        end
      end
      FTCS = T;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
                   %//DuFort-Frankel Equation loop//%
      for y=2:N-1
          T(2,y) = d*T(1,y+1) + (1-2*d)*T(1,y) + d*T(1,y-1);
      end    
      for n=2:M-1
        for i=2:N-1
            T(n+1,i) = ((1-2*d)/(1+2*d))*T(n-1,i) + ((2*d)/(1+2*d))*(T(n,i+1) + T(n,i-1));
        end
      end
      DF = T;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                   %//Laasonen(BTCS) Equation loop//%
A = zeros(N-2,N-2);
B = zeros(N-2,1);
A(1,1) = 2*(1+d);
A(1,2) = -d;
A(N-2,N-3) = -d;
A(N-2,N-2) = 2*(1+d);
for p=2:N-3
    A(p,p-1) = -d;
    A(p,p) = (1+2*d);
    A(p,p+1) = -d;
end    
for n=1:M-1
    B(1) = (T(n,2)+T(n+1,1));
    B(N-2) = (T(n,N-1)+T(n+1,N));
    for i=2:N-3
        B(i) = T(n,i+1);
    end   
    W = A\B;  %solution matrix
    T(n+1,2:N-1) = W;
    BTCS = T;
end         
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 
                   %//Crank-Nicolson Equation loop//%
A = zeros(N-2,N-2);
B = zeros(N-2,1);
A(1,1) = 2*(1+d);
A(1,2) = -d;
A(N-2,N-3) = -d;
A(N-2,N-2) = 2*(1+d);
for p=2:N-3
    A(p,p-1) = -d;
    A(p,p) = (1+2*d);
    A(p,p+1) = -d;
end    
for n=1:M-1
    B(1) = d*T(n,1)+2*(1-d)*T(n,2)+d*T(n,3)+d*T(n+1,1);
    B(N-2) = d*T(n,2)+2*(1-d)*T(n,3)+d*T(n,4)+d*T(n+1,N);
    for i=2:N-3
        B(i) = T(n,i+1);
    end
    W = A\B;  %solution matrix
    T(n+1,2:N-1) = W;
    CN = T;
end 
%}      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %//FTCS/Analytic Equation loop//%
                     T(1,2:N-1) = 90; 
      for n=1:M
        for i=2:N-1
            T(n+1,i) = d*T(n,i+1) + (1-2*d)*T(n,i) + d*T(n,i-1);
        end
      end
      AN = T;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                          %//Error//%
E_FTCS = (FTCS - AN)/AN(M,N);
E_DF = (DF - AN)/AN(M,N);
E_BTCS = (BTCS - AN)/AN(M,N);
E_CN = (CN - AN)/AN(M,N);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %Plots for first 2 questions%
plot(E_FTCS(M,:),x,"-o",'LineWidth',2); 
hold on;
plot(E_DF(M,:),x,"-^",'LineWidth',2); 
plot(E_BTCS(M,:),x,"-square",'LineWidth',2); 
plot(E_CN(M,:),x,"-diamond",'LineWidth',2); 
xlabel('Relative Error')
ylabel('x(m)')
title(sprintf('x vs. Error with dt=%.2f,time=0.5',dt))
set(gca,'FontSize',20)
legend('FTCS','DuFort-Frankel','Laasonen','Crank-Nicolson','Analytical')
title(legend,'Methods')
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
                      %Plots for 3rd question%
r = 1;
while (r<(N/2+1))
plot(t,T(:,r),'LineWidth',2);
hold on;
r = r + 1;
end
xlabel('time(h)')
ylabel('Temperature(°C)')
title(sprintf('Temperature vs. time with dt=%.2f(Crank-Nicolson)',dt))
set(gca,'FontSize',20)
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%end %end of row8 for loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%{                      
                      %//Analitic Equation loop//%            
        
      time=0;
      for n=2:M-1  
          f=0;
        for i=2:N-1      
                syms m  ;
                v = exp((-((m*pi)/L)^2)*alpha*time)*((1-(-1).^m)/(m*pi))*sind((m*pi*f)/L);
                sum_final = symsum(v,m,1,10);
           
            T(n,i) = 300 + 2*(100 - 300)*sum_final;
            %{       
            T(n,i) = 1/(sqrt(4*pi*d*time))*exp(-(f^2)/(4*d*time));
            %}
            f = f + dx;
        end
        
         time = time + dt;
      end
plot(x,T(M,:),'LineWidth',2); 
xlabel('x(m)')
ylabel('Temperature(°C)')
title(sprintf('Temperature vs. x with dt=%.2f(FTCS)',dt))
set(gca,'FontSize',20)
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%