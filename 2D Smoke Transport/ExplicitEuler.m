clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pi = 4*atan(1); wide = 2; high = 2;
dx = 0.02; dy = 0.02; dt = 0.005;
t_max=0.5;

M = round(t_max/dt + 1);     %# of nodes in time
Nx = round(wide/dx + 1);     %# of nodes in distance (x-axis)
Ny = round(high/dy + 1);     %# of nodes in distance (y-axis)

t = zeros(M,1);
x = zeros(Nx,1); 
y = zeros(Ny,1); 

for i=1:M                %time node values
    t(i) = 0 + (i-1)*dt;
end   

for i=1:Nx  
    for j=1:Ny  
    x(i) = 0 + (i-1)*dx;
    y(j) = 0 + (j-1)*dy;
    u(i,j) = cos(pi*x(i)).*sin(pi*y(j));
    v(i,j) = -sin(pi*x(i)).*cos(pi*y(j));
    end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 0.5*ones(Nx,Ny);   %Initial conditions
c(:,1) = 0;        %Boundary conditions
c(:,Ny) = 0;        %Boundary conditions
c_initial = c;       %Initial concentration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                        %//FTCS Equation loop//%
for j=2:Ny-2
   for i=2:Nx-2                
      c(i+1,j+1) = c(i,j) - ((u(i,j)*dt)/(2*dx))*(c(i+1,j)-...
          c(i-1,j)) - ((v(i,j)*dt)/(2*dy))*(c(i,j+1)-c(i,j-1));                
   end
end
T=c;      
surf(x,y,T)
colorbar;
xlabel('x')
ylabel('y')
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 