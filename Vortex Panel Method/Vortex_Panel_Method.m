clc; clear;

%% Parameters
rho = 1.225; % Air density at sea level (kg/m^3)
mu = 1.789e-5; % Dynamic viscosity at sea level (Pa·s)
Re = 500000; % Reynolds number
chord = 1; % Airfoil chord length (m)
Umag = Re * mu / (rho * chord); % Freestream velocity magnitude

% Airfoil coordinates
coordinates = load('NASA 64-012.txt');
coordinates = flip(coordinates);

NumPanels = length(coordinates) - 1; % Number of panels

%% Preprocessing

for i=1:NumPanels
    point1(i,1) = coordinates(i,1);
    point1(i,2) = coordinates(i,2);
    point2(i,1) = coordinates(i+1,1);
    point2(i,2) = coordinates(i+1,2);
end
for i=1:NumPanels
    dx = point2(i,1)-point1(i,1);
    dy = point2(i,2)-point1(i,2);
    dl(i) = sqrt( (point2(i,1)-point1(i,1))^2 + (point2(i,2)-point1(i,2))^2); % Panel lengths
    th(i) = atan2(dy,dx); % Panel angles
    tnx(i) = cos(th(i)); % Tangential unit vectors (x-component)
    tny(i) = sin(th(i)); % Tangential unit vectors (y-component)
    vnx(i) = -tny(i); % Normal unit vectors (x-component)
    vny(i) = tnx(i); % Normal unit vectors (y-component)
end
for i=1:NumPanels
    midpoints(i,1) = 0.5*(point1(i,1)+point2(i,1));
    midpoints(i,2) = 0.5*(point1(i,2)+point2(i,2));
end

%% Influence Coefficients

for i=1:NumPanels
    
    for j=1:NumPanels
        
        xt = midpoints(i,1) - point1(j,1);
        yt = midpoints(i,2) - point1(j,2);
        x = xt*tnx(j) + yt*tny(j);
        y = -xt*tny(j) + yt*tnx(j);
        x1 = 0.0; y1 = 0.0;
        x2t = point2(j,1) - point1(j,1);
        y2t = point2(j,2) - point1(j,2);
        x2 = x2t*tnx(j) + y2t*tny(j);
        r1 = sqrt(x^2+y^2);
        r2 = sqrt((x-x2)*(x-x2)+y^2);
        th1 = atan2(y,x);
        th2 = atan2(y,x-x2);
        if(i==j) % self-induced velocity
            ax1 = 0.5*(x/x2-1.0);
            ay1 = 1.0/(2*pi);
            ax2 =-0.5*x/x2;
            ay2 =-1.0/(2*pi);
        else
            dth = th2-th1;
            rrt = r2/r1;
            rrtl = log(rrt);
            fcc = 1/(2*pi*x2);
            ax1 = fcc*( y*rrtl + (x-x2)*dth );
            ay1 = fcc*((x-x2)*rrtl - y*dth + x2);
            ax2 = -fcc*(y*rrtl + x*dth );
            ay2 = -fcc*(x*rrtl - y*dth + x2);
        end
        
        ux1 = ax1*tnx(j) - ay1*tny(j);
        uy1 = ax1*tny(j) + ay1*tnx(j);
        ux2 = ax2*tnx(j) - ay2*tny(j);
        uy2 = ax2*tny(j) + ay2*tnx(j);
        
        if(j==1)
            a(i,1)= ux1*vnx(i) + uy1*vny(i);
            holda = ux2*vnx(i) + uy2*vny(i);
        elseif(j==NumPanels)
            a(i,NumPanels) = ux1*vnx(i) + uy1*vny(i) + holda;
            a(i,NumPanels+1) = ux2*vnx(i) + uy2*vny(i);
        else
            a(i,j)= ux1*vnx(i) + uy1*vny(i) + holda;
            holda = ux2*vnx(i) + uy2*vny(i);
        end
        
        if(j==1)
            b(i,1)= ux1*tnx(i) + uy1*tny(i);
            holdb = ux2*tnx(i) + uy2*tny(i);
        elseif(j==NumPanels)
            b(i,NumPanels) = ux1*tnx(i) + uy1*tny(i) + holdb;
            b(i,NumPanels+1) = ux2*tnx(i) + uy2*tny(i);
        else
            b(i,j)= ux1*tnx(i) + uy1*tny(i) + holdb;
            holdb = ux2*tnx(i) + uy2*tny(i);
        end
    end
end
a(NumPanels+1,1) = 1.0;%the Kutta condition
a(NumPanels+1,NumPanels+1) = 1.0;%the Kutta condition

%% Solving for a given AoA range
alpha = -5:1:15; % angle of attack in degrees
for k = 1:length(alpha)
    alpha_rad = alpha(k)*pi/180.0;
    coalf = cos(alpha_rad); sialf = sin(alpha_rad);
    Ux = Umag*coalf; Uy = Umag*sialf;
    for i=1:NumPanels
        rhs(i) = -Ux*vnx(i)-Uy*vny(i);
    end
    rhs(NumPanels+1)=0.0;
    gamma = rhs/a';
    circ = 0.0; % circulation
    circgam = 0.0; % circulation in terms of gamma
    for i=1:NumPanels
        tnvel = Ux*tnx(i)+Uy*tny(i); % tangential velocity
        for j=1:NumPanels+1
            tnvel = tnvel + b(i,j)*gamma(j);
        end
        circ = circ - tnvel*dl(i); % circulation
        circgam = circgam +0.5*(gamma(i)+gamma(i+1))*dl(i); % circulation
        cp(i) = 1.0-tnvel*tnvel/(Umag*Umag); % pressure coefficient
    end
    cp(NumPanels+1)=cp(1);
    cy = 0.0;
    cx = 0.0;
    for i=1:NumPanels
        cy = cy-cp(i)*dl(i)*vny(i);
        cx = cx-cp(i)*dl(i)*vnx(i);    
    end
    cy = cy/chord;
    cx = cx/chord;
    cd(k) =  cx*coalf+cy*sialf; %drag coefficient
    cl(k) = -cx*sialf+cy*coalf; %lift coefficient
end    

% Plot Lift Coefficient (Cl)
figure;
plot(alpha, cl, 'b-', 'LineWidth', 2);
grid on;
title('Lift Coefficient vs. Angle of Attack', 'FontSize', 14); 
xlabel('Angle of Attack, \alpha (degrees)', 'FontSize', 12);
ylabel('Lift Coefficient, C_l', 'FontSize', 12);
legend('C_l', 'FontSize', 10, 'Location', 'best');
set(gca, 'FontSize', 12);

% Plot Drag Coefficient (Cd)
figure;
plot(alpha, cd, 'r-', 'LineWidth', 2); 
grid on;
title('Drag Coefficient vs. Angle of Attack', 'FontSize', 14);
xlabel('Angle of Attack, \alpha (degrees)', 'FontSize', 12); 
ylabel('Drag Coefficient, C_d', 'FontSize', 12);
legend('C_d', 'FontSize', 10, 'Location', 'best');
set(gca, 'FontSize', 12);

