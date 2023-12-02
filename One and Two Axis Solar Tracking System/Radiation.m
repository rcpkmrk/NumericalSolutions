clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%finding day number algorithm
dateStr = input('Enter the date (dd/mm/yyyy): ', 's'); %date
dateVec = datevec(dateStr, 'dd/mm/yyyy'); %datevec: MATLAB function that converts date to vector
day = dateVec(3) ; month = dateVec(2); year = dateVec(1);
daysInMonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
if is29februaryYear(year)
        daysInMonth(2) = 29;
end
dateNum = sum(daysInMonth(1:month-1)) + day;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startSolarTime = input('Enter the start solar time (HH:MM): ', 's'); %start time
endSolarTime = input('Enter the end solar time (HH:MM): ', 's'); %end time
latitude = input('Enter the latitude (in degrees): '); 
reflectance = input('Enter the reflectance value of the site: ');

startHour = str2double(startSolarTime(1:2)); %parsing start hour
endHour = str2double(endSolarTime(1:2)); %parsing end hour
timeHours = startHour:endHour;

A = 1160 + 75 * sind((360 / 365) * (dateNum - 275)); %coefficient A
k = 0.174 + 0.035 * sind((360 / 365) * (dateNum - 100)); %coefficient k
delta = 23.45 * sind((360 / 365) * (dateNum + 284)); %coefficient δ

% Initial arrays
totalInsolation_2axis = zeros(endHour-startHour+1, 1);
beamInsolation_2axis = zeros(endHour-startHour+1, 1);
diffuseInsolation_2axis = zeros(endHour-startHour+1, 1);
reflectedInsolation_2axis = zeros(endHour-startHour+1, 1);

totalInsolation_1axis = zeros(endHour-startHour+1, 1);
beamInsolation_1axis = zeros(endHour-startHour+1, 1);
diffuseInsolation_1axis = zeros(endHour-startHour+1, 1);
reflectedInsolation_1axis = zeros(endHour-startHour+1, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation for Two-Axis Tracking Solar Insolation
disp('********************************************');
disp('Two-Axis Tracking Solar Insolation:');
disp('********************************************');
disp('|--------|--------|--------|--------|');
disp('|--Ic----|---Ibc--|---Idc--|---Irc--|');
disp('|--------|--------|--------|--------|');
for i = 1:(endHour-startHour+1)
    H = (12 - timeHours(i)) * 15;
    B = cosd(latitude) * cosd(delta) * cosd(H) + sind(latitude) * sind(delta);
    Bn = asind(B);
    m = 1 / sind(Bn);
    Ibc = A * exp(-k * m);
    C = 0.095 + 0.04 * sind((360 / 365) * (dateNum - 100));
    Idc = C * Ibc * ((1 + cosd(90 - Bn)) / 2);
    Irc = reflectance * (Ibc * cosd(delta) + Idc);
    Ic = Ibc + Idc + Irc;
    
    % Results for Two-Axis
    totalInsolation_2axis(i) = Ic;
    beamInsolation_2axis(i) = Ibc;
    diffuseInsolation_2axis(i) = Idc;
    reflectedInsolation_2axis(i) = Irc;

    fprintf('|%8.3f|%8.3f|%8.3f|%8.3f|\n', Ic, Ibc, Idc, Irc);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Calculation for One-Axis Tracking Solar Insolation
disp('********************************************');
disp('One-Axis Tracking Solar Insolation:');
disp('********************************************');
disp('|--------|--------|--------|--------|');
disp('|--Ic----|---Ibc--|---Idc--|---Irc--|');
disp('|--------|--------|--------|--------|');
for i = 1:(endHour-startHour+1)
    H = (12 - timeHours(i)) * 15;
    B = cosd(latitude) * cosd(delta) * cosd(H) + sind(latitude) * sind(delta);
    Bn = asind(B);
    m = 1 / sind(Bn);
    Ibc = beamInsolation_2axis(i)*cosd(delta);
    C = 0.095 + 0.04 * sind((360 / 365) * (dateNum - 100));
    Idc = C * beamInsolation_2axis(i) * ((1 + cosd(latitude)) / 2);
    Irc = reflectance * (Ibc * cosd(delta) + Idc);
    Ic = Ibc + Idc + Irc;
    
    % Results for One-Axis
    totalInsolation_1axis(i) = Ic;
    beamInsolation_1axis(i) = Ibc;
    diffuseInsolation_1axis(i) = Idc;
    reflectedInsolation_1axis(i) = Irc;

    fprintf('|%8.3f|%8.3f|%8.3f|%8.3f|\n', Ic, Ibc, Idc, Irc);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(timeHours,totalInsolation_1axis,'LineWidth',2);
hold on;
plot(timeHours, totalInsolation_2axis,'LineWidth',2);
hold off;
xlabel('Time(HH)')
ylabel('Radiation(W/m^2)')
title('Total Insolation')
grid on;
legend('One-Axis','Two-Axis')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(timeHours,beamInsolation_1axis,'LineWidth',2);
hold on;
plot(timeHours, beamInsolation_2axis,'LineWidth',2);
hold off;
xlabel('Time(HH)')
ylabel('Radiation(W/m^2)')
title('Beam Insolation')
grid on;
legend('One-Axis','Two-Axis')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(timeHours,diffuseInsolation_1axis,'LineWidth',2);
hold on;
plot(timeHours, diffuseInsolation_2axis,'LineWidth',2);
hold off;
xlabel('Time(HH)')
ylabel('Radiation(W/m^2)')
title('Diffuse Insolation')
grid on;
legend('One-Axis','Two-Axis')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(timeHours,reflectedInsolation_1axis,'LineWidth',2);
hold on;
plot(timeHours, reflectedInsolation_2axis,'LineWidth',2);
hold off;
xlabel('Time(HH)')
ylabel('Radiation(W/m^2)')
title('Reflected Insolation')
grid on;
legend('One-Axis','Two-Axis')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function is29february = is29februaryYear(year)
    is29february = (mod(year, 4) == 0 && mod(year, 100) ~= 0) || (mod(year, 400) == 0);
end
