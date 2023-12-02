clc;
clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                                    DATA                                           %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[r, T_liq, X_liq] = LiqDat();
save('LiqDat.mat', 'r', 'T_liq', 'X_liq');
load('LiqDat.mat');

[r, T_vap, X_vap] = VapDat();
save('VapDat.mat', 'r', 'T_vap', 'X_vap');
load('VapDat.mat');

r_liq_matrix = repmat(r,1,6);
r_vap_matrix = repmat(r,1,15);
X_liq_one_row = reshape(X_liq', 1, []);
X_vap__one_row = reshape(X_vap', 1, []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                COEFFICIENTS FOR LIQUID/WATER REACTION MODEL                       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X_design_liq = zeros(126,6); %initial 126x6 matrix
j = 1;
for i = 1:length(r)*length(T_liq)
    X_design_liq(i,:) = [ones(size(r_liq_matrix(i))), exp(-2.25 * r_liq_matrix(i)), r_liq_matrix(i),...
        r_liq_matrix(i) .* T_liq(j), log(r_liq_matrix(i)), erf(r_liq_matrix(i))];
    if mod(i, 21) == 0
        j = j + 1;
    end 
end
coefficients_liq = (X_design_liq' * X_design_liq) \ (X_design_liq' * X_liq_one_row');
X_numerical_liq_one_row = coefficients_liq' * X_design_liq';
X_numerical_liq = reshape(X_numerical_liq_one_row, 21, 6).';

fprintf('  COEFFICIENTS FOR LIQUID/WATER REACTION MODEL\n');
fprintf('------------------------------------------------------------------------------------------\n');
coeff_names = {'a0:', 'a1:', 'a2:', 'a3:', 'a4:', 'a5:'};

for j = 1:6
    fprintf('%s %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f\n', coeff_names{j}, coefficients_liq(j,:));
end
fprintf('\n------------------------------------------------------------------------------------------\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                 COEFFICIENTS FOR VAPOR/WATER REACTION MODEL                       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X_design_vap = zeros(315,6); %initial 315x6 matrix
j = 1;
for i = 1:length(r)*length(T_vap)
    X_design_vap(i,:) = [ones(size(r_vap_matrix(i))), exp(-2.25 * r_vap_matrix(i)), r_vap_matrix(i), ...
        r_vap_matrix(i) .* T_vap(j), log(r_vap_matrix(i)), erf(r_vap_matrix(i))];
    if mod(i, 21) == 0
        j = j + 1;
    end 
end
coefficients_vap = (X_design_vap' * X_design_vap) \ (X_design_vap' * X_vap__one_row');
X_numerical_vap_one_row = coefficients_vap' * X_design_vap';
X_numerical_vap = reshape(X_numerical_vap_one_row, 21, 15).';

fprintf('  COEFFICIENTS FOR VAPOR/WATER REACTION MODEL\n');
fprintf('------------------------------------------------------------------------------------------\n');
coeff_names = {'a0:', 'a1:', 'a2:', 'a3:', 'a4:', 'a5:'};

for j = 1:6
    fprintf('%s %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f\n', coeff_names{j}, coefficients_vap(j,:));
end
fprintf('\n------------------------------------------------------------------------------------------\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                                  PLOTTING                                         %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
% Loop to create isoquant lines
for i = 1:length(T_liq)       
    plot(r, X_liq,'kx', 'LineWidth', 2); % experimental data
    hold on;
    plot(r, X_numerical_liq,"b", 'LineWidth', 2); % model data  
end
for i = 1:length(T_vap)       
    plot(r, X_vap,'kx', 'LineWidth', 2); % experimental data
    plot(r, X_numerical_vap,"r", 'LineWidth', 2); % model data  
end

xlabel('Reforming Steam/Water Rate (kmol H₂O/kmol Natural Gas)');
ylabel('% Conversion of CH₄ in Reformer');
ylim([0.5,1])
title('Isoquant Plot for Steam Methane Reforming (SMR)');
grid on;
hold off;
