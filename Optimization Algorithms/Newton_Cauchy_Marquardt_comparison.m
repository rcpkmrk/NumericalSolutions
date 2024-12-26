clc;
clear;

%% Function
f = @(X) (X(1) - 1)^2 + 5 * (X(2) - 1)^2 + 10; % Objective function
grad_f = @(X) [2 * (X(1) - 1); 10 * (X(2) - 1)]; % Gradient
hessian_f = @(X) [2, 0; 0, 10]; % Hessian (constant)

%% Parameters
X0 = [0; 0];         % Initial guess
tol = 1e-6;          % Convergence tolerance
max_iter = 100;      % Maximum iterations
lambda0 = 1;         % Initial damping parameter for Marquardt

%% Results

% a) Newton's Method
[X_min_Newton, f_min_Newton, iter_Newton] = newtons_method_multivar(f, grad_f, hessian_f, X0, tol);

% b) Steepest Descent (Cauchy) Method
[X_min_Steepest, f_min_Steepest, iter_Steepest] = steepest_descent(f, grad_f, X0, tol, max_iter);

% c) Marquardt Method
[X_min_Marquardt, f_min_Marquardt, iter_Marquardt] = marquardt_method(f, grad_f, hessian_f, X0, tol, max_iter, lambda0);

MethodNames = {'Newton'; 'Cauchy'; 'Marquardt'};
Iterations = [iter_Newton; iter_Steepest; iter_Marquardt];
Final_X1 = [X_min_Newton(1); X_min_Steepest(1); X_min_Marquardt(1)];
Final_X2 = [X_min_Newton(2); X_min_Steepest(2); X_min_Marquardt(2)];
Final_F = [f_min_Newton; f_min_Steepest; f_min_Marquardt];

ResultsTable = table(MethodNames, Iterations, Final_X1, Final_X2, Final_F, ...
    'VariableNames', {'Method', 'Iterations', 'X1', 'X2', 'F(X)'});
disp(ResultsTable);

%% Functions

% Newton's Method
function [X_min, f_min, iter] = newtons_method_multivar(f, grad_f, hessian_f, X0, tol)
    X = X0; % Initialize
    iter = 0;

    while true
        grad = grad_f(X); % Gradient
        hess = hessian_f(X); % Hessian

        % Newton's Step
        delta_X = -hess \ grad;
        X = X + delta_X;

        iter = iter + 1;
        if norm(delta_X) < tol % Convergence check
            break;
        end
    end

    X_min = X; % Result
    f_min = f(X_min); % Minimum value
end

% Steepest Descent (Cauchy) Method
function [X_min, f_min, iter] = steepest_descent(f, grad_f, X0, tol, max_iter)
    X = X0; % Initialize
    iter = 0;

    while iter < max_iter
        grad = grad_f(X); % Gradient
        if norm(grad) < tol % Convergence check
            break;
        end

        % Line Search
        alpha = fminbnd(@(a) f(X - a * grad), 0, 1); % Optimal step size
        X = X - alpha * grad; % Update

        iter = iter + 1;
    end

    X_min = X; % Result
    f_min = f(X_min); % Minimum value
end

% Marquardt Method
function [X_min, f_min, iter] = marquardt_method(f, grad_f, hessian_f, X0, tol, max_iter, lambda0)
    X = X0; % Initialize
    lambda = lambda0; % Initial damping
    iter = 0;

    while iter < max_iter
        grad = grad_f(X); % Gradient
        hess = hessian_f(X); % Hessian

        if norm(grad) < tol % Convergence check
            break;
        end

        % Modified Hessian
        H_mod = hess + lambda * eye(size(hess));
        delta_X = -H_mod \ grad; % Step

        % Update based on function improvement
        if f(X + delta_X) < f(X)
            X = X + delta_X; % Accept step
            lambda = lambda / 10; % Decrease damping
        else
            lambda = lambda * 10; % Increase damping
        end

        iter = iter + 1;
    end

    X_min = X; % Result
    f_min = f(X_min); % Minimum value
end
