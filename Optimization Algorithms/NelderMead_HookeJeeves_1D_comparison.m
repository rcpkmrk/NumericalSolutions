clc; clear;

% Define the objective function
f = @(x) (x - 1).^4 - x + 2 * exp(x);

% Initialize parameters for optimization methods
simplex = [0; 0.5; 1];  % Initial simplex points (for Nelder-Mead method)
X0 = 0;  % Initial guess for Hooke-Jeeves method
tol = 1e-7;  % Convergence tolerance
max_iter = 100000;  % Maximum number of iterations
step_size = 1e-6;  % Step size for Hooke-Jeeves method
alpha = 850;  % Scaling factor for exploration (Hooke-Jeeves method)
% important: playing with alpha might decrease iteration

%% Results

% a) Nelder-Mead Method
[x_min_Nelder, f_min_Nelder, iter_Nelder] = nelder_mead(f, simplex, tol, max_iter);

% b) Hooke-Jeeves Pattern Search Method
[x_min_Hooke, f_min_Hooke, iter_Hooke] = hooke_jeeves(f, X0, tol, max_iter, step_size, alpha);

results = table( ...
    {'Nelder-Mead'; 'Hooke-Jeeves'}, ...  % Method names
    [iter_Nelder; iter_Hooke], ...        % Iteration counts
    [x_min_Nelder; x_min_Hooke], ...      % Final x values
    [f_min_Nelder; f_min_Hooke], ...      % Final function values
    'VariableNames', {'Method', 'Iterations', 'x', 'f(x)'} ...
);

% Display the results table
disp(results);

%% Nelder-Mead Method Function
function [x_min, f_min, iter] = nelder_mead(f, simplex, tol, max_iter)
    iter = 0;  % Iteration counter
    
    % Loop until convergence or max iterations
    while iter < max_iter
        % Evaluate function values at simplex points
        f_values = arrayfun(f, simplex);
        
        % Sort points by their function values (ascending)
        [f_values, order] = sort(f_values);
        simplex = simplex(order);
        
        % Check for convergence: Simplex diameter
        if max(abs(simplex(1:end-1) - simplex(end))) < tol
            break;
        end
        
        % Calculate centroid of all but the worst point
        centroid = mean(simplex(1:end-1));
        
        % Reflection step
        reflected_point = centroid + (centroid - simplex(end));  % Reflection formula
        f_reflected = f(reflected_point);
        
        if f_reflected < f_values(1)  % Reflected point is better than the best point
            % Expansion step
            expanded_point = centroid + 2 * (reflected_point - centroid);  % Expansion formula
            f_expanded = f(expanded_point);
            
            if f_expanded < f_reflected
                simplex(end) = expanded_point;  % Use expanded point
            else
                simplex(end) = reflected_point;  % Use reflected point
            end
        elseif f_reflected >= f_values(end)  % Reflected point is worse than the worst
            % Contraction step
            contracted_point = centroid + 0.5 * (simplex(end) - centroid);  % Contraction formula
            f_contracted = f(contracted_point);
            
            if f_contracted < f_values(end)
                simplex(end) = contracted_point;  % Use contracted point
            else
                % Shrink the simplex towards the best point
                simplex(2:end) = simplex(1) + 0.5 * (simplex(2:end) - simplex(1));
            end
        else
            % Reflected point is better than the worst but not better than the best
            simplex(end) = reflected_point;  % Just replace the worst point
        end
        
        iter = iter + 1;  % Increment iteration count
    end
    
    % Return the best point and its function value
    x_min = simplex(1);  % The first point in the sorted simplex is the best
    f_min = f(x_min);
end

%% Hooke-Jeeves Pattern Search Method Function
function [x_min, f_min, iter] = hooke_jeeves(f, X0, tol, max_iter, step_size, alpha)
    x = X0;  % Initial point
    iter = 0;  % Iteration counter
    f_x = f(x);  % Initial function value
    prev_x = x;  % Previous point for convergence check

    while iter < max_iter
        % Exploration Step: Search along the axes (positive and negative)
        x_new = x + step_size;  % Move in the positive direction
        f_new = f(x_new);  % Evaluate the function
        
        if f_new < f_x  % If improvement, accept the new point
            x = x_new;
            f_x = f_new;  % Update function value
        else
            x_new = x - step_size;  % Move in the negative direction
            f_new = f(x_new);  % Evaluate the function
            if f_new < f_x  % If improvement, accept the new point
                x = x_new;
                f_x = f_new;  % Update function value
            end
        end
        
        % Pattern Search: Try to move in the direction of improvement
        direction = x - X0;  % Direction of improvement
        if norm(direction) > tol  % Check if there's any improvement
            x_new = x + alpha * direction;  % Move along the pattern direction
            f_new = f(x_new);  % Evaluate the function
            
            if f_new < f_x  % If improvement, accept the new point
                x = x_new;
                X0 = x;  % Update the base point for the pattern search
                f_x = f_new;  % Update function value
            end
        end
        
        % Check for convergence
        if norm(x - prev_x) < tol
            break;
        end
        
        prev_x = x;  % Update the previous point
        iter = iter + 1;  % Increment iteration counter
    end
    
    % Return the minimum point and its function value
    x_min = x;
    f_min = f_x;
end
