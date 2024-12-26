clc; clear;

% the multi-variable function f(x, y)
f = @(x) (11 - x(1) - x(2))^2 + (1 + x(1) + 10*x(2) - x(1)*x(2))^2;

% Initial guesses
simplex = [0 0; 1 0; 0 1];  % Simplex vertices
X0 = [0; 0];  % Initial guess for Hooke-Jeeves method
tol = 1e-7;  % Convergence tolerance
max_iter = 10000;  % Maximum number of iterations
step_size = 1e-5;  % Step size for Hooke-Jeeves method
alpha = 850;  % Scaling factor for exploration (Hooke-Jeeves method)

% Nelder-Mead Method for multi-variable optimization
[x_min_Nelder, f_min_Nelder, iter_Nelder] = nelder_mead(f, simplex, tol, max_iter);

% Hooke-Jeeves Pattern for multi-variable optimization
[x_min_Hooke, f_min_Hooke, iter_Hooke] = hooke_jeeves(f, X0, tol, max_iter, step_size, alpha);

% Table
results = table( ...
    {'Nelder-Mead'; 'Hooke-Jeeves'}, ...  % Method names
    [iter_Nelder; iter_Hooke], ...  % Iteration counts
    {x_min_Nelder'; x_min_Hooke'}, ...  % Final x and y values
    [f_min_Nelder; f_min_Hooke], ...  % Final function values
    'VariableNames', {'Method', 'Iterations', 'x, y', 'f(x, y)'} ...
);

disp(results);

%% Nelder-Mead Method Function for multi-variable optimization
function [x_min, f_min, iter] = nelder_mead(f, simplex, tol, max_iter)
    iter = 0;  % Iteration counter
    
    % Loop until convergence or max iterations
    while iter < max_iter
        % Evaluate function values at simplex points
        f_values = arrayfun(@(i) f(simplex(i, :)'), 1:size(simplex, 1));
        
        % Sort points by their function values (ascending)
        [f_values, order] = sort(f_values);
        simplex = simplex(order, :);
        
        % Check for convergence: Standard deviation of function values
        if std(f_values) < tol
            break;
        end
        
        % Centroid of the best n points (excluding the worst)
        centroid = mean(simplex(1:end-1, :), 1);
        
        % Reflection
        reflected_point = centroid + (centroid - simplex(end, :));  % Reflection formula
        f_reflected = f(reflected_point');
        
        if f_reflected < f_values(1)
            % Expansion
            expanded_point = centroid + 2 * (reflected_point - centroid);  % Expansion formula
            f_expanded = f(expanded_point');
            
            if f_expanded < f_reflected
                simplex(end, :) = expanded_point;  % Use expanded point
            else
                simplex(end, :) = reflected_point;  % Use reflected point
            end
        elseif f_reflected < f_values(end-1)
            simplex(end, :) = reflected_point;  % Accept reflection
        else
            % Contraction
            contracted_point = centroid + 0.5 * (simplex(end, :) - centroid);  % Contraction formula
            f_contracted = f(contracted_point');
            
            if f_contracted < f_values(end)
                simplex(end, :) = contracted_point;  % Use contracted point
            else
                % Shrink simplex
                for i = 2:size(simplex, 1)
                    simplex(i, :) = simplex(1, :) + 0.5 * (simplex(i, :) - simplex(1, :));
                end
            end
        end
        
        iter = iter + 1;  % Increment iteration count
    end
    
    % Return the best point and its function value
    x_min = simplex(1, :)';  % Best point
    f_min = f(x_min);  % Best function value
end

%% Hooke-Jeeves Pattern Search Method Function for multi-variable optimization
function [x_min, f_min, iter] = hooke_jeeves(f, X0, tol, max_iter, step_size, alpha)
    x = X0;  % Initial point
    iter = 0;  % Iteration counter
    f_x = f(x);  % Initial function value
    prev_x = x;  % Previous point for convergence check

    while iter < max_iter
        % Exploration Step: Search along the axes (positive and negative)
        for i = 1:length(x)
            % Try both positive and negative steps
            x_new = x; 
            x_new(i) = x_new(i) + step_size;  % Move in the positive direction
            f_new = f(x_new);
            
            if f_new < f_x
                x = x_new;
                f_x = f_new;  % Update function value
            else
                x_new(i) = x_new(i) - 2 * step_size;  % Move in the negative direction
                f_new = f(x_new);
                if f_new < f_x
                    x = x_new;
                    f_x = f_new;  % Update function value
                end
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
        
        % Check for convergence: If the change in solution vector is small enough
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
