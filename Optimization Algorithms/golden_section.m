clc; clear;

f = @(x) (x - 1).^4 - x + 2 * exp(x);  % Function definition
a = 0;                                % Lower bound
b = 0.6;                              % Upper bound
tol = 1e-6;                           % Convergence tolerance
x = linspace(a, b, 100);
y = f(x);

[xmin, fmin, iter] = golden_section_min(f, a, b, tol);
fprintf('Minimum at x = %.6f, f(x) = %.6f, after %d iterations.\n', xmin, fmin, iter);

figure;
plot(x, y, 'b-', 'LineWidth', 1.5); 
hold on;
plot(xmin, fmin, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('x');
ylabel('f(x)');
title('Function Plot with Minimum Point');
grid on;
legend('f(x)', 'Minimum Point');

function [xmin, fmin, iter] = golden_section_min(f, a, b, tol)
    
    phi = (1 + sqrt(5)) / 2; % Golden ratio   
    x1 = a + (2 - phi) * (b - a);
    x2 = b - (2 - phi) * (b - a);
    f1 = f(x1);
    f2 = f(x2);    
    iter = 0;
    while (b - a) > tol
        iter = iter + 1;
        if f1 < f2
            b = x2; % Narrow interval from the right
            x2 = x1;
            f2 = f1;
            x1 = a + (2 - phi) * (b - a);
            f1 = f(x1);
        else
            a = x1; % Narrow interval from the left
            x1 = x2;
            f1 = f2;
            x2 = b - (2 - phi) * (b - a);
            f2 = f(x2);
        end
    end
    xmin = (a + b) / 2;
    fmin = f(xmin);
end
