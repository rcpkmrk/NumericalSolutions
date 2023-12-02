clc; clear;

populationSize = 20;
population = initialpop();
initial_population = population;
generations = 2000;
tolerance = 0.000001;
m3contour;
% Main loop for generations
for gen = 1:generations
    phi_prev = GPF(population(:,1),population(:,2));
    average_prevgen_phi = mean(phi_prev);
    pairs = randperm(populationSize);    
    offspring = zeros(populationSize / 2, 2);
    offspringIndex = 1;
    for i = 1:2:populationSize
        parent1 = population(pairs(i), :);
        parent2 = population(pairs(i+1), :);     
        r1 = (rand() - 0.5) * 0.5; 
        r2 = (rand() - 0.5) * 0.5; 
        r = [r1, r2];
        m = r*(1/(10*gen));
        offspring(offspringIndex, :) = (parent1+parent2)/2 + m;        
        offspringIndex = offspringIndex + 1;
    end
    population = [population; offspring];
    phi = GPF(population(:,1),population(:,2));
    population_phi = [population, phi];
    population_phi_sorted = sortinghat(population_phi,3);

    next_gen = population_phi_sorted(11:end, :);
    average_nextgen_phi = mean(next_gen(:, 3));
    population = next_gen(:, [1, 2]);
    if(abs(average_nextgen_phi-average_prevgen_phi) < tolerance)
        break;
    end    
    scatter(next_gen(:, 1), next_gen(:, 2), 'bo');
    hold on;
end
averageGPF = [mean(next_gen(:, 1)),mean(next_gen(:, 2)),mean(next_gen(:, 3))];
fprintf('Average of Points = (%.4f, %.4f), Average phi = %.4f\n', averageGPF(1, 1), averageGPF(1, 2), averageGPF(1, 3));
fprintf('Generation = %.4f\n', gen);