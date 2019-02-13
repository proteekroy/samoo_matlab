function pop = combine_and_reduce_model(opt, pop1, pop2)

[N, k] = size(pop1); k = k-2; % N = 2 * N_initial

pop2 = [pop2, zeros(N, 2)];
pop2 = evaluate_fitness_model(opt, pop2(:, 1:k));

pop3 = [pop1; pop2];
pop3 = pop3( randperm(2*N), :); % Shuffle pop3

tour = randperm(2*N);

% START THE COMPETITION
pop = zeros(N, k+2);
ID = zeros(N,1);
count = 0;
for i=1:2:2*N
    fitness = pop3( tour(i:i+1), k+1 );
    cons_violation = pop3( tour(i:i+1), k+2 );

    if (cons_violation(1) <= 0) && (cons_violation(2) <= 0)
        % If both solutions are feasible: Select the one with the minimum fitness value
        min_member = find(fitness == min(fitness));
    elseif (cons_violation(1) >= 0) && (cons_violation(2) <= 0)
        % If 2nd solution is feasible and 1st solution is infeasible: Select 2nd solution
        min_member = 2;
    elseif (cons_violation(1) <= 0) && (cons_violation(2) >= 0)
        % If 1st solution is feasible and 2nd solution is infeasible: Select 1st solution
        min_member = 1;
    elseif (cons_violation(1) >= 0) && (cons_violation(2) >= 0)
        % If both solutions are feasible: Select the one with the minimum constraint violation
        min_member = find(cons_violation == min(cons_violation));
    end
    count = count + 1;
    ID(count) = tour(i-1+min_member(1));
    pop(count, :) = pop3( ID(count), :);
end


