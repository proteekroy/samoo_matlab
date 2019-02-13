function F = find_pareto(pop)

N = size( pop, 1 );
F = [];
ncon = 0;

for i = 1:N

    individual(i).n = 0;
    individual(i).S = [];
    for j = 1:N
        if j ~= i
            f = check_dominance(pop(i,:), pop(j,:), ncon);
            if f == 1
                individual(i).S = [individual(i).S j];
            elseif f == -1
                individual(i).n = individual(i).n+1;
                break
            end
        end
    end

    if individual(i).n == 0
        F = [F; i];
    end
    
end
end

function f = check_dominance(a, b, ncon)
% f = 1 if a dominates b
% f = -1 if b dominates a
% f = 0 if both a and b are non-dominated
flag1 = 0;
flag2 = 0;
M = size(a,2);
if ncon > 0 % Constrained Dominance Check
    
    if a(end) > b(end) && all( [a(end), b(end)] ) % If CV(a) > CV(b): Select b
        f = -1;
    elseif (a(end) <= b(end)) && all( [a(end), b(end)] )  % If CV(a) <= CV(b): Select a
        f = 1;
    else % If both solutions are feasible: Apply usual domination check
        for i = 1:M
            if a(i) < b(i)
                flag1 = 1;
            elseif a(i) > b(i)
                flag2 = 1;
            end
        end
        if flag1 == 1 && flag2 == 0
            f = 1;
        else
            if flag1 == 0 && flag2 == 1
                f = -1;
            else
                f = 0;
            end
        end
    end
    
else % Unconstrained Dominance Check
    
    for i = 1:M
        if a(i) < b(i)
            flag1 = 1;
        elseif a(i) > b(i)
            flag2 = 1;
        end
    end
    if flag1 == 1 && flag2 == 0
        f = 1;
    else
        if flag1 == 0 && flag2 == 1
            f = -1;
        else
            f = 0;
        end
    end
    
end
end

