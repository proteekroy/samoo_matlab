function [pop_xover] = sbx_sampling( pop_sel, N, p_cross, cum_dist,  eta_c, Xmin, Xmax )

%[N, nreal] = size(pop_sel); % Population size & Number of variables

nreal = size(pop_sel,2);
% N should be even to compare pairs
if mod(N,2) ~= 0
    N = N - 1;
end

pop_xover = zeros(N, size(pop_sel,2));



for ind=1:2:N
    
    r = rand;
    p1 = find(r<=cum_dist);
    p1 = p1(1);
    r = rand;
    p2 = find(r<=cum_dist);
    p2 = p2(1);
    
    %Parent Selection
    parent1 = pop_sel(p1,:);
    parent2 = pop_sel(p2,:);
    
    %Child Initialization
    child1 = zeros(1, nreal);
    child2 = zeros(1, nreal);
    if rand <= p_cross

        for i=1:nreal
            if rand <= 0.5
                if abs( parent1(i)-parent2(i) ) > 1e-14
                    if parent1(i) < parent2(i)
                        y1 = parent1(i);
                        y2 = parent2(i);
                    else
                        y1 = parent2(i);
                        y2 = parent1(i);
                    end
                    yl = Xmin(i);
                    yu = Xmax(i);
                    beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
                    alpha = 2.0 - beta^(-(eta_c+1.0));
                    rand_var = rand;
                    if rand_var <= (1.0/alpha)
                        betaq = (rand_var*alpha)^(1.0/(eta_c+1.0));
                    else
                        betaq = (1.0/(2.0 - rand_var*alpha))^(1.0/(eta_c+1.0));
                    end
                    c1 = 0.5*((y1+y2) - betaq*(y2-y1));
                    beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
                    alpha = 2.0 - beta^(-(eta_c+1.0));
                    if rand_var <= (1.0/alpha)
                        betaq = (rand_var*alpha)^(1.0/(eta_c+1.0));
                    else
                        betaq = (1.0/(2.0 - rand_var*alpha))^(1.0/(eta_c+1.0));
                    end
                    c2 = 0.5*((y1+y2)+betaq*(y2-y1));
                    if (c1 < yl), c1 = yl; end
                    if (c2 < yl), c2 = yl; end
                    if (c1 > yu), c1 = yu; end
                    if (c2 > yu), c2 = yu; end
                    if rand <= 0.5
                        child1(i) = c2;
                        child2(i) = c1;
                    else
                        child1(i) = c1;
                        child2(i) = c2;
                    end
                else
                    child1(i) = parent1(i);
                    child2(i) = parent2(i);
                end
            else
                child1(i) = parent1(i);
                child2(i) = parent2(i);
            end
        end
    else
        child1 = parent1;
        child2 = parent2;
    end
    pop_xover(ind,:) = child1;
    pop_xover(ind+1,:) = child2;
end
