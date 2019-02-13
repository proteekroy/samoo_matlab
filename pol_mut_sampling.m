
function [pop_mut] = pol_mut_sampling(pop_sel, N, p_mut, cum_dist, eta_m, Xmin, Xmax )

nreal = size(pop_sel,2);
pop_mut = zeros(N, size(pop_sel,2));

for j=1:N
    r = rand;
    p1 = find(r<=cum_dist);
    p1 = p1(1);
    x = pop_sel(p1,:);%take a solution according to roulette wheel selection
    
    for i = 1:nreal
        y = x(i);
        if rand <= p_mut
            yl = Xmin(i);
            yu = Xmax(i);
            delta1 = (y-yl) / (yu-yl);
            delta2 = (yu-y) / (yu-yl);
            rand_var = rand;
            mut_pow = 1.0/(eta_m+1.0);
            if rand_var <= 0.5
                xy = 1.0 - delta1;
                val = 2.0*rand_var + (1.0 - 2.0*rand_var) * xy^(eta_m+1.0);
                deltaq =  val^mut_pow - 1.0;
            else
                xy = 1.0 - delta2;
                val = 2.0*(1.0 - rand_var) + 2.0*(rand_var-0.5) * xy^(eta_m+1.0);
                deltaq = 1.0 - val^mut_pow;
            end
            y = y + deltaq*(yu - yl);
            if (y<yl), y = yl; end
            if (y>yu), y = yu; end
            pop_mut(j,i) = y;
        end
    end

end