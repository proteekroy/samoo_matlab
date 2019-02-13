function [lm] = lebesgue_measure(F, ub)
% [lm] = lebesgue_measure(F, ub)
%
% Computes the hypervolume (or Lebesgue measure) of the M x l
% matrix F of l vectors of M objective function values.
%
% Implementation of the Lebesgue Measure Algorithm as described by:
%
% 'M. Fleischer. The measure of Pareto Optima Applications to 
%  Multi-objective Metaheuristics. EMO 2003, LNCSS 2632
%  519-533, 2003.'
%
% IMPORTANT:
%   Considers Minimization of the objective function values!
%   This function assumes that all solutions of F are non-dominated!
%
% Input:
% - F              - An M x l matrix where each of the l columns 
%                    represents a vector of M objective function values
%                    (note that this function assumes that all 
%                    solutions of F are non-dominated).
% - ub             - Optional: Upper bound reference point (default:
%                    the boundary point containing the maximum of F
%                    for each objective).
%
% Output:
% - lm             - The hypervolume (or lebesgue measure) of F.
%
% Author: Johannes W. Kruisselbrink
% Last modified: March 17, 2011

    F = F(paretoFront(F),:);
    
    F = F';

	%if (nargin < 2), ub = (max(F,[],2)); end
    %ub = max(F,[],2) ;%+ offset_fraction * (max(F,[],2) - min(F,[],2));
    
	[M, l] = size(F);
	if (M == 2)
		lm = lebesgue_2D(F, ub);
	else
		%lm = lebesgue_ND(F, ub);
        lm = hypeIndicatorExact( F', ub', 1);
        lm = sum(lm);
        %lm = approximate(F', ub, 1000);
	end

end

function f = hypeIndicatorExact( points, bounds, k )
    Ps = size(points,1);
    if( k < 0 )
        k = Ps;
    end
    actDim = size(points,2);
    if( length(bounds) == 1 )
        bounds = repmat(bounds, actDim, 1);
    end
    pvec = 1 : 1 : size(points,1);
    
	alpha = [];
    for i = 1 : k
        j = 1:i-1;
        alpha(i) = prod( (k-j) ./ (Ps - j ) )./i;
    end  
    
    f = hypesub( size(points,1), points, actDim, bounds, pvec, alpha, k );
end
    
function [h] = hypesub( l, A, actDim, bounds, pvec, alpha, k )
        h = zeros(1,l);
        [S,i] = sortrows(A,actDim);
        pvec = pvec(i);
        for i = 1 : size(S,1)
            if( i < size(S,1) )
                extrusion = S( i+1, actDim ) - S( i, actDim );
            else
                extrusion = bounds(actDim) - S( i, actDim );
            end
            if( actDim == 1 )
                if( i > k )
                    break;
                end                
                if( alpha >= 0 )
                    h( pvec(1:i) ) = h( pvec(1:i) ) + extrusion * alpha(i);
                end                
            elseif( extrusion > 0 )
                h = h + extrusion*hypesub( l, S(1:i,:), actDim - 1, ...
                        bounds, pvec(1:i), alpha, k );            
            end        
        end
end

function v = approximate(P, r, N)

    P=P*diag(1./r);
    [n,d]=size(P);
    C=rand(N,d);
    
    fDominated=false(N,1);
    lB=min(P);
    fcheck=all(bsxfun(@gt, C, lB),2);

    for k=1:n
        if any(fcheck)
            f=all(bsxfun(@gt, C(fcheck,:), P(k,:)),2);
            fDominated(fcheck)=f;
            fcheck(fcheck)=~f;
        end
    end

    v=sum(fDominated)/N;

end


function [lm] = lebesgue_2D(F, ub)
% Efficient method for 2D objective function values
	L = sortrows(F',1)';
	l = length(L(1,:));
	lm = 0;
	for i = 1:l
		lm = lm + ((L(1,i) - ub(1)) * (L(2,i) - ub(2)));
		ub(2) = L(2,i);
	end
end


function [lm] = lebesgue_ND(F, ub)
% Method for ND objective function values as described in:
%
% 'M. Fleischer. The measure of Pareto Optima Applications to 
%  Multi-objective Metaheuristics. EMO 2003, LNCSS 2632
%  519-533, 2003.'
%
% Author: Johannes W. Kruisselbrink
% Last modified: March 17, 2011

	lm = 0;

	[M, l] = size(F);

	% Remove the duplicates from F and compute Lebesque measure of L
	L = unique(F','rows')';

	while l >= 1
		if (l > 1)
			b = zeros(M,1);
			spawn_vector = repmat(L(:,1), 1, M);
			for i = 1:M
				% Bound b(i) is either the least upper bound of the i-th value of the
				% other points or it is the value of the absolute upper bound ub(i)
				difL = (L(i,2:end) - L(i,1));
				lub = find(difL > 0);
				if (length(lub) > 0)
					b(i) = min(L(i,lub+1));
				else
					b(i) = ub(i);
				end

				b(i) = min((difL > 0) .* L(i,2:end) + (difL <= 0) * ub(i));
				% Update i-th spawn vector
				spawn_vector(i,i) = b(i);
			end

			% Compute lop-Off volume and update lebesgue measure
			lov = prod(b - L(:,1));
			lm = lm + lov;

			% Remove L(:,1) from L
			L = L(:,2:end);

			% Add the spawn_vector to L, but first filter dominated
			% solutions and the solutions that touch the upper bounds
			% from the spawn_vector
			L = nd_filter(L, spawn_vector, ub);

		else
			lov = prod(ub - L(:,1));
			lm = lm + lov;
			L = [];
		end
		% Update l
		[M, l] = size(L);
	end

end


function L = nd_filter(L, spawn_vector, ub)
% Implementation of the filter routine as described in:
%
% 'M. Fleischer. The measure of Pareto Optima Applications 
%  to Multi-objective Metaheuristics. EMO 2003, LNCSS 2632
%  519-533, 2003.'
%
% Author: Johannes W. Kruisselbrink
% Last modified: March 17, 2011

	[M, l_L] = size(L);
	[M, l_sp] = size(spawn_vector);
	do_assign = zeros(1, l_sp);

	for i = 1 : l_sp

		% Find if the spawnvector hits the upper bound
		at_ub = false;
		for j = 1:M
			if (spawn_vector(j,i) == ub(j))
				at_ub = true;
				break;
			end
		end

		% For this if statement, the following would be more elegant
		% (replacing the loop above), but less efficient:
		% if (all(spawn_vector(:,i) ~= ub))
		if (at_ub == false)
			do_assign(i) = 1;
			for j = 1 : l_L
				if (weakly_dominates(L(:,j), spawn_vector(:,i)))
					do_assign(i) = 0;
					break;
				end
			end
		end
	end
	L = [spawn_vector(:,find(do_assign == 1)), L];

end
