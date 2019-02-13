% HYPEINDICATOREXACT calculates the HypE fitness
%   HYPEINDICATOREXACT( POINTS, BOUNDS, K  ) calculates the HypE
%   fitness for objective vectors POINTS. 
%
%   POINTS: objective vectors as rows, all to be minimized
%   BOUNDS: reference point
%   K:      parameter of HypE
%
%   Example: f =  hypeIndicatorExact( [1 3; 3 1], [4 4], 1 )

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