% wfgextend.M
%
% Scalable test problem toolkit.
%
% F = wfgextend(Z, M, k, l, testNo,S_bounds)
%
% Inputs: Z      - candidate solutions
%         M      - no. of objectives
%         k      - no. of position-related parameters
%         l      - no. of distance-related parameters
%         testNo - test problem no.
%         S_bounds - the range of PFs  e.g., S_bounds = 2*(1:M)
%
% Output: F      - objective vectors
%

%
% wfgextend is an extended version of WFG test problems. 
% in the wfgextend a set of new problems with various PF geometries are
% constructed. WFG41-WFG48. These problems are created by applying
% different shape functions provided in the WFG suite to the original WFG4
% problem.

% wfg41 is WFG4
% wfg42 is WFG4 with a convex shape
% wfg43 is WFG4 with a sharp concave shape
% wfg44 is WFG4 with a sharp convex shape
% wfg45 is WFG4 with a mixed shape
% wfg46 is WFG4 with a linear shape
% wfg47 is WFG7 with a disconnected convex shape
% wfg48 is WFG8 with a disconnected concave shape
% The above WFG4X problems can be used to test the performance of
% algorithms on various PF geometries. 

% Please cite the following papers when you are using this code
% Ref: Huband S, Hingston P, Barone L, While L, 2006, A review
%      of multiobjective test problems and a scalable test problem
%      toolkit. IEEE Transactions on Evolutionary Computation,
%      10(5), pp477-506.
% Ref: Wang, R., Purshouse, R. C., Fleming, P. J., Preference-inspired 
%      co-evolutionary algorithms using weight vectors, European Journal
%      of Operational Research, 243(2), 423-441, 2015.
% Ref: Wang, R., Zhang, Q. F., Zhang, T., Decomposition based algorithms 
%      using Pareto adaptive scalarizing methods IEEE Transactions on 
%      Evolutionary Computation, IEEE, 2015, To appear


%
% Robin Purshouse, Rui Wang, 27 May 2012
% Contact: ruiwangnudt@gmail.com

function F = wfgextend(Z, M, k, l, testNo,S_bounds)

% Check for correct number of inputs.
if nargin < 6
    % Could also check data input range z_i in [0, 2i].
    S = NaN * ones(1, M);
    for i = 1:M
        S(i) = 2*i;
    end
    S_bounds = S;
end
if size(S_bounds,2)~=M
    error('Inconsistent: S_bounds with M.');
end
% Get total number of decision variables and no. of candidates.
[noSols, n] = size(Z);

% Data input checks.
if n ~= (k + l)
    error('Inconsistent number of variabes.');
end
if rem(k,M-1) ~= 0
    error('k must be divisible by M-1.');
end
if (testNo == 2 || testNo == 3) && (rem(l,2) ~= 0)
    error('l must be a multiple of 2 for WFG2 and WFG3.');
end
% Initialise function-wide constants.
NO_TESTS = 48;
S = S_bounds;

D = 1;
A = ones(NO_TESTS, M-1);
A(3,2:M-1) = 0;

% Transform all variable ranges to [0 1].
Y = Z;
for i = 1:n
    Y(:,i) = Z(:,i) ./ (2*i);
end

if testNo == 1

    % Apply first transformation.
    Ybar = Y;
    lLoop = k + 1;
    shiftA = 0.35;
    Ybar(:, lLoop:n) = s_linear(Ybar(:, lLoop:n), shiftA);

    % Apply second transformation.
    Ybarbar = Ybar;
    biasA = 0.8;
    biasB = 0.75;
    biasC = 0.85;
    Ybarbar(:, lLoop:n) = b_flat(Ybarbar(:, lLoop:n), biasA, biasB, ...
        biasC);

    % Apply third transformation.
    Ybarbarbar = Ybarbar;
    biasA = 0.02;
    Ybarbarbar = b_poly(Ybarbarbar, biasA);

    % Apply fourth transformation.
    T = NaN * ones(noSols, M);
    uLoop = M - 1;
    for i = 1:uLoop
        lBnd = 1+(i-1)*k/(M-1);
        uBnd = i*k/(M-1);
        weights = 2*(lBnd:uBnd);
        T(:,i) = r_sum(Ybarbarbar(:,lBnd:uBnd), weights);
    end
    T(:,M) = r_sum(Ybarbarbar(:,lLoop:n), 2*(lLoop:n));

    % Apply degeneracy constants.
    X = T;
    for i = 1:M-1
        X(:,i) = max(T(:,M), A(testNo, i)) .* (T(:,i) - 0.5) + 0.5;
    end

    % Generate objective values.
    fM = h_mixed(X(:,1), 1, 5);
    F = h_convex(X(:,1:uLoop));
    F(:,M) = fM;
    F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* F;

elseif testNo == 2 || testNo == 3

    % Apply first transformation.
    Ybar = Y;
    lLoop = k + 1;
    shiftA = 0.35;
    Ybar(:, lLoop:n) = s_linear(Ybar(:, lLoop:n), shiftA);

    % Apply second transformation.
    Ybarbar = Ybar;
    uLoop = k + l/2;
    for i = lLoop:uLoop
        lBnd = k+2*(i-k)-1;
        uBnd = k+2*(i-k);
        Ybarbar(:,i) = r_nonsep(Ybar(:,lBnd:uBnd), 2);
    end

    % Apply third transformation.
    T = NaN * ones(noSols, M);
    uLoop = M - 1;
    weights = ones(1, k/(M-1));
    for i = 1:uLoop
        lBnd = 1+(i-1)*k/(M-1);
        uBnd = i*k/(M-1);
        T(:,i) = r_sum(Ybarbar(:,lBnd:uBnd), weights);
    end
    T(:,M) = r_sum(Ybarbar(:,lLoop:k+l/2), ones(1, (k+l/2)-lLoop+1));

    % Apply degeneracy constants.
    X = T;
    for i = 1:M-1
        X(:,i) = max(T(:,M), A(testNo, i)) .* (T(:,i) - 0.5) + 0.5;
    end

    % Generate objective values.
    if testNo == 2
        fM = h_disc(X(:,1), 1, 1, 5);
        F = h_convex(X(:,1:uLoop));
        F(:,M) = fM;
        F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* F;
    else
        F = h_linear(X(:,1:uLoop));
        F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* F;
    end

elseif testNo == 4 || testNo==47 || testNo==48 ||  testNo == 41 || testNo == 42 || testNo == 43 || testNo == 44 || testNo == 45 || testNo == 46 || testNo == 5

    % Apply first transformation.
    if testNo == 4 || testNo==47 || testNo==48 || testNo == 41 || testNo == 42 || testNo == 43 || testNo == 44 || testNo == 45 || testNo == 46
        shiftA = 30;
        shiftB = 10;
        shiftC = 0.35;
        Ybar = s_multi(Y, shiftA, shiftB, shiftC);
    elseif testNo == 5
        shiftA = 0.35;
        shiftB = 0.001;
        shiftC = 0.05;
        Ybar = s_decep(Y, shiftA, shiftB, shiftC);
    end
    % Apply second transformation.
    T = NaN * ones(noSols, M);
    lLoop = k + 1;
    uLoop = M - 1;
    weights = ones(1, k/(M-1));
    for i = 1:uLoop
        lBnd = 1+(i-1)*k/(M-1);
        uBnd = i*k/(M-1);
        T(:,i) = r_sum(Ybar(:,lBnd:uBnd), weights);
    end
    T(:,M) = r_sum(Ybar(:,lLoop:n), ones(1, n-lLoop+1));

    % Apply degeneracy constants.
    X = T;
    for i = 1:M-1
        X(:,i) = max(T(:,M), A(testNo, i)) .* (T(:,i) - 0.5) + 0.5;
    end
    % Generate objective values.
    if  testNo == 4 || testNo == 5
        F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* h_concave(X(:,1:uLoop));
    elseif testNo==41
        F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* h_concave(X(:,1:uLoop));
    elseif testNo==42
        F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* h_convex(X(:,1:uLoop));
    elseif testNo==43
        F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* h_concave_sharp(X(:,1:uLoop));
    elseif testNo==44
        F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* h_convex_sharp(X(:,1:uLoop));
    elseif testNo == 45
        fM = h_mixed(X(:,1), 1, 2);
        F = h_concave(X(:,1:uLoop));
        F(:,M) = fM;
        F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* F;
    elseif testNo==46
        F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* h_linear(X(:,1:uLoop));
    elseif testNo == 47
        fM = h_disc(X(:,1), 1/2, 1/2, 2);
        F = h_concave(X(:,1:uLoop));
        F(:,M) = fM;
        F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* F;
    elseif testNo==48
        %F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .*  h_disc(X(:,1:uLoop));
        %fM = h_disc(X(:,1), 1, 1, 2);
        fM = h_disc(X(:,1), 1/2, 1/2, 2);
        F = h_convex(X(:,1:uLoop));
        F(:,M) = fM;
        F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* F;
    
    end
elseif testNo == 6
    % Apply first transformation.
    Ybar = Y;
    lLoop = k + 1;
    shiftA = 0.35;
    Ybar(:, lLoop:n) = s_linear(Ybar(:, lLoop:n), shiftA);

    % Apply second transformation.
    T = NaN * ones(noSols, M);
    uLoop = M - 1;
    for i = 1:uLoop
        lBnd = 1+(i-1)*k/(M-1);
        uBnd = i*k/(M-1);
        T(:,i) = r_nonsep(Ybar(:,lBnd:uBnd), k/(M-1));
    end
    T(:,M) = r_nonsep(Ybar(:,k+1:k+l), l);

    % Apply degeneracy constants.
    X = T;
    for i = 1:M-1
        X(:,i) = max(T(:,M), A(testNo, i)) .* (T(:,i) - 0.5) + 0.5;
    end

    % Generate objective values.
    
    F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* h_concave(X(:,1:uLoop));
    

elseif testNo == 7

    % Apply first transformation.
    Ybar = Y;
    biasA = 0.98 / 49.98;
    biasB = 0.02;
    biasC = 50;
    for i = 1:k
        Ybar(:,i) = b_param(Y(:,i), r_sum(Y(:,i+1:n), ones(1, n-i)), ...
            biasA, biasB, biasC);
    end

    % Apply second transformation.
    Ybarbar = Ybar;
    lLoop = k + 1;
    shiftA = 0.35;
    Ybarbar(:, lLoop:n) = s_linear(Ybar(:, lLoop:n), shiftA);

    % Apply third transformation.
    T = NaN * ones(noSols, M);
    lLoop = k + 1;
    uLoop = M - 1;
    weights = ones(1, k/(M-1));
    for i = 1:uLoop
        lBnd = 1+(i-1)*k/(M-1);
        uBnd = i*k/(M-1);
        T(:,i) = r_sum(Ybarbar(:,lBnd:uBnd), weights);
    end
    T(:,M) = r_sum(Ybarbar(:,lLoop:n), ones(1, n-lLoop+1));

    % Apply degeneracy constants.
    X = T;
    for i = 1:M-1
        X(:,i) = max(T(:,M), A(testNo, i)) .* (T(:,i) - 0.5) + 0.5;
    end

    % Generate objective values.
    
    F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* h_concave(X(:,1:uLoop));

elseif testNo == 8

    % Apply first transformation.
    Ybar = Y;
    lLoop = k + 1;
    biasA = 0.98 / 49.98;
    biasB = 0.02;
    biasC = 50;
    for i = lLoop:n
        Ybar(:,i) = b_param(Y(:,i), r_sum(Y(:,1:i-1), ones(1, i-1)), ...
            biasA, biasB, biasC);
    end

    % Apply second transformation.
    Ybarbar = Ybar;
    shiftA = 0.35;
    Ybarbar(:, lLoop:n) = s_linear(Ybar(:, lLoop:n), shiftA);

    % Apply third transformation.
    T = NaN * ones(noSols, M);
    lLoop = k + 1;
    uLoop = M - 1;
    weights = ones(1, k/(M-1));
    for i = 1:uLoop
        lBnd = 1+(i-1)*k/(M-1);
        uBnd = i*k/(M-1);
        T(:,i) = r_sum(Ybarbar(:,lBnd:uBnd), weights);
    end
    T(:,M) = r_sum(Ybarbar(:,lLoop:n), ones(1, n-lLoop+1));

    % Apply degeneracy constants.
    X = T;
    for i = 1:M-1
        X(:,i) = max(T(:,M), A(testNo, i)) .* (T(:,i) - 0.5) + 0.5;
    end
    % Generate objective values.
    
    F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* h_concave(X(:,1:uLoop));
    
elseif testNo == 9

    % Apply first transformation.
    Ybar = Y;
    uLoop = n - 1;
    biasA = 0.98 / 49.98;
    biasB = 0.02;
    biasC = 50;
    for i = 1:uLoop
        Ybar(:,i) = b_param(Y(:,i), r_sum(Y(:,i+1:n), ones(1, n-i)), ...
            biasA, biasB, biasC);
    end

    % Apply second transformation.
    Ybarbar = Ybar;
    biasA = 0.35;
    biasB = 0.001;
    biasC = 0.05;
    Ybarbar(:,1:k) = s_decep(Ybar(:,1:k), biasA, biasB, biasC);
    biasA = 30;
    biasB = 95;
    biasC = 0.35;
    Ybarbar(:,k+1:n) = s_multi(Ybar(:,k+1:n), biasA, biasB, biasC);

    % Apply third transformation.
    T = NaN * ones(noSols, M);
    uLoop = M - 1;
    for i = 1:uLoop
        lBnd = 1+(i-1)*k/(M-1);
        uBnd = i*k/(M-1);
        T(:,i) = r_nonsep(Ybarbar(:,lBnd:uBnd), k/(M-1));
    end
    T(:,M) = r_nonsep(Ybarbar(:,k+1:k+l), l);

    % Apply degeneracy constants.
    X = T;
    for i = 1:M-1
        X(:,i) = max(T(:,M), A(testNo, i)) .* (T(:,i) - 0.5) + 0.5;
    end

    % Generate objective values.
    
    F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* h_concave(X(:,1:uLoop));
    
 
else
    error('Requested test function not implemented.');
end
end

% Transformation functions.

% Reduction: weighted sum.
function ybar = r_sum(y, weights)

[noSols noY] = size(y);
wgtMatrix=rep(weights,[noSols 1]);
ybar = y .* wgtMatrix;
ybar = sum(ybar, 2) ./ sum(wgtMatrix, 2);

end

% Reduction: non-separable.
function y_bar = r_nonsep(y, A)

[noSols noY] = size(y);

y_bar = 0;
for j = 1:noY
    innerSum = 0;
    for k = 0:(A-2)
        innerSum = innerSum + abs(y(:,j) - y(:,1+mod(j+k,noY)));
    end
    y_bar = y_bar + y(:,j) + innerSum;
end
y_bar = y_bar / ( (noY/A) * ceil(A/2) * (1+2*A-2*ceil(A/2)) );

end

% Bias: polynomial.
function y_bar = b_poly(y, alpha)

y_bar = y.^alpha;

end

% Bias: flat region
function y_bar = b_flat(y, A, B, C)

[noSols noY] = size(y);
min1 = min(0, floor(y - B));
min2 = min(0, floor(C - y));
y_bar = A + min1*A.*(B-y)/B - min2*(1-A).*(y-C)/(1-C);

% Machine precision problems can cause y_bar to go slightly negative so
% force >=0 condition.
y_bar=max(0,y_bar);

end

% Bias: parameter dependent.
function ybar = b_param(y, uy, A, B, C)

[noSols noY] = size(y);
v = A - (1 - 2*uy) .* abs(floor(0.5 - uy) + A);
v = rep(v, [1 noY]);
ybar = y.^(B + (C-B)*v);

end

% Shift: linear.
function ybar = s_linear(y, A)

ybar = abs(y - A) ./ abs(floor(A - y) + A);

end

% Shift: deceptive.
function ybar = s_decep(y, A, B, C)

y1 = floor(y - A + B) * (1 - C + (A - B)/B) / (A - B);
y2 = floor(A + B - y) * (1 - C + (1 - A - B)/B) / (1 - A - B);
ybar = 1 + (abs(y - A) - B) .* (y1 + y2 + 1/B);

end

% Shift: multi-modal.
function ybar = s_multi(y, A, B, C)

y1 = abs(y-C) ./ (2*(floor(C-y)+C));
ybar = (1 + cos((4*A+2)*pi*(0.5 - y1)) + 4*B*y1.^2) / (B+2);

end

% Shape functions.
% Linear.
function f = h_linear(x)

[noSols mMinusOne] = size(x);

M = mMinusOne + 1;
f = NaN * ones(noSols, M);

f(:,1) = prod(x, 2);
for i = 2:mMinusOne
    f(:,i) = prod(x(:,1:M-i), 2) .* (1 - x(:,M-i+1));
end
f(:,M) = 1 - x(:,1);
end

% Convex.
function f = h_convex(x)
[noSols mMinusOne] = size(x);

M = mMinusOne + 1;
f = NaN * ones(noSols, M);

f(:,1) = prod(1-cos(x*pi/2), 2);
for i = 2:mMinusOne
    f(:,i) = prod(1-cos(x(:,1:M-i)*pi/2), 2) .* (1-sin(x(:,M-i+1)*pi/2));
end
f(:,M) = 1-sin(x(:,1)*pi/2);
end

% Concave.
function f = h_concave(x)

[noSols mMinusOne] = size(x);

M = mMinusOne + 1;
f = NaN * ones(noSols, M);

f(:,1) = prod(sin(x*pi/2), 2);
for i = 2:mMinusOne
    f(:,i) = prod(sin(x(:,1:M-i)*pi/2), 2) .* cos(x(:,M-i+1)*pi/2);
end
f(:,M) = cos(x(:,1)*pi/2);

end

function f = h_concave_sharp(x)
[noSols mMinusOne] = size(x);
M = mMinusOne + 1;
f = NaN * ones(noSols, M);
wRate = 1/4;
f(:,1) = prod(sin(x*pi/2), 2);
for i = 2:mMinusOne
    f(:,i) = prod( sin( x(:,1:M-i)*pi/2), 2 ) .* cos(x(:,M-i+1)*pi/2).^wRate;
end
f(:,M) = cos(x(:,1)*pi/2).^wRate;
end

% Added by Rui: to generate a sharp shape
% Concave 21/12/2011.
function f = h_convex_sharp(x)
[noSols mMinusOne] = size(x);
M = mMinusOne + 1;
f = NaN * ones(noSols, M);
wRate = 1/4;
f(:,1) = prod(1-cos(x*pi/2), 2);
for i = 2:mMinusOne
    f(:,i) = prod(1-cos(x(:,1:M-i)*pi/2), 2) .* (1-sin(x(:,M-i+1)*pi/2).^wRate);
end
f(:,M) = (1-sin(x(:,1)*pi/2).^wRate);
end



% Mixed.
function f = h_mixed(x,alpha,A)

f = (1 - x(:,1) - cos(2*A*pi*x(:,1) + pi/2) / (2*A*pi)).^alpha;

end

% Disconnected.
function f = h_disc(x,alpha,beta,A)
f = 1 - x(:,1).^alpha .* cos(A * x(:,1).^beta * pi).^2;
end