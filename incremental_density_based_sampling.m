% Copyright [2016] [Proteek Chandan Roy]
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

%
% Proteek Chandan Roy, 2016
% Contact: royprote@msu.edu, proteek_buet@yahoo.com

%n = the number of sampling
%opt = option for sampling space

function [pop, popObj, popCons] = incremental_density_based_sampling(n, opt)

%-------IF FUNCTION EVALUATION SMALL, USE LHS------------------------------
%{
n = 200;
opt.objfunction = 'zdt6';
opt.methodology = 7;
opt.r = 1;
opt = basic_parameters(opt);
%}
if n<10
    pop = lhsamp_model(n, opt);
    [popObj, popCons] = evaluate_pop(opt, pop);
    return;
end

%----------ELSE USE INCREMENTAL DENSITY BASED SAMPLING---------------------

random_pop = zeros(n,opt.V);
random_pop_obj = zeros(n, opt.M);
if opt.C>0
    random_pop_cons = zeros(n, opt.C);
else
    random_pop_cons = zeros(n, 1);
end

p = floor(n/10);
K = 5;%2*opt.V;%K nearest neighbor for density calculation
random_pop(1:p,:) = lhsamp_model(p, opt);%create 25% random solutions
[random_pop_obj(1:p,:), random_pop_cons(1:p,:)] = evaluate_pop(opt, random_pop(1:p,:));


%{
figure
hold all;
plot(random_pop_obj(1:p,1),random_pop_obj(1:p,2),'sb');
%}



for i = p+1:2:n
    
    
    %[D,~] = pdist2(random_pop_obj(1:i-1,:), random_pop_obj(1:i-1,:),'euclidean','Smallest',K);%K nearest to calculate the distance
    %D = D(2:end,:);
    [~,D] = knnsearch(random_pop_obj(1:i-1,:), random_pop_obj(1:i-1,:),'k', K);%find nearest solution in decision space
    D = D(:,2:end);
    d = mean(D,2)';% average distance of K nearest neighbor
    d = 3*d;
    %d = exp(5*d);
    normConst = sum(d);%normalization constant

    d = d./repmat(normConst,size(d,1),1);%probability for choosing the solution, will be used for crossover and mutation

    cd = cumsum(d);%cumulative distribution for roulette wheel selection
    
    %if(rand<0.5)
        %pop_child = sbx_sampling( random_pop(1:i-1,:), 2,  1.0, cd, opt.eta_c, opt.bound(1,:), opt.bound(2,:));%produce child by crossover
        pop_child = sbx_sampling( random_pop(1:i-1,:), 2,  1.0, cd, 30, opt.bound(1,:), opt.bound(2,:));%produce child by crossover
        %pop_child = pol_mut(pop_child, opt.pmut, opt.nrealmut, opt.eta_m,  opt.bound(1,:), opt.bound(2,:) );%produce child by mutation
        pop_child = pol_mut(pop_child, opt.pmut, opt.nrealmut, 30,  opt.bound(1,:), opt.bound(2,:) );%produce child by mutation
        
    %else
    %    pop_child = pol_mut_sampling(random_pop(1:i-1,:), 2, 1.0, cd, opt.eta_m,  opt.bound(1,:), opt.bound(2,:) );%produce child by mutation
    %end
    
    %disp(pop_child');
    random_pop(i:i+1,:) = pop_child;%create 25% random solutions

    %if(opt.metamodelOption==1)
    %    random_pop_obj(i:i+1,:) = evaluate_pop_metamodel(opt, pop_child);%evaluate by metamodel
    %else
        [random_pop_obj(i:i+1,:), random_pop_cons(i:i+1,:)] = evaluate_pop(opt, pop_child);
    %end

    %plot(random_pop_obj(i:i+1,1),random_pop_obj(i:i+1,2),'sr','MarkerFaceColor','r');
    
end


if(i<n)
    p = n - i;
    random_pop(i+1:n,:) = rand(p, opt.V);
    [random_pop_obj(i+1:n,:), random_pop_cons(i+1:n,:)] = evaluate_pop(opt, random_pop(i+1:n,:));
    %plot(random_pop_obj(i+1:n,1),random_pop_obj(i+1:n,2),'go','MarkerFaceColor','g');
end

pop = random_pop;
popObj = random_pop_obj;
popCons = random_pop_cons;

%{
figure

plot(popObj(:,1),popObj(:,2),'bo','MarkerFaceColor','b');
xlabel('Objective 1')
ylabel('Objective 2')
title('ZDT6 with Incremental Initialization')
set(gca,'fontsize',18);
figure
pop = lhsamp_model(n, opt);
[popObj, popCons] = evaluate_pop(opt, pop);
plot(popObj(:,1),popObj(:,2),'bo','MarkerFaceColor','b');
title('ZDT6 with LHS')
xlabel('Objective 1')
ylabel('Objective 2')
set(gca,'fontsize',18);
%}

end

