function sampling_test()%, opt)%n is the number of sampling

%%{
n = 500;
opt.objfunction = 'zdt6';
opt.methodology = 3;
opt = basic_parameters(opt);

opt.metamodelOption = 2;
opt.diffusionOption = 1;


if n<10
    pop = lhsamp_model(n, opt);
    pop_obj = evaluate_pop(opt, pop);
    return;
end


random_pop = zeros(n,opt.V);
random_pop_obj = zeros(n,opt.M);

p = ceil(n/10);

K = 2*opt.V;
random_pop(1:p,:) = lhsamp_model(p, opt);%create 25% random solutions

if(opt.metamodelOption==1)
    random_pop_obj(1:p,:) = evaluate_pop_metamodel(opt, random_pop(1:p,:));%evaluate by metamodel
else
    random_pop_obj(1:p,:) = evaluate_pop(opt, random_pop(1:p,:));
end

figure
hold all;
plot(random_pop_obj(1:p,1),random_pop_obj(1:p,2),'sb');

pareto = load('../Metamodel Classification/ZDT/ZDT6.pf');
plot(pareto(:,1),pareto(:,2),'k.'); 

for i = p+1:2:n
    
    %{
    IDX = knnsearch(value,value,'k',k+1,'distance','euclidean');%k-d tree, this is faster if [n,d]=knnsearch(x,y,'k',10,'distance','minkowski','p',5);
    
    [D,~] = pdist2(random_pop_obj(1:i-1,:), random_pop_obj(1:i-1,:),'euclidean','Smallest',K);%K nearest to calculate the distance
    D = D(2:end,:);
    d = mean(D,1)';% average distance of K nearest neighbor
    %}
    
    d = diffusion_diversity_metric(opt, random_pop_obj(1:i-1,:), 5);
    d = exp(5*d);
    normConst = sum(d);%normalization constant

    d = d./normConst;%probability for choosing the solution, will be used for crossover and mutation

    cd = cumsum(d);%cumulative distribution for roulette wheel selection
    
    %if(rand<0.5)
        pop_child = sbx_sampling( random_pop(1:i-1,:), 2,  1.0, cd, opt.eta_c, opt.bound(1,:), opt.bound(2,:));%produce child by crossover
        
    %else
    %    pop_child = pol_mut_sampling(random_pop(1:i-1,:), 2, 1.0, cd, opt.eta_m,  opt.bound(1,:), opt.bound(2,:) );%produce child by mutation
    %end
    
    %disp(pop_child');
    random_pop(i:i+1,:) = pop_child;%create 25% random solutions

    %if(opt.metamodelOption==1)
    %    random_pop_obj(i:i+1,:) = evaluate_pop_metamodel(opt, pop_child);%evaluate by metamodel
    %else
        random_pop_obj(i:i+1,:) = evaluate_pop(opt, pop_child);
    %end

    plot(random_pop_obj(i:i+1,1),random_pop_obj(i:i+1,2),'sr','MarkerFaceColor','r');
    
end


if(i<n)
    p = n - i;
    random_pop(i+1:n,:) = rand(p, opt.V);
    random_pop_obj(i+1:n,:) = evaluate_pop(opt, random_pop(i+1:n,:));
    plot(random_pop_obj(i+1:n,1),random_pop_obj(i+1:n,2),'go','MarkerFaceColor','g');
end




