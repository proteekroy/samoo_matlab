function pop = kktpm_combine_and_reduce_model(opt, pop1, pop2)

[N, k] = size(pop1); k = k-2; % N = 2 * N_initial


pop2 = [pop2, zeros(N, 2)];
pop2 = evaluate_fitness_model(opt, pop2(:, 1:k));%FIND KKTPM PREDICTION OF CHILDREN POPULATION
pop3 = [pop1; pop2];


%FIND CLUSTER CENTERS
closest_point_index = dsearchn(opt.activearchive(:,1:opt.V),pop3(:,1:opt.V));%for each element in pop, find closest solution in x-space
pop_cluster = opt.activearchiveCluster(closest_point_index,:); %cluster values of those solutions
    





    % START THE COMPETITION
    pop = [];%zeros(N, k+2);
    
    
    
    temp = zeros(opt.numdir,1);%saves number of elements for each direction
    member = (1:opt.numdir)';
    for i=1:2*N
        temp(pop_cluster(i)) =  temp(pop_cluster(i))+1;%find number of members for each cluster
    end
    
    [~,index] = sort(temp);
    temp = temp(index,:);%total size
    member = member(index,:);%cluster number
    
    s=1;
    remain = N;
    for i=1:opt.numdir
        if (opt.numdir-i+1)>remain
            h = input('hello');
        end
        t = floor(remain/(opt.numdir-i+1));
        if temp(i)>0 && temp(i)<=t %size of the
            r = 0;
            for j=1:2*N
               if  pop_cluster(j)==member(i)
                   pop(s,1:opt.V+2) = pop3(j,1:opt.V+2);
                   s = s+1;
                   r = r+1;
               end
            end
            remain = remain -r; 
        elseif temp(i)>t %find smallest kktpm values
            collect = []; 
            kktpm = [];
            q = 1;
            for j=1:2*N
                if  pop_cluster(j)==member(i)
                   collect(q,1:opt.V+2) = pop3(j,1:opt.V+2);
                   kktpm(q,1) =  pop3(j,opt.V+1);
                   q = q+1;
               end
            end
            
            [~,index] = sort(kktpm);
            if i== opt.numdir
                collect = collect(index(1:remain),:);
            else
                collect = collect(index(1:t),:);
            end
            pop(s:s+t-1,1:opt.V+2)  = collect;
            s = s + size(collect,1);
            remain = remain - t; 
        end
    
    end
    
    if size(pop,1)<N
        disp('here');
    end
end


