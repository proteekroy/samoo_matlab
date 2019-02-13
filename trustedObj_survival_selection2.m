function [selected_pop, selected_pop_obj, selected_pop_cv] = trustedObj_survival_selection2(opt, pop, popObj, popCV)

    num = opt.N;
    %--------------RETURN VALUES-------------------------------------------
    selected_size = 0;
    selected_pop = zeros(opt.N, opt.V);
    selected_pop_obj = zeros(opt.N, size(popObj,2));
    selected_pop_cv = zeros(opt.N, size(popCV,2));
    
    %------------FIND WHO IS WITHIN TRUST REGION---------------------------
    
    %{
    %leader = opt.archive(opt.ParetoIndex,:);
    leader = opt.archive(opt.LeaderIndex{opt.curcluster},:);
    %normalize leaeders and population in x-space with bound
    leader = normalize(opt, leader, opt.bound(1,:), opt.bound(2,:));
    norm_pop = normalize(opt, pop, opt.bound(1,:), opt.bound(2,:)); 
    div = pdist2(norm_pop, leader);
    div = min(div,[],2);
    %}
    
    %%{
    leader = normalize(opt, opt.archive(opt.archiveCluster==opt.curcluster), opt.bound(1,:), opt.bound(2,:));
    if isempty(leader)
       leader =  normalize(opt, opt.archive, opt.bound(1,:), opt.bound(2,:));
    end
    norm_pop = normalize(opt, pop, opt.bound(1,:), opt.bound(2,:)); 
    div = pdist2(norm_pop, leader);%distance in normalized space
    div = min(div,[],2);
    I1 = find(div<=opt.TrustDistVar);
    I2 = find(div>opt.TrustDistVar);
    %}
    %{
    leader = normalize(opt, opt.activeArchiveObj, opt.min_val, opt.max_val);
        
    norm_pop = normalize(opt, popObj, opt.min_val, opt.max_val); 
    div = pdist2(norm_pop, leader);%distance in normalized space
    div = min(div,[],2);
    I1 = find(div<=opt.TrustDistObj);
    I2 = find(div>opt.TrustDistObj);
    %}
    
    
    if size(I1,1)>opt.N
        pop = pop(I1,:);
        popObj = popObj(I1, :);
        popCV = popCV(I1,:);
    else
        selected_pop(1:size(I1,1), :) = pop(I1,1:opt.V); 
        selected_pop_obj(1:size(I1,1),:) = popObj(I1,:);
        selected_pop_cv(1:size(I1,1),:) = popCV(I1,:);
        selected_size = size(I1,1);
        pop = pop(I2,:);
        popObj = popObj(I2,:);
        popCV = popCV(I2,:);
    end
    
    
    %-------------FIND CLUSTER NUMBERS-------------------------------------
    
    archive = normalize(opt, opt.archive, opt.bound(1,:), opt.bound(2,:));
    normalized_pop = normalize(opt, pop, opt.bound(1,:), opt.bound(2,:));
    
    closest_point_index = dsearchn(archive,normalized_pop);
    pop_cluster = opt.archiveCluster(closest_point_index);
    
    %---------------COLLECT CLUSTER MEMBERS AND SORT THEM------------------
    index = cell(opt.numdir,1);
    
    for i = 1:opt.numdir %number of clusters
        
        index{i} = find(pop_cluster == i); %all solutions in cluster i
        if ~isempty(index{i})
            obj = popObj(index{i},:); %objectives of cluster i
            [~,I] = sort(obj(:,1),'ascend'); %should be non-dominated sort for multiple objectives
            index{i} = index{i}(I,:);%store in a sorted order    
        end
    end

    
    %---------SELECTION BASED ON CLUSTERING AND MIN ASF VALUE--------------
    k = selected_size+1;
    while(k<=num)
        
        for j=1:opt.numdir
            if ~isempty(index{j})
                selected_pop(k, :) = pop(index{j}(1),1:opt.V); 
                selected_pop_obj(k,:) = popObj(index{j}(1),:);
                selected_pop_cv(k,:) = popCV(index{j}(1),:);
                index{j}(1) = [];%delete from the cluster as taken
                k = k + 1;
            %else
            %    disp('Infinite Loop');
            end
            if(k > num)
                break;
            end
        end
    end

end

%---------------------------END OF FILE------------------------------------