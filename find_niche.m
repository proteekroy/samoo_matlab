function [cluster] = find_niche(opt, pop, popArchive)

    
    
    normalizedpop = normalize(opt, pop, opt.bound(1,:), opt.bound(2,:));%normalize decision space
    
    [Idx,~] = knnsearch(popArchive, normalizedpop,'k',1);%find nearest solution in decision space
    
    cluster = opt.archiveCluster(Idx);%cluster is assigned by the cluster of nearest neighbor
    
    
end