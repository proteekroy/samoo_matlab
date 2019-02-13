
%cluster in the decision space
function [pop_cluster] = niching(opt, pop)

    archive = normalize(opt, opt.archive, opt.bound(1,:), opt.bound(2,:));
    normalized_pop = normalize(opt, pop, opt.bound(1,:), opt.bound(2,:));
    
    closest_point_index = dsearchn(archive,normalized_pop);
    pop_cluster = opt.archiveCluster(closest_point_index);
    

end