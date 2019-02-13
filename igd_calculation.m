function val = igd_calculation(ParetoObj, FeasibleArchiveObj)

    dist = pdist2(ParetoObj, FeasibleArchiveObj,'Euclidean');%distance in normalized space
    
    dist = min(dist,[],2);
    
    val = mean(dist);


end