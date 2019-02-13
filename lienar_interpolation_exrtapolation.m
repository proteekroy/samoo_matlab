

function [pop_interpolation, pop_extrapolation] =  lienar_interpolation_exrtapolation(X, gen)
      
    
    [n, ~] = size(X);
    neighbor_size = min(5, n);%consider at least 5 edges
    
    [idx, D] = knnsearch(X, X, 'K', neighbor_size);
    idx(:, 1) = []; %delete the first column
    D(:, 1) = [];% detele firt column distances
    s = [];
    for i=1:size(idx, 1)
        s = horzcat(s, repmat(i,1,neighbor_size-1));
    end
    
    %build minimum spanning tree
    t = reshape(idx', [1, size(idx,1)*size(idx,2)]);
    D = reshape(D', [1, size(D,1)*size(D,2)]);
    
    SampleEdgeList = vertcat(s,t)';
    [uniqueEdgeList,ia,~] = unique(sort(SampleEdgeList,2),'rows');
    s = uniqueEdgeList(:,1);
    t = uniqueEdgeList(:,2);
    weights = D(ia);
    G = graph(s,t,weights);
    [T, ~] = minspantree(G);
    
    %collect interpolating and extrapolating points with pre-defined ratios
    pop_interpolation = [];
    pop_extrapolation = [];
    ratio1 = gen;
    ratio2 = 1;
    
    
    fig_flag = 2;
    if fig_flag==1
        figure;
        hold on;
        box on;
        grid on;
    end
    
    for i=1:size(T.Edges.EndNodes,1) %for each edge
        s = T.Edges.EndNodes(i, 1);
        t = T.Edges.EndNodes(i, 2);
        point_interpolate = (X(s,:)+X(t,:))./2;
        point_extrapolate = (ratio1*X(s,:)-ratio2*X(t,:))./(ratio1-ratio2);
        pop_interpolation = vertcat(pop_interpolation, point_interpolate);
        pop_extrapolation = vertcat(pop_extrapolation, point_extrapolate);
        if fig_flag==1
            plot([X(s,1), X(t,1)],[X(s,2), X(t,2)],'b-','LineWidth',2);%[1.0000 0.6725 0.2627]);
            plot([X(s,1), point_extrapolate(1)],[X(s,2), point_extrapolate(2)],'r-','LineWidth',2);%[1.0000 0.6725 0.2627]);
        end
        
    end
    
    %figure
    if fig_flag==1
        plot(X(:,1),X(:,2),'bo','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',10);%[1.0000 0.6725 0.2627]);
        plot(pop_interpolation(:,1), pop_interpolation(:,2),'o','MarkerFaceColor',[1.0000 0.6725 0.2627], 'MarkerSize', 7);
        plot(pop_extrapolation(:,1), pop_extrapolation(:,2),'go','MarkerFaceColor','g', 'MarkerSize', 7);
        hold off;
        %figure;
        %p = plot(G,'EdgeLabel',G.Edges.Weight);
        %highlight(p,T)
    end
end