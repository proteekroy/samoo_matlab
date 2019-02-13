function I = choose_active_set_TrustRegion(opt)


    %I1 = find(opt.div<opt.delta);%all the elements which has
    nd_set = opt.archive(opt.ParetoIndex,:);
    nd_size = size(nd_set,1);
    if nd_size~=0
        temp_n = floor(opt.activeSetSize/size(nd_set,1));%number of solution per non-dominated point
        
        if(temp_n<opt.SampleSize)
            temp_n = opt.SampleSize;
        end
        
        idx = knnsearch(opt.archive, nd_set, 'K', temp_n);%for each element in nd_set, find      
        idx = reshape(idx, [1, temp_n*nd_size]);%make one row
        I = unique(idx);%find unique indices
        %color = [rand rand rand];%RGB  color to show selection
        %hold all;
        %plot_points(opt);
        %plot(opt.archive(I,1), opt.archive(I,2),'*','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize', 5);
        %disp('hello');
    else
        %find the nearest
        I = randsample(opt.activeSetSize, size(opt.archive,1));
        %disp('hello');
        
    end
    

end