

%Find diverse set of solutions of size 'opt.NumExplorationPoints'

function selected_index = find_diverse_points(opt)


    temp = opt.div(opt.ParetoIndex);
    [~,ia] = sort(temp,'ascend');%sort by diversity value
    if size(temp,1)<opt.NumExplorationPoints
        selected_index = opt.ParetoIndex;
    else
        selected_index = opt.ParetoIndex(ia(1:opt.NumExplorationPoints));
    end
    %{
    hold all;
    currentPareto = opt.archiveObj(opt.ParetoIndex,:);
    plot(opt.archiveObj(:,1), opt.archiveObj(:,2), 'bo', 'MarkerSize', 5);
    plot(currentPareto(:,1), currentPareto(:,2),'mo','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize', 10);
    plot(opt.archive(selected_index, 1), opt.archive(selected_index, 2),'o','MarkerFaceColor',[.500 0.6725 0.2627], 'MarkerSize', 10);
    hold off;
    disp('hello');
    %}
end