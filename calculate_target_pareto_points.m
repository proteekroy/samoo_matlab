
function point = calculate_target_pareto_points(filename, dirs)

    %opt.dirs = initweight(2, 21)';
    %pareto_front = load('GD/TNK.2D.jmetal.pf');
    pareto_front = load(filename);

    point = zeros(21,2);
    %hold all;
    for i = 1:21
        w = dirs(i,:);
        w = w./norm(w);
        %plot([0 w(1)*2],[0 w(2)*2],'-k');
        dist = intmax;
        index = 1;
        for j = 1:size(pareto_front,1)
            obj = pareto_front(j,:);
            obj = obj./norm(obj);
            d = norm(obj-((w*obj')*w)/(norm(w)*norm(w)));
            if d<dist
                index = j;
                dist = d;
            end

            %point(i,:) = pareto_front(index,:)./norm(pareto_front(index,:));
            point(i,:) = pareto_front(index,:);
        end

    end


    %plot(point(:,1), point(:,2),'ro', 'MarkerFaceColor','r', 'MarkerSize', 5);
    %plot(pareto_front(:,1), pareto_front(:,2),'ro', 'MarkerFaceColor','r', 'MarkerSize', 5);
end
