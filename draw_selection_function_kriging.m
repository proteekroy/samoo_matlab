function draw_selection_function_kriging(dmodel)


    N = 50;
    z = [0 0];
    x1 = linspace(0,1,N);
    x2 = linspace(0,1,N);
    [X,Y] = meshgrid(x1,x2);
    Z = zeros(size(X));




for j=1:size(X,1)
    for k=1:size(X,2)
        Z(j,k) = predictor([X(j,k), Y(j,k), zeros(1,8)], dmodel);
    end
end

figure
hold all;
surfc(X,Y,Z);
% v = [0:2:30];
% contour(Z,v);
title('NSGA-II');
xlabel('x1')
ylabel('x2')
zlabel('Selection(x)')
set(gca,'fontsize',18);
drawnow;


end