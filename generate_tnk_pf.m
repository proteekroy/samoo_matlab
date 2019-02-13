function generate_tnk_pf()

N = 1000;
x1 = linspace(0,pi, N);
x2 = linspace(0,pi, N);
[x1,x2] = meshgrid(x1, x2);

x1 = reshape(x1, [N*N 1]);
x2 = reshape(x2, [N*N 1]);

X = horzcat(x1,x2);

g = zeros(N*N, 2);

g(:, 1) = (-1)*(X(:,1).^2+X(:,2).^2-1.0-0.1*cos(16*atan(X(:,1)./X(:,2))));
g(:, 2) = (1/0.5)*((X(:,1)-0.5).^2+(X(:,2)-0.5).^2-0.5);

for i=1:2
    g(g(:,i)<0,i)=0;
end
cv = sum(g, 2);

I = cv<=0;
F = X(I,:);


index = paretofront(F);
F = F(index,:);

pf_web = load('IGD Calculation/GD/TNK.2D.pf');
hold all;
plot(F(:,1), F(:,2), 'bo', 'MarkerFaceColor','b', 'MarkerSize', 5);
plot(pf_web(:,1), pf_web(:,2), 'ro', 'MarkerFaceColor','r', 'MarkerSize', 5);

F2 = vertcat(F,pf_web);
index = paretofront(F2);
F2 = F2(index,:);
plot(F2(:,1), F2(:,2), 'ko', 'MarkerFaceColor',[0.5451    0.2706    0.0745], 'MarkerSize', 5);
dlmwrite('IGD Calculation/GD/TNK.2D.new.pf', F, 'delimiter',' ','precision','%.10f','-append');

end