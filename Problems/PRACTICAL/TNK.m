function varargout = TNK(Operation,Global,input)

    switch Operation
        case 'init'
            Global.M        = 2;
            Global.D        = 2;
            Global.C        = 2;
            
            bound = zeros(2,Global.D);
            bound(2,:) = ones(1, Global.D);
            bound(1,1:end) = bound(1,1:end)+1e-12;
            bound(2,1:end) = bound(2,1:end)*pi;
            
            Global.lower    = bound(1,:);
            Global.upper    = bound(2,:);
            Global.operator = @EAreal;
            
            varargout = {Global};
        case 'value'
            x = input;
            
            PopObj(:, 1) = x(:,1);
            PopObj(:, 2) = x(:,2);
            PopCon(:,1) = (-1)*(x(:,1).^2+x(:,2).^2-1.0-0.1*cos(16*atan(x(:,1)./x(:,2))));
            PopCon(:,2) = (1/0.5)*((x(:,1)-0.5).^2+(x(:,2)-0.5).^2-0.5);            
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = load('TNK.2D.pf');
            varargout = {f};
    end
end