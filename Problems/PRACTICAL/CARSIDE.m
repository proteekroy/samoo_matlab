function varargout = CARSIDE(Operation,Global,input)

    switch Operation
        case 'init'
            Global.M        = 3;
            Global.D        = 7;
            Global.C        = 10;
            
            bound = zeros(2, Global.D);
            bound(1,1:end) = [0.5 0.45 0.5 0.5 0.875 0.4 0.4];
            bound(2,1:end) = [1.5 1.35 1.5 1.5 2.625 1.2 1.2];
            
            Global.lower    = bound(1,:);
            Global.upper    = bound(2,:);
            Global.operator = @EAreal;
            
            varargout = {Global};
        case 'value'
            x = input;
            
            PopObj(:, 1) = 4*x(:,1).^2+4*x(:,2).^2;
            PopObj(:, 2) = (x(:,1)-5).^2+(x(:,2)-5).^2;
            PopCon(:,1) = (1/25)*((x(:,1)-5).^2 + x(:,2).^2-25);   
            PopCon(:,2) = -1/7.7*((x(:,1)-8).^2 + (x(:,2)+3).^2-7.7);            
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = load('BNH.2D.pf');
            varargout = {f};
    end
end