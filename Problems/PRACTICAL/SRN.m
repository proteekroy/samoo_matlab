function varargout = SRN(Operation,Global,input)

    switch Operation
        case 'init'
            Global.M        = 2;
            Global.D        = 2;
            Global.C        = 2;
            
            bound = zeros(2,Global.D);
            bound(2,:) = ones(1, Global.D);
            bound(1,1:end) = bound(1,1:end)+(-20);
            bound(2,1:end) = bound(2,1:end)*20;
            
            Global.lower    = bound(1,:);
            Global.upper    = bound(2,:);
            Global.operator = @EAreal;
            
            varargout = {Global};
        case 'value'
            x = input;
            
            
            PopObj(:, 1) = 2.0+(x(:,1)-2.0).^2+(x(:,2)-1.0).^2;
            PopObj(:, 2) = 9.0*x(:,1)-(x(:,2)-1).^2;
            PopCon(:, 1) = (1/225)*(x(:,1).^2+x(:,2).^2-225.0);
            PopCon(:, 2) = 0.1*(x(:,1)-3.0*x(:,2)+10.0);
            
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = load('SRN.2D.pf');
            varargout = {f};
    end
end