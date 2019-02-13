function varargout = OSY(Operation,Global,input)

    switch Operation
        case 'init'
            Global.M        = 2;
            Global.D        = 6;
            Global.C        = 6;
            
            bound = zeros(2,Global.D);
            bound(2,:) = ones(1, Global.D);
            bound(2,1)=10;
            bound(2,2)=10;
            bound(2,6)=10;
            bound(1,3)=1;
            bound(1,5)=1;
            bound(2,3)=5;
            bound(2,5)=5;
            bound(2,4)=6;
            
            Global.lower    = bound(1,:);
            Global.upper    = bound(2,:);
            Global.operator = @EAreal;
            varargout = {Global};
        case 'value'
            x = input;
            
            PopObj(:, 1) = -(25*(x(:,1)-2).^2+(x(:,2)-2).^2+(x(:,3)-1).^2+(x(:,4)-4).^2+(x(:,5)-1).^2);
            PopObj(:, 2) = x(:,1).^2+x(:,2).^2+x(:,3).^2+x(:,4).^2+x(:,5).^2+x(:,6).^2;
            PopCon(:,1) = (-1/2)*(x(:,1)+x(:,2)-2);
            PopCon(:,2) = (-1/6)*(6-x(:,1)-x(:,2));
            PopCon(:,3) = (-1/2)*(2-x(:,2)+x(:,1));
            PopCon(:,4) = (-1/2)*(2-x(:,1)+3*x(:,2));
            PopCon(:,5) = (-1/4)*(4-(x(:,3)-3).^2 -x(:,4));
            PopCon(:,6) = (-1/4)*((x(:,5)-3).^2+x(:,6)-4);
            
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = load('OSY.2D.pf');
            varargout = {f};
    end
end