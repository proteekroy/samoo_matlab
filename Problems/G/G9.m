function varargout = G9(Operation, Global, input)
% Proteek Chandan Roy 2018, email:proteek_buet@yahoo.com

    switch Operation
        case 'init'
            Global.M        = 1;
            Global.D        = 7;
            Global.C        = 4;
            
            bound = zeros(2, Global.D);
            bound(1,1:Global.D) = -10;
            bound(2,1:Global.D) = 10;
                    
            Global.lower    = bound(1,:);
            Global.upper    = bound(2,:);
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            x = input;
            PopObj = (x(:,1)-10).^2+5*(x(:,2)-12).^2+x(:,3).^4+3*(x(:,4)-11).^2+...
                10*x(:,5).^6+7*x(:,6).^2+x(:,7).^4-4*x(:,6).*x(:,7)-10*x(:,6)-8*x(:,7);

            % Constraints
            v1 = 2*x(:,1).^2;
            v2 = x(:,2).^2;
            y(:,1) = v1+3*v2.^2+x(:,3)+4*x(:,4).^2+5*x(:,5)-127;
            y(:,2) = 7*x(:,1)+3*x(:,2)+10*x(:,3).^2+x(:,4)-x(:,5)-282;
            y(:,3) = 23*x(:,1)+v2+6*x(:,6).^2-8*x(:,7)-196;
            y(:,4) = 2*v1+v2-3*x(:,1).*x(:,2)+2*x(:,3).^2+5.*x(:,6)-11*x(:,7);
            
            PopCon = y; 
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = 680.630057374402;
            varargout = {f};
        case 'PS'
            varargout  = {[2.33049935147405174, 1.95137236847114592, -0.477541399510615805, 4.36572624923625874, -0.624486959100388983, 1.03813099410962173, 1.5942266780671519]}; 
    end
end