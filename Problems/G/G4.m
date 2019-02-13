function varargout = G4(Operation, Global, input)
% Proteek Chandan Roy 2018, email:proteek_buet@yahoo.com

    switch Operation
        case 'init'
            Global.M        = 1;
            Global.D        = 5;
            Global.C        = 6;
            
            bound = zeros(2, Global.D);
            bound(1,1)  = 78;
            bound(2,1)  = 102;
            bound(1,2)  = 33;
            bound(2,2)  = 45;
            bound(1,3:5)  = 27;
            bound(2,3:5)  = 45;
                    
            Global.lower    = bound(1,:);
            Global.upper    = bound(2,:);
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            x = input;
            PopObj = 5.3578547*x(:,3).^2+0.8356891*x(:,1).*x(:,5)+37.293239.*x(:,1)-40792.141;

            % Constraints
            u = 85.334407+0.0056858*x(:,2).*x(:,5)+0.0006262*x(:,1).*x(:,4)-0.0022053*x(:,3).*x(:,5);
            y(:,1) = -u;
            y(:,2) = u-92;
            v = 80.51249+0.0071317*x(:,2).*x(:,5)+0.0029955*x(:,1).*x(:,2)+0.0021813*x(:,3).^2;
            y(:,3) = -v+90;
            y(:,4) = v-110;
            w = 9.300961+0.0047026*x(:,3).*x(:,5)+0.0012547*x(:,1).*x(:,3)+0.0019085.*x(:,3).*x(:,4);
            y(:,5) = -w+20;
            y(:,6) = w-25; 
            PopCon = y; 
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = -3.066553867178332e+004;
            varargout = {f};
        case 'PS'
            varargout  = {[ 78, 33, 29.9952560256815985, 45, 36.7758129057882073 ]}; 
    end
end