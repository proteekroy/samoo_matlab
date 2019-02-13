function varargout = G2(Operation, Global, input)
% Proteek Chandan Roy 2018, email:proteek_buet@yahoo.com

    switch Operation
        case 'init'
            Global.M        = 1;
            Global.D        = 20;
            Global.C        = 2;
            
            bound = zeros(2, Global.D);
            bound(2,1:Global.D) = 10*ones(1, Global.D);
                    
            Global.lower    = bound(1,:);
            Global.upper    = bound(2,:);
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            n = size(input, 2); 
            x = input;
            sum_jx = zeros(size(input,1),1);
            for j=1:n
                sum_jx = sum_jx+j*x(:,j).^2; 
            end
            PopObj = -abs((sum(cos(x).^4,2)-2*prod(cos(x).^2,2))./sqrt(sum_jx));

            % Constraints
            y(:,1) = -prod(x,2)+0.75;
            y(:,2) = sum(x,2)-7.5*n;      
            PopCon = y; 
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = -0.80361910412559;
            varargout = {f};
        case 'PS'
            varargout  = {[3.16246061572185, 3.12833142812967, 3.09479212988791, 3.06145059523469, 3.02792915885555, 2.99382606701730,...
                            2.95866871765285, 2.92184227312450, 0.49482511456933, 0.48835711005490, 0.48231642711865, 0.47664475092742,...
                            0.47129550835493, 0.46623099264167, 0.46142004984199, 0.45683664767217, 0.45245876903267, 0.44826762241853,...
                            0.44424700958760, 0.44038285956317]}; 
    end
end