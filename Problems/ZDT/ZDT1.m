function varargout = ZDT1(Operation,Global,input)

    switch Operation
        case 'init'
            Global.M=2;
            Global.D=10;
            Global.C=0;
            Global.lower=zeros(1,Global.D);
            Global.upper=ones(1,Global.D);
            Global.operator=@EAreal;
            varargout = {Global};
        case 'value'
            PopDec = input;
            
            PopObj(:,1) = PopDec(:,1);
            g = 1+9*sum(PopDec(:,2:end),2)./(size(input,2)-1);
            h = 1-(PopObj(:,1)./g).^0.5;
            PopObj(:,2) = g.*h;
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f(:,1)    = (0:1/(input-1):1)';
            f(:,2)    = 1-f(:,1).^0.5;
            varargout = {f};
    end
end