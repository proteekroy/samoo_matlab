        
function opt = buildModel(opt)
    
    [~,ia,~] = unique(opt.nsga2.hifipop,'rows');
    temp.hifipop = opt.nsga2.hifipop(ia,:);
    temp.hifipopObj = opt.nsga2.hifipopObj(ia,:);
    temp.hifipopCons = opt.nsga2.hifipopCons(ia,:);
    temp.hifipopCV = opt.nsga2.hifipopCV(ia, :);
    %%{
    for i=1:opt.M
        [opt.dmodel_obj{i}, ~] = dacefit(temp.hifipop, temp.hifipopObj(:,i), @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model asf
    end
    
    if opt.C>0
        switch(opt.methodology)
            
            case {12} 
                for i=1:opt.C
                    [opt.dmodel_cons{i}, ~] = dacefit(temp.hifipop, temp.hifipopCons(:,i), @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model cv
                end
            case {22}
                [opt.dmodel_cv, ~] = dacefit(temp.hifipop, temp.hifipopCV, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model cv
                
        end
    end
    %}   
    %{
    opt.regC{1} = quadratic_model(temp.hifipop, temp.hifipopCons(:,1),1);%model constraints
    opt.regC{2} = quadratic_model(temp.hifipop, temp.hifipopCons(:,2),1);%model constraints
    opt.regC{3} = quadratic_model(temp.hifipop, temp.hifipopCons(:,3),1);%model constraints
    opt.regC{4} = quadratic_model(temp.hifipop, temp.hifipopCons(:,4),1);%model constraints
    opt.regC{5} = quadratic_model(temp.hifipop, temp.hifipopCons(:,5),2);%model constraints
    opt.regC{6} = quadratic_model(temp.hifipop, temp.hifipopCons(:,6),2);%model constraints
        
    %}
%     for i=1:opt.M
%         [opt.dmodel_asf{i}, ~] = dacefit(opt.hifipop, opt.hifipopObj(:,i), @regpoly1, @corrgauss, opt.theta, opt.lob, opt.upb);%model asf
%     end
    
    %{
    opt.regO{1} = quadratic_model(temp.hifipop, temp.hifipopObj(:,1),2);%model objectives
    opt.regO{2} = quadratic_model(temp.hifipop, temp.hifipopObj(:,2),2);%model objectives
    %}
    %{
    [opt.dmodel_cons, ~] = dacefit(opt.hifipop, opt.hifipopCons(:,1), @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model cv
    opt.regC{1} = quadratic_model(opt.hifipop, opt.hifipopCons(:,2),2);%model constraints
    %}
end