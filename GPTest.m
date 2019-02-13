
%testing java code compare to matlab DACE fit 

function GPTest()

    clear java;
    javaaddpath(fullfile(pwd,'gp.jar'));
    %javaaddpath(pwd);
    %import gp.*;
    obj = gp.GaussianProcessJava;
    %create points
    n = 10;
    size = 200;
    theta = 10;% starting value of Kriging correlation parameter
    lob = 1e-2;%lower bound of Kriging correlation parameter
    upb = 20;%20;%upper bound of Kriging correlation parameter
    
    pop = rand(size,n);
    
    opt.objfunction = 'zdt1';
    opt.M = 2;
    opt.C = 0;
    [popObj, ~] = evaluate_pop(opt, pop);
    f1 = popObj(:,1);
    f2 = popObj(:,2);
    
    GaussianProcessJava.trainGP(pop, f1, 1);
    
    
    [model1, ~] = dacefit(pop, f1, @regpoly1, @corrgauss, theta, lob, upb);
    [model2, ~] = dacefit(pop, f2, @regpoly1, @corrgauss, theta, lob, upb);
        

    pop2 = rand(size,n);
    y_predict = zeros(size,opt.M);
    
    for i=1:size
        [y_predict(i,1),~,mse_asf,~] = predictor(pop2(i,:), model1);
        [y_predict(i,2),~,mse_asf,~] = predictor(pop2(i,:), model2);
    end
    
    


end