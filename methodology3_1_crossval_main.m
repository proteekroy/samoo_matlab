% Copyright [2018] [Proteek Chandan Roy]
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

%
% Proteek Chandan Roy, 2018
% Contact: proteek_buet@yahoo.com

%This file performs k-fold cross valiaditon for Framework 3-1

function  [returnval] = methodology3_1_crossval_main(opt, TrainIndex, TestIndex, k)

    returnval = cell(1,2);
    sep = zeros(1,10);
    yerror = zeros(1,10);
    
    for i=1:k
        [sep(i), yerror(i)] = methodology3_1_crossval(opt, opt.archive(TrainIndex{i},:), opt.archive(TestIndex{i},:), ...
            opt.archiveASFAll, TrainIndex{i}, TestIndex{i}, opt.archiveCons(TrainIndex{i},:), opt.archiveCons(TestIndex{i},:));
    end
    
    returnval{1} = sep; 
    returnval{2} = yerror;
    
end


function [sep, yerror] =  methodology3_1_crossval(opt, xtrain, xtest, ASF, TrainIndex, TestIndex, cons_train, cons_test)

    [~,ia,~] = unique(xtrain,'rows');
    xtrain = xtrain(ia,:);
    cons_train = cons_train(ia,:);
    
    y_asf = cell2mat(ASF);
    y_asf_train = y_asf(TrainIndex,:);
    y_asf_train = y_asf_train(ia,:);
    y_asf_test = y_asf(TestIndex,:);

    cons_error = zeros(size(xtest,1), 1);
    
    [original_cv, ~] = evaluateCV(cons_test);
    index1 = original_cv<=0;
    
    %----------------MODEL CONSTRAINTS-------------------------------------
    if opt.C>0
        [opt.dmodel_cons, ~] = dacefit(xtrain, cons_train, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model cv
    end
    %--------------------PREDICT CONSTRAINTS-------------------------------
    n = size(xtest,1);
    if opt.C>0
        [cons_pred,~] = predictor(xtest, opt.dmodel_cons);
    end
    
    if opt.C>0
        [test_cv, ~] = evaluateCV(cons_test);
        [pred_cv, ~] = evaluateCV(cons_pred);
    else
        test_cv = zeros(n, 1);
        pred_cv = zeros(n, 1);
    end
        
    if opt.C>0
        index2 = pred_cv<=0;
        index = index1 & index2;
        if sum(index)>0
            cons_pred(index,:) = zeros(sum(index),opt.C);   
        end
        cons_error = sqrt(sum((cons_pred-cons_test).^2,2));
    end
    
    %----------------COMPUTE CONSTRAINTS ERROR-----------------------------
    
    selection_err = 0;
    counter = 0;
    n = size(xtest,1);
    
    active_size = min(opt.activeSetSize, size(xtrain, 1));
    rnd = randperm(size(xtrain,1));
    rnd = rnd(1:active_size);
    
    [opt.dmodel_asf, ~] = dacefit(xtrain(rnd,:), y_asf_train(rnd,:), @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model asf
    [y_pred,~] = predictor(xtest, opt.dmodel_asf);
    
    
    for dir = 1:1: size(opt.dirs,1) 
        
        %---------------------COMPUTE SEP----------------------------------     
        for i=1:n-1
            for j=i+1:n
                d1 = constrained_domination(y_asf_test(i, dir), y_asf_test(j, dir), test_cv(i), test_cv(j));%1, 2, 3
                d2 = constrained_domination(y_pred(i,dir), y_pred(j,dir), pred_cv(i), pred_cv(j));
                if d1~=d2
                    selection_err = selection_err + 1;
                end
                counter = counter + 1;
            end
        end

    end
    %----------------COMPUTE ASF ERROR---------------------------------
    obj_error = sqrt(sum(((y_pred-y_asf_test)./y_asf_test).^2,2));
    
    sep = selection_err/counter;
    
    obj_error  = mean(obj_error);
    cons_error = mean(cons_error);
    yerror = sqrt(obj_error.^2 + cons_error.^2);

end


