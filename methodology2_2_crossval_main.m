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

%This file performs k-fold cross valiaditon for Framework 2-2

function  [returnval] = methodology2_2_crossval_main(opt, TrainIndex, TestIndex, k)

    returnval = cell(1,2);
    sep = zeros(1,10);
    yerror = zeros(1,10);
    
    for i=1:k
        [sep(i), yerror(i)] = methodology2_2_crossval(opt, opt.archive(TrainIndex{i},:), opt.archive(TestIndex{i},:), ...
            opt.archiveObj(TrainIndex{i},:), opt.archiveObj(TestIndex{i},:), opt.archiveCVModified(TrainIndex{i},:), opt.archiveCVModified(TestIndex{i},:));
    end

    returnval{1} = sep; 
    returnval{2} = yerror;
end


function [sep, yerror] =  methodology2_2_crossval(opt, xtrain, xtest, train_obj, test_obj, cons_train, cons_test)

    [~,ia,~] = unique(xtrain,'rows');
    xtrain = xtrain(ia,:);
    train_obj = train_obj(ia,:);
    cons_train = cons_train(ia,:);
    

    [original_cv, ~] = evaluateCV(cons_test);
    index1 = original_cv<=0;
        
    %------------------MODEL OBJECTIVES------------------------------------
    [opt.dmodel_obj, ~] = dacefit(xtrain, train_obj, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model asf
    
        
    %----------------MODEL CONSTRAINTS-------------------------------------
    if opt.C>0
        [opt.dmodel_cv, ~] = dacefit(xtrain, cons_train, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model cv
    end
       
    %---------------------PREDICT OBJECTIVES-------------------------------
    [y_pred,~] = predictor(xtest, opt.dmodel_obj);
         
    %--------------------PREDICT CONSTRAINTS-------------------------------
    n = size(xtest,1);
    if opt.C>0
        [cons_pred,~] = predictor(xtest, opt.dmodel_cv);
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
        index = index1 & index2; % solutions where both actual and predicted solutions are feasible
        if sum(index)>0 % when there are solutions where both are feasible, it occurs no error 
            cons_pred(index,:) = zeros(sum(index), 1);   
        end        
        cons_error = sqrt(sum(((cons_pred-cons_test)./cons_test).^2,2));
        cons_error = mean(cons_error); 
    else
        cons_error = 0;
    end
    
    %--------------------COMPUTE SEP---------------------------------------
    selection_err = 0;
    counter = 0;
    n = size(xtest,1);
    
    for i=1:n-1
        for j=i+1:n
            d1 = constrained_domination(test_obj(i,:), test_obj(j,:), test_cv(i,:), test_cv(j,:));%1, 2, 3
            d2 = constrained_domination(y_pred(i,:), y_pred(j,:), pred_cv(i,:), pred_cv(j,:));
            if d1~=d2
                selection_err = selection_err + 1;
            end
            counter = counter + 1;
        end
    end
    sep = selection_err/counter;
   
    %----------------COMPUTE MSE ERROR-------------------------------------
    obj_error = sqrt(sum(((y_pred-test_obj)./test_obj).^2,2));
    obj_error  = mean(obj_error);
    yerror = sqrt(obj_error.^2 + cons_error.^2);
    
end


