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

%This file performs k-fold cross valiaditon for Framework 5

function  [returnval] = methodology5_crossval_main(opt, TrainIndex, TestIndex, k)

    returnval = cell(1,2);
    sep = zeros(1,10);
    yerror = zeros(1,10);
    
    for i=1:k
        [sep(i), yerror(i)] = methodology5_crossval(opt, opt.archive(TrainIndex{i},:), opt.archive(TestIndex{i},:), ...
            opt.archiveASFCVAll, TrainIndex{i}, TestIndex{i});
    end
    
    returnval{1} = sep; 
    returnval{2} = yerror;
    
end


function [sep, yerror] =  methodology5_crossval(opt, xtrain, xtest, ASFCV, TrainIndex, TestIndex)

    [~,ia,~] = unique(xtrain,'rows');
    xtrain = xtrain(ia,:);
    
    y_asfcv = cell2mat(ASFCV);
    y_asfcv_train = y_asfcv(TrainIndex,:);
    y_asfcv_train = y_asfcv_train(ia,:);
    y_asfcv_test = y_asfcv(TestIndex,:);
    
    
    selection_err = 0;
    counter = 0;
    n = size(xtest,1);
    
    active_size = min(opt.activeSetSize, size(xtrain, 1));
    rnd = randperm(size(xtrain,1));
    rnd = rnd(1:active_size);
    
    [opt.dmodel_asfcv, ~] = dacefit(xtrain(rnd,:), y_asfcv_train(rnd,:), @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model asf
    [y_pred,~] = predictor(xtest, opt.dmodel_asfcv);
    
    
    for dir = 1:1: size(opt.dirs,1) 
        
        %---------------------COMPUTE SEP----------------------------------
        for i=1:n-1
            for j=i+1:n
                d1 = constrained_domination(y_asfcv_test(i, dir), y_asfcv_test(j,dir), 0, 0);%1, 2, 3
                d2 = constrained_domination(y_pred(i, dir), y_pred(j, dir), 0, 0);
                if d1~=d2
                    selection_err = selection_err + 1;
                end
                counter = counter + 1;
            end
        end
        
        
               
    end
    %-----------------COMPUTE ASF+CV ERROR-----------------------------
    obj_error = sqrt(sum(((y_pred-y_asfcv_test)./y_asfcv_test).^2,2));
    sep = selection_err/counter;
    yerror  = mean(obj_error);
end
