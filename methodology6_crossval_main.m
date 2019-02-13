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

%This file performs k-fold cross valiaditon for Framework 6

function  [returnval] = methodology6_crossval_main(opt, TrainIndex, TestIndex, k)

    returnval = cell(1,2);
    sep = zeros(1,10);
    yerror = zeros(1,10);
    m6_net = cell(1, 10);
    
    for i=1:k
        [sep(i), yerror(i), m6_net{i}] = methodology6_crossval(opt, opt.archive(TrainIndex{i},:), opt.archive(TestIndex{i},:), ...
            opt.archiveSF_VAL, TrainIndex{i}, TestIndex{i});
    end
    
    returnval{1} = sep; 
    returnval{2} = yerror;
    returnval{3} = m6_net;
end


function [sep, yerror, net] =  methodology6_crossval(opt, xtrain, xtest, SF_VAL, TrainIndex, TestIndex)

    [~,ia,~] = unique(xtrain,'rows');
    xtrain = xtrain(ia,:);
    TrainIndex = TrainIndex(ia,:);
    
    [~,ia,~] = unique(xtest,'rows');
    xtest = xtest(ia,:); 
    TestIndex = TestIndex(ia,:);
       
    selection_err = 0;
    counter = 0;
    n = size(xtest,1);
    
    ytrain = SF_VAL(TrainIndex);
    ytest = SF_VAL(TestIndex);
        
    %----------------MODEL SELECTION FUNCTION------------------------------    
    net = train_neural_network(xtrain, ytrain); 
    
    %-----------------PREDICT SELECTION FUNCTION---------------------------
    y_pred = net(xtest(:,1:opt.V)'); %prediction
    y_pred = y_pred';  
    
    %[opt.dmodel_asf, ~] = dacefit(xtrain, ytrain, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model asf
    %[y_pred,~] = predictor(xtest, opt.dmodel_asf);
    
    %---------------------COMPUTE SEP--------------------------------------
    for i=1:n-1
        for j=i+1:n
            d1 = constrained_domination(ytest(i), ytest(j), 0, 0);%1, 2, 3
            d2 = constrained_domination(y_pred(i), y_pred(j), 0, 0);
            if d1~=d2
                selection_err = selection_err + 1;
            end
            counter = counter + 1;
        end
    end
        
    %-----------------COMPUTE SELECTION FUNCTION ERROR---------------------
    obj_error = sqrt(sum(((y_pred-ytest)./ytest).^2,2));
               
    sep = selection_err/counter;
    yerror  = mean(mean(obj_error, 2));
end
