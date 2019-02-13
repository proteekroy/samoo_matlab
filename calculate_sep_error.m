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

%This function calculates SEP for all frameworks

function [opt, sep] = calculate_sep_error(opt, TestIndex)
    
    xtest = opt.archive(TestIndex,:);
    test_obj = opt.archiveObj(TestIndex,:);
    cons_test = opt.archiveCons(TestIndex,:);
    SF_VAL_test = opt.archiveSF_VAL(TestIndex,:);

    cv_actual = zeros(size(xtest,1), 1);
    ASF = cell(1, size(opt.dirs,1) );
    ASFCV = cell(1, size(opt.dirs,1) );
    
    %return value
    sep = cell(1, 10);  
    count = cell(1, 10);
    Z_rank = cell(1, 10); 
    Z_rank_pred = cell(1, 10);
    
    for i=1:size(opt.dirs,1) 
        ASF{i} = opt.archiveASFAll{i}(TestIndex,:);
        ASFCV{i} = opt.archiveASFCVAll{i}(TestIndex,:);
    end
    
    %---------------------PREDICT OBJECTIVES-------------------------------
    [y_pred_normalized,~] = predictor(xtest,opt.dmodel_normalized_obj);
    [y_pred,~] = predictor(xtest, opt.dmodel_obj);%framework M1-2 and M2-2
        
    %---------------------PREDICT CONSTRAINTS------------------------------
    
    if opt.C>0
        [cons_pred,~] = predictor(xtest, opt.dmodel_cons);%framework M1-1, M1-2, M3-1, M3-2
        [cv_pred_M1,~] = evaluateCV(cons_pred);
        [cv_actual, ~] = evaluateCV(cons_test);
        [cv_pred_M2,~] = predictor(xtest, opt.dmodel_cv);%framework M2-1, M2-2, M4-1, M4-2
    else
        cv_pred_M1 = zeros(size(xtest,1), 1);
        cv_pred_M2 = zeros(size(xtest,1), 1);
    end
    
    %---------------------PREDICT ASF, ASFCV-------------------------------
    
    asf_pred_M1_1 = zeros(size(xtest,1), size(opt.dirs, 1));
    cv_M5 = zeros(size(xtest,1), 1);
    cv_M6 = zeros(size(xtest,1), 1);
    
    asf_pred_M3_1 = predictor(xtest, opt.dmodel_asf);%prediction
    asf_pred_M3_2 = predictor(xtest, opt.dmodel_asf);%prediction
    asf_pred_M5 = predictor(xtest, opt.dmodel_asfcv);%prediction
    
    asf_actual_M1_1 = cell2mat(ASF);
    asf_actual_M3_2 = asf_actual_M1_1;%actual      
    asf_actual_M5 = cell2mat(ASFCV);
    
    for dir = 1:1: size(opt.dirs,1) 
        %framework M1-1 and M2-1
        asf_pred_M1_1(:, dir) = calculate_Asf(y_pred_normalized, opt.dirs(dir,:), []);%prediction
    end
    
    %framework M3-2
    asf_pred_M3_2 = min(asf_pred_M3_2, [], 2);
    asf_actual_M3_2 = min(asf_actual_M3_2, [], 2);
    %asf_actual_M3_2 = calculate_m32_rank(asf_actual_M3_2, cv_actual);
    %asf_pred_M3_2 = calculate_m32_rank(asf_pred_M3_2, cv_pred_M2); 
    
    %framework M6
    %sv_pred = opt.net(xtest(:,1:opt.V)')';%prediction
    net = opt.prev_net{opt.crossval_index};
    sv_pred = net(xtest(:,1:opt.V)')';
    sv_actual = SF_VAL_test;%actual
    
    %compute and plot selection function
    dir = 1;
    [Z_rank{1}, Z_rank_pred{1}, sep{1}, count{1}] = compute_sf_framework(asf_actual_M1_1(:,dir), cv_actual, asf_pred_M1_1(:,dir), cv_pred_M1);%M1-1
    [Z_rank{3}, Z_rank_pred{3}, sep{3}, count{3}] = compute_sf_framework(asf_actual_M1_1(:,dir), cv_actual, asf_pred_M1_1(:,dir), cv_pred_M2);%M2-1
    [Z_rank{5}, Z_rank_pred{5}, sep{5}, count{5}] = compute_sf_framework(asf_actual_M1_1(:,dir), cv_actual, asf_pred_M3_1(:,dir), cv_pred_M1);%M3-1
    [Z_rank{7}, Z_rank_pred{7}, sep{7}, count{7}] = compute_sf_framework(asf_actual_M1_1(:,dir), cv_actual, asf_pred_M3_1(:,dir), cv_pred_M2);%M4-1
    [Z_rank{9}, Z_rank_pred{9}, sep{9}, count{9}] = compute_sf_framework(asf_actual_M5(:,dir), cv_M5, asf_pred_M5(:,dir), cv_M5);%M5                 
    
    for dir=2:size(opt.dirs,1)
        
        %framework M1-1
        [Z_rank_temp, Z_rank_pred_temp, sep_temp, count_temp] = compute_sf_framework(asf_actual_M1_1(:,dir), cv_actual, asf_pred_M1_1(:,dir), cv_pred_M1);
        Z_rank{1} = Z_rank_temp + Z_rank{1};
        Z_rank_pred{1} = Z_rank_pred_temp + Z_rank_pred{1};
        sep{1} = sep{1} + sep_temp; 
        count{1} = count{1} + count_temp;
        
        %framework M2-1
        [Z_rank_temp, Z_rank_pred_temp, sep_temp, count_temp] = compute_sf_framework(asf_actual_M1_1(:,dir), cv_actual, asf_pred_M1_1(:,dir), cv_pred_M2);
        Z_rank{3} = Z_rank_temp + Z_rank{3};
        Z_rank_pred{3} = Z_rank_pred_temp + Z_rank_pred{3};
        sep{3} = sep{3} + sep_temp; 
        count{3} = count{3} + count_temp;
        
        %framework M3-1
        [Z_rank_temp, Z_rank_pred_temp, sep_temp, count_temp] = compute_sf_framework(asf_actual_M1_1(:,dir), cv_actual, asf_pred_M3_1(:,dir), cv_pred_M1);
        Z_rank{5} = Z_rank_temp + Z_rank{5};
        Z_rank_pred{5} = Z_rank_pred_temp + Z_rank_pred{5};
        sep{5} = sep{5} + sep_temp; 
        count{5} = count{5} + count_temp;
        
        %framework M4-1
        [Z_rank_temp, Z_rank_pred_temp, sep_temp, count_temp] = compute_sf_framework(asf_actual_M1_1(:,dir), cv_actual, asf_pred_M3_1(:,dir), cv_pred_M2);
        Z_rank{7} = Z_rank_temp + Z_rank{7};
        Z_rank_pred{7} = Z_rank_pred_temp + Z_rank_pred{7};
        sep{7} = sep{7} + sep_temp; 
        count{7} = count{7} + count_temp;
        
        %framework M5
        [Z_rank_temp, Z_rank_pred_temp, sep_temp, count_temp] = compute_sf_framework(asf_actual_M5(:,dir), cv_M5, asf_pred_M5(:,dir), cv_M5);
        Z_rank{9} = Z_rank_temp + Z_rank{9};
        Z_rank_pred{9} = Z_rank_pred_temp + Z_rank_pred{9};
        sep{9} = sep{9} + sep_temp; 
        count{9} = count{9} + count_temp;
        
    end
    
    %framework M1-2
    [Z_rank{2}, Z_rank_pred{2}, sep{2}, count{2}] = compute_sf_framework(test_obj, cv_actual, y_pred, cv_pred_M1);
    %framework M2-2
    [Z_rank{4}, Z_rank_pred{4}, sep{4}, count{4}] = compute_sf_framework(test_obj, cv_actual, y_pred, cv_pred_M2);
    %framework M3-2
    [Z_rank{6}, Z_rank_pred{6}, sep{6}, count{6}] = compute_sf_framework(asf_actual_M3_2, cv_actual, asf_pred_M3_2, cv_pred_M1);
    %framework M4-2
    [Z_rank{8}, Z_rank_pred{8}, sep{8}, count{8}] = compute_sf_framework(asf_actual_M3_2, cv_actual, asf_pred_M3_2, cv_pred_M2);
    %framework M6
    [Z_rank{10}, Z_rank_pred{10}, sep{10}, count{10}] = compute_sf_framework(sv_actual, cv_M6, sv_pred, cv_M6);
    
    %find sep and average rank    
    for i=1:10
        sep{i} = sep{i}/count{i};
        Z_rank{i} = Z_rank{i}/count{i};
        Z_rank_pred{i} = Z_rank_pred{i}/count{i};
    end
    
    %{
    %plot model space created by metamodeling frameworks
    framework_name = ['M1-1';'M1-2';'M2-1';'M2-2';'M3-1';'M3-2';'M4-1';'M4-2';'M5  ';'M6  '];
    for i=1:10
        plot_sf_framework(test_obj, Z_rank_pred, framework_name(i,:));
    end
    %}
    opt.z_rank = Z_rank;
    opt.z_rank_pred = Z_rank_pred;
    sep = cell2mat(sep);

    
end



%This function computes selection function values for frameworks
function [Z_rank, Z_rank_pred, selection_err, count] = compute_sf_framework(ytest, test_cv, y_pred, pred_cv)

    n = size(ytest,1);
    Z_rank = zeros(n, 1);
    Z_rank_pred = zeros(n, 1);

    %---------------------COMPUTE SELECTION FUNCTION-----------------------
    selection_err = 0;
    count = 0;
    
    for i=1:n-1
        for j=i+1:n
            d1 = constrained_domination(ytest(i,:), ytest(j,:), test_cv(i,:), test_cv(j,:));%1, 2, 3
            d2 = constrained_domination(y_pred(i,:), y_pred(j,:), pred_cv(i,:), pred_cv(j,:));

            if d1==1
                Z_rank(j) = Z_rank(j)+1;
            elseif d1==2
                Z_rank(i) = Z_rank(i)+1;
            end

            if d2==1
                Z_rank_pred(j) = Z_rank_pred(j)+1;
            elseif d2==2
                Z_rank_pred(i) = Z_rank_pred(i)+1;
            end
            if d1~=d2
                selection_err = selection_err + 1;
            end
            count = count+1;
        end
    end

end



