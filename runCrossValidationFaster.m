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


function [method_number] = runCrossValidationFaster(opt)

    %------------Switching Framework-III-----------------------------------
    num_of_methodology = 10;
    methods = cell(1, num_of_methodology);
    methods{1} = 'methodology1_1_crossval_main';
    methods{2} = 'methodology1_2_crossval_main';
    methods{3} = 'methodology2_1_crossval_main';
    methods{4} = 'methodology2_2_crossval_main';
    methods{5} = 'methodology3_1_crossval_main';
    methods{6} = 'methodology3_2_crossval_main';
    methods{7} = 'methodology4_1_crossval_main';
    methods{8} = 'methodology4_2_crossval_main';
    methods{9} = 'methodology5_crossval_main';
    methods{10} = 'methodology6_crossval_main';

    %------------Switching Framework-I-------------------------------------
%     num_of_methodology = 6;
%     methods = cell(1, num_of_methodology);
%     methods{1} = 'methodology1_1_crossval_main';
%     methods{2} = 'methodology1_2_crossval_main';
%     methods{3} = 'methodology2_1_crossval_main';
%     methods{4} = 'methodology3_1_crossval_main';
%     methods{5} = 'methodology4_1_crossval_main';
%     methods{6} = 'methodology5_crossval_main';

    %------------Switching Framework-II------------------------------------
%     num_of_methodology = 5;
%     methods = cell(1, num_of_methodology);
%     methods{1} = 'methodology1_2_crossval_main';
%     methods{2} = 'methodology2_2_crossval_main';
%     methods{3} = 'methodology3_2_crossval_main';
%     methods{4} = 'methodology4_2_crossval_main';
%     methods{5} = 'methodology6_crossval_main';
    
    
    %----------------CROSS_VALIDATION_SET----------------------------------
    K = 10;
    Indices = crossvalind('Kfold', size(opt.archive,1), K);
    TrainIndex = cell(1,10);
    TestIndex = cell(1,10);
    
    for i=1:10
        TrainIndex{i} = find(Indices~=i);
        TestIndex{i} = find(Indices==i);
    end
    
    %create all the necessary models----can be parallelized
    temp_opt = cell(1, K);
    R_opt = cell(1, K);
    SEP = zeros(num_of_methodology, K);
    
    %--------------RUN-CROSSVALIDATION IN PARALLEL-------------------------
    %{
    p = gcp();%('nocreate');
    
    for i=1:K
        fh = @collect_models;
        temp_opt{i} = parfeval(p, fh, 1, opt, TrainIndex{i});
    end
    for i=1:10
        R_opt{i} = fetchOutputs(temp_opt{i});
        [~, temp] = calculate_sep_error(R_opt{i}, TrainIndex{i}, TestIndex{i});
        SEP(i,:) = temp';
    end  
    %}
    %%{
    for i=1:K
        fh = @collect_models;
        R_opt = feval(fh, opt, TrainIndex{i});
        [~, temp] = calculate_sep_error(R_opt, TestIndex{i});
        SEP(:,i) = temp';
    end
    %}
    
    mean_SEP = mean(SEP, 2);
    %mean_MSE = mean(MSE, 2);
       
    %dlmwrite('error_data_sep.txt', ['Epoch: ' num2str(opt.iter)], 'delimiter',' ','precision','%.10f','-append');
    %dlmwrite('error_data_mse.txt', ['Epoch: ' num2str(opt.iter)], 'delimiter',' ','precision','%.10f','-append');
    %dlmwrite(opt.sep_filename, SEP, 'delimiter',' ','precision','%.10f','-append');
    %dlmwrite(opt.mse_filename, MSE, 'delimiter',' ','precision','%.10f','-append');
    
    [~, I_SEP] = min(mean_SEP);
    minMeanRowSEP = SEP(I_SEP,:);
    %[~, I_MSE] = min(mean_MSE);
    %minMeanRowMSE = MSE(I_MSE,:);

    % Loop over all other rows and add those that are not significantly
    % worse.
    indices = I_SEP;
    for i = 1 : num_of_methodology
        if i ~= I_SEP
            p = ranksum(minMeanRowSEP, SEP(i,:));
            if p >= 0.05 || isnan(p)
                indices = [indices, i];
            end
        end
    end
    % Pick one row (methodology) randomly
    method_no = indices(randi(size(indices, 2)));
    
    switch(method_no)
        case 1
            method_number = 11;
        case 2
            method_number = 12;
        case 3
            method_number = 21;
        case 4
            method_number = 22;
        case 5
            method_number = 31;
        case 6
            method_number = 32;
        case 7
            method_number = 41;
        case 8
            method_number = 42;
        case 9
            method_number = 5;
        case 10
            method_number = 6;
        otherwise
            method_number = 1000;
    end
    
    a = 1:num_of_methodology;
    equivalent_methods = int8(ismember(a, indices));
    equivalent_methods(method_no) = 2;
    dlmwrite(opt.selected_method_filename, equivalent_methods, 'delimiter',' ','-append');
        
    
    
end