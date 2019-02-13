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


function [method_number] = runCrossValidation(opt)

    num_of_methodology = opt.num_of_frameworks;
    methods = cell(1, num_of_methodology);
    k = 1;
    if ismember(11,opt.framework_list)
        methods{k} = 'methodology1_1_crossval_main'; 
        k = k+1;
    end
    if ismember(12,opt.framework_list)
        methods{k} = 'methodology1_2_crossval_main'; 
        k = k+1;
    end
    if ismember(21,opt.framework_list)
        methods{k} = 'methodology2_1_crossval_main'; 
        k = k+1;
    end
    if ismember(22,opt.framework_list)
        methods{k} = 'methodology2_2_crossval_main'; 
        k = k+1;
    end
    if ismember(31,opt.framework_list)
        methods{k} = 'methodology3_1_crossval_main'; 
        k = k+1;
    end
    if ismember(32,opt.framework_list)
        methods{k} = 'methodology3_2_crossval_main'; 
        k = k+1;
    end
    if ismember(41,opt.framework_list)
        methods{k} = 'methodology4_1_crossval_main'; 
        k = k+1;
    end
    if ismember(42,opt.framework_list)
        methods{k} = 'methodology4_2_crossval_main'; 
        k = k+1;
    end
    if ismember(5,opt.framework_list)
        methods{k} = 'methodology5_crossval_main'; 
        k = k+1;
    end
    if ismember(6,opt.framework_list)
        methods{k} = 'methodology6_crossval_main'; 
    end
      
    %{
    if opt.switching_option==3
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
    elseif opt.switching_option==1
        methods{1} = 'methodology1_1_crossval_main';
        methods{2} = 'methodology1_2_crossval_main';
        methods{3} = 'methodology2_1_crossval_main';
        methods{4} = 'methodology3_1_crossval_main';
        methods{5} = 'methodology4_1_crossval_main';
        methods{6} = 'methodology5_crossval_main';
    elseif opt.switching_option==2
        methods{1} = 'methodology1_2_crossval_main';
        methods{2} = 'methodology2_2_crossval_main';
        methods{3} = 'methodology3_2_crossval_main';
        methods{4} = 'methodology4_2_crossval_main';
        methods{5} = 'methodology6_crossval_main';
    end
    %}
    
    fh = cell(1,num_of_methodology);
    for i=1:num_of_methodology
        fh{i} = str2func(methods{i});
    end
    %----------------CROSS_VALIDATION_SET----------------------------------
    K = 10;
    Indices = crossvalind('Kfold', size(opt.archive,1), K);
    TrainIndex = cell(1,10);
    TestIndex = cell(1,10);
    
    for i=1:10
        TrainIndex{i} = find(Indices~=i);
        TestIndex{i} = find(Indices==i);
    end
    
    %--------------RUN-CROSSVALIDATION IN PARALLEL-------------------------
    
    F = cell(1, num_of_methodology);
    R = cell(1, num_of_methodology);
    
    %%{
    %parpool(10);
    %delete(p);
    %delete(gcp('nocreate'))
    
    p = gcp();%('nocreate');
    %t = cputime;
    method_counter = 0;
    for j=method_counter+1:num_of_methodology
        working_cores = min(num_of_methodology-method_counter, p.NumWorkers);
        for i=1:working_cores
            F{i} = parfeval(p, fh{method_counter+i}, 1, opt, TrainIndex, TestIndex, K);
        end
        
        for i=1:working_cores
            R{method_counter+i} = fetchOutputs(F{i});
            %disp(R{i}{1});
            if method_counter+i==10
                prev_net = R{method_counter+i}{3};
            end
        end
        method_counter = method_counter+working_cores;
    end

    %elapse_time_prev = cputime-t;
    %}
    
    
    %======================RUN CROSSVAL PARALLEL===========================    
    %{
    t = cputime;
    crossval_counter = 0;
    fh_model = @collect_models;
    R_opt = cell(1, K);
    G = cell(1, K);
    SEP2 = zeros(num_of_methodology, K);
    
    for j=crossval_counter+1:K
        working_cores = min(K-crossval_counter, p.NumWorkers);
        for i=1:working_cores
            G{i} = parfeval(p, fh_model, 1, opt, TrainIndex{crossval_counter+i});
        end
        
        for i=1:working_cores
            R_opt{crossval_counter+i} = fetchOutputs(G{i});
            R_opt{crossval_counter+i}.prev_net = prev_net;
            R_opt{crossval_counter+i}.crossval_index = crossval_counter+i;
            [~, temp] = calculate_sep_error(R_opt{crossval_counter+i}, TestIndex{crossval_counter+i});
            SEP2(:,crossval_counter+i) = temp';
        end
        crossval_counter = crossval_counter+working_cores;
    end

    elapse_time_current = cputime-t;
    disp(elapse_time_prev);
    disp(elapse_time_current);
    
    %}
    
    %=======================RUN ONE BY ONE=================================
    %{
    for i=1:num_of_methodology%for each metamodel
        fh = str2func(methods{i});
        R{i} = cell(1, 2);
        R{i} = feval(fh, opt, TrainIndex, TestIndex, K);
        %disp(R{i}{1})
        %disp(R{i}{2})
        if i==10
            opt.prev_net = R{i}{3};
        end
    end
    
    %}
    %create all the necessary models----can be parallelized
    %{
    SEP2 = zeros(num_of_methodology, K);
    
    for i=1:K%for each cross-val partition
        fh = @collect_models;
        R_opt = feval(fh, opt, TrainIndex{i});
        R_opt.prev_net = opt.prev_net;
        R_opt.crossval_index = i;
        [~, temp] = calculate_sep_error(R_opt, TestIndex{i});
        SEP2(:,i) = temp';
    end
    %}     
    
    
    %=====================CHOOSE A MODEL===================================
    SEP = []; 
    MSE = [];
    for i=1:num_of_methodology
        SEP = vertcat(SEP, R{i}{1});
        MSE = vertcat(MSE, R{i}{2});
    end
    
    mean_SEP = mean(SEP, 2);
    mean_MSE = mean(MSE, 2);
    
    %disp(norm(SEP-SEP2));
       
    %dlmwrite('error_data_sep.txt', ['Epoch: ' num2str(opt.iter)], 'delimiter',' ','precision','%.10f','-append');
    %dlmwrite('error_data_mse.txt', ['Epoch: ' num2str(opt.iter)], 'delimiter',' ','precision','%.10f','-append');
    dlmwrite(opt.sep_filename, SEP, 'delimiter',' ','precision','%.10f','-append');
    dlmwrite(opt.mse_filename, MSE, 'delimiter',' ','precision','%.10f','-append');
       
    
    [~, I_SEP] = min(mean_SEP);
    minMeanRowSEP = SEP(I_SEP,:);
    [~, I_MSE] = min(mean_MSE);
    minMeanRowMSE = MSE(I_MSE,:);

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
    method_number = opt.framework_list(method_no);
    
    a = 1:num_of_methodology;
    equivalent_methods = int8(ismember(a, indices));
    equivalent_methods(method_no) = 2;
    dlmwrite(opt.selected_method_filename, equivalent_methods, 'delimiter',' ','-append');
        
    
    
end