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


function runKfoldParallel()%(k, opt)

    methods = cell(1,5);
    methods{1} = 'methodology1_2_crossval_main';
    methods{2} = 'methodology2_2_crossval_main';
    methods{3} = 'methodology3_crossval_main';
    methods{4} = 'methodology4_crossval_main';
    methods{5} = 'methodology5_crossval_main';
    

    %----------------CROSS_VALIDATION_SET----------------------------------
    k = 10;
    opt.archive = rand(100, 10);
    Indices = crossvalind('Kfold', size(opt.archive,1), k);
    c = cvpartition(size(opt.archive,1), 'KFold', k);
    disp(c);
    TrainIndex = cell(1,10);
    TestIndex = cell(1,10);
    
    for i=1:10
        TrainIndex{i} = find(Indices~=i);
        TestIndex{i} = find(Indices==i);
    end
    
    %----------------DO_CROSS_VALIDATION-----------------------------------
    opt.r = 1;
    opt.objfunction = strtrim('zdt1');%remove whitespaces
    opt.func_no = 1;
    opt.methodology = 2; 
    [opt] = basic_parameters(opt);
    opt.dirs = initweight(2, 21)';
    opt.numdir = 21;
    opt.M = 2;
    
    opt.archiveASF = cell(1,21);
    for i=1:21
        opt.archiveASF{i} = rand(100,1);
    end
    opt.archiveCons = rand(100,2);
    index =  randperm(100);
    oneMatrix = -1*ones(1,60);
    opt.archiveCons(index(1:60),:) = repmat(oneMatrix',1,2);
    opt.archive = rand(100, 3);
    opt.archiveObj = rand(100,2);
    opt.C = 1;
    
    %--------------RUN-CROSSVALIDATION IN PARALLEL-------------------------
    num_of_methodology = 5;
    F = cell(1, 5);
    %f = cell(1, num_of_methodology);
    R = cell(1, num_of_methodology);
    p = gcp();
    for i=1:num_of_methodology
        fh = str2func(methods{i});
        F{i} = parfeval(p, fh, 1, opt, TrainIndex, TestIndex, k);
    end
    
    for i=1:num_of_methodology
        R{i} = fetchOutputs(F{i});
        disp(R{i}{1});
    end
    
    disp('done');
    
    S = [];
    
    for i=1:num_of_methodology
        S = vertcat(S, R{i}(1));
    end
    S = cell2mat(S);
    mean_S = mean(S, 2);
    [~, I] = min(mean_S);
    minMeanRow = S(I,:);

    % Loop over all other rows and add those that are not significantly
    % worse.
    indices = I;
    for i = 1 : num_of_methodology
        if i ~= I
            p = ranksum(minMeanRow, S(i,:));
            if p >= 0.05
                indices = [indices, i];
            end
        end
    end
    % Pick one row (methodology) randomly
    selectedMethodology = indices(randi(size(indices, 2)));
    opt.methodology = selectedMethodology;
    %--------------WAIT UNTIL ALL THE PROCESS IS FINISHED------------------
%     for i=1:num_of_methodology
%         [completedIdx,value] = fetchNext(f);
%         F{completedIdx} = value;
%         fprintf('Got result with index: %d.\n', completedIdx);
%     end
    
    %----------------DISPLAY_ERROR-----------------------------------------
%     for i=1:num_of_methodology
%         error{i} = mean(F{i});
%         disp(error{i});
%     end
    
    
    
end