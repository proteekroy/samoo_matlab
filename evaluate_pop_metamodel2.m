% Copyright [2016] [Proteek Chandan Roy]
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
% Proteek Chandan Roy, 2016
% Contact: royprote@msu.edu, proteek_buet@yahoo.com


%--------------------EVALUATE BY METAMODEL---------------------------------

function [popObj, popCons] = evaluate_pop_metamodel2(opt, pop)


    popCons = zeros(size(pop,1), opt.C);
    popObj = zeros(size(pop,1), opt.M);
    %%{
    popCons(:,1) = quadratic_predictor(opt.regC{1}, pop(:,1:opt.V));
    popCons(:,2) = quadratic_predictor(opt.regC{2}, pop(:,1:opt.V));
    popCons(:,3) = quadratic_predictor(opt.regC{3}, pop(:,1:opt.V));
    popCons(:,4) = quadratic_predictor(opt.regC{4}, pop(:,1:opt.V));
    popCons(:,5) = quadratic_predictor(opt.regC{5}, pop(:,1:opt.V));
    popCons(:,6) = quadratic_predictor(opt.regC{6}, pop(:,1:opt.V));
                    
    popObj(:,1) = quadratic_predictor(opt.regO{1}, pop(:,1:opt.V));
    popObj(:,2) = quadratic_predictor(opt.regO{2}, pop(:,1:opt.V));
    %}
    
    %{
    popObj(:,1) = quadratic_predictor(opt.regO{1}, pop(:,1:opt.V));
    popObj(:,2) = quadratic_predictor(opt.regO{2}, pop(:,1:opt.V));
    
    [popCons(:,1),~,~,~] = predictor(pop(:,1:opt.V), opt.dmodel_cons);
    popCons(:,2) = quadratic_predictor(opt.regC{1}, pop(:,1:opt.V));
    %}
    
end



