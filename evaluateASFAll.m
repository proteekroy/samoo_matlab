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
% Contact: royprote@egr.msu.edu, proteek_buet@yahoo.com

%FIND ASF ACCORDING TO EACH DIRECTION

function [opt] = evaluateASFAll(opt)
    
    opt.archiveASFAll = cell(1, opt.numdir);
    opt.archiveASFCVAll = cell(1, opt.numdir);
    for i=1:size(opt.dirs,1)
        w = opt.dirs(i,:);
        asf = calculate_Asf(opt.normalizedObj, w);
        opt.archiveASFAll{i} = asf; 
    end
    
    opt.archiveASFCVAll = opt.archiveASFAll;
    feasbile_index = opt.archiveCV<=0;
    
    for i=1:size(opt.dirs,1)
        if opt.C>0
            feasibleASF =  opt.archiveASFAll{i}(feasbile_index,:);
            
            if ~isempty(feasibleASF)
                fmax = max(feasibleASF);%maximum feasible ASF
            else
                fmax = 0;
            end
            infeasible_index = find(opt.archiveCV>0);%infeasible indices
            if ~isempty(infeasible_index)
                opt.archiveASFCVAll{i}(infeasible_index) = fmax + opt.archiveCV(infeasible_index);
            end
        end
    end
    
end

function asf = calculate_Asf(obj, dir)   
    asf = max( (obj-repmat(dir, size(obj,1),1)),[],2);
end

