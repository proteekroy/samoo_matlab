% Copyright [2017] [Proteek Chandan Roy]
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
% Proteek Chandan Roy, 2017
% Contact: royprote@egr.msu.edu, proteek_buet@yahoo.com



function [opt] = survival_selection(opt)

    switch(opt.nsga2.survivalselectionOption)
        
        case 1 %------------------NSGA-II----------------------------------
            
            %opt = nsga2_selection(opt);
            opt = trusted_nsga2_selection(opt);
            
        case 2%-------------------NSGA-III---------------------------------
            
            opt = nsga3_selection(opt);
            
        case 3%-----------------PageRank-NDS-------------------------------
            
            opt = pagerank_selection(opt);
        
        case 4%------------------NSGA-II with TRUST REGION-----------------
            
            opt = trust_region_nsga2_selection(opt);
            
            
        otherwise
            
            
            disp('Selection is not defined');
            
    end

end



