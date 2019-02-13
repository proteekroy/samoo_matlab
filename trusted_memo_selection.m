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


function opt = trusted_memo_selection(opt)
    
    %---------DIVIDE POPULATION INSIDE AND OUTSIDE TRUST REGION------------
    if opt.trust_region_option_nsga2==1
        allpop = normalize(opt, opt.nsga2.hifipop, opt.bound(1,:), opt.bound(2,:));
        norm_pop = normalize(opt, opt.nsga2.totalpop, opt.bound(1,:), opt.bound(2,:)); 
        div = pdist2(norm_pop, allpop);%distance in normalized space
        div = min(div,[],2);
        I1 = find(div<=opt.TrustDistVar);
        I2 = find(div>opt.TrustDistVar);
    elseif opt.trust_region_option_nsga2==2
        leader = normalize(opt, opt.archive(opt.ParetoIndex,:), opt.bound(1,:), opt.bound(2,:));
        if isempty(leader)
            leader = normalize(opt, opt.nsga2.hifipop, opt.bound(1,:), opt.bound(2,:));
        end
        norm_pop = normalize(opt, opt.nsga2.totalpop, opt.bound(1,:), opt.bound(2,:));
        div = trust_distance(opt, norm_pop, leader);
        div = min(div,[],2);
        I1 = find(div<=opt.delta);
        I2 = find(div>opt.delta);
    elseif opt.trust_region_option_nsga2==3
        I1 = (1:1:size(opt.nsga2.totalpop, 1))';
        I2 = [];
    end  

    %--------------FIND CLUSTER MEMBERS------------------------------------
    [asf, pop_cluster] = find_memo_cluster(opt, opt.nsga2.totalpopObj);
     
    %---------------SELECTION FOR NEXT GENERATION--------------------------
    if size(I1,1)<=opt.nsga2.N
        selectedPopIndex = I1;
        temp_N = opt.nsga2.N - size(I1,1);
        
        if temp_N>0
            temp_selectedPopIndex = apply_memo_selection(opt, asf(I2,:),pop_cluster(I2), temp_N);
            
            if ~isempty(selectedPopIndex)
                selectedPopIndex = vertcat(selectedPopIndex, I2(temp_selectedPopIndex));
            else
                selectedPopIndex = I2(temp_selectedPopIndex);
            end
            
        end
    else %size(I1,1)>opt.nsga2.N
        temp_selectedPopIndex = apply_memo_selection(opt, asf(I1,:),pop_cluster(I1), opt.nsga2.N);
        selectedPopIndex = I1(temp_selectedPopIndex);
    end
       
    opt.nsga2.pop =  opt.nsga2.totalpop(selectedPopIndex,:);
    opt.nsga2.popObj = opt.nsga2.totalpopObj(selectedPopIndex,:);
    opt.nsga2.popCV = opt.nsga2.totalpopCV(selectedPopIndex,:);
    opt.nsga2.popCons = opt.nsga2.totalpopCons(selectedPopIndex,:);  
    opt.nsga2.pop_cluster = pop_cluster(selectedPopIndex);
    opt.nsga2.asf = asf(selectedPopIndex);
       
end

