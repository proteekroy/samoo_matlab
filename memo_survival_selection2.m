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


function [opt] = memo_survival_selection2(opt)

    
    %--------------TRUST REGION--------------------------------------------
    n  = size(opt.nsga2.hifipop, 1);
    allpop = normalize(opt, opt.nsga2.hifipop, opt.bound(1,:), opt.bound(2,:));
    norm_pop = normalize(opt, opt.nsga2.totalpop, opt.bound(1,:), opt.bound(2,:)); 
    div = pdist2(norm_pop, allpop);%distance in normalized space 
    index = any(   div  <=  repmat(opt.TrustRadiusDeltaK(1:n), size(norm_pop,1),1),      2);
    I1 = find(index==1);
    I2 = find(index==0);

    %-------------FIND CLUSTER NUMBERS-------------------------------------
    [ASF, pop_cluster] = find_memo_cluster(opt, opt.nsga2.totalpopObj);
    index1 = find(opt.nsga2.totalpopCV<=0);
    index2 = find(opt.nsga2.totalpopCV>0);
    if ~isempty(index1)
        temp = max(ASF(index1));
    else
        temp = 0;
    end
        
    ASF(index2) = temp + opt.nsga2.totalpopCV(index2);
    
    
    %---------------Select for Next Generation-----------------------------
    if size(I1,1)<=opt.nsga2.N
        selectedPopIndex = I1;
        temp_N = opt.nsga2.N - size(I1,1);
        if temp_N>0
            %[temp_selectedPopIndex, ~] = apply_nsga2_selection(opt, opt.nsga2.totalpop(I2,:), opt.nsga2.totalpopObj(I2,:), opt.nsga2.totalpopCV(I2,:), temp_N);
            temp_selectedPopIndex = apply_memo_selection(opt, ASF(I2,:), pop_cluster(I2), temp_N);
            
            if ~isempty(selectedPopIndex)
                selectedPopIndex = vertcat(selectedPopIndex, I2(temp_selectedPopIndex));
            else
                selectedPopIndex = I2(temp_selectedPopIndex);
            end
            
        end
    else %size(I1,1)>opt.nsga2.N
        %[temp_selectedPopIndex, ~] = apply_nsga2_selection(opt, opt.nsga2.totalpop(I1,:), opt.nsga2.totalpopObj(I1,:), opt.nsga2.totalpopCV(I1,:), opt.nsga2.N);
        temp_selectedPopIndex = apply_memo_selection(opt, ASF(I1,:), pop_cluster(I1), opt.nsga2.N);
        selectedPopIndex = I1(temp_selectedPopIndex);
    end
    
    [~,opt.nsga2.CD] = apply_nsga2_selection(opt, opt.nsga2.totalpop, opt.nsga2.totalpopObj, opt.nsga2.totalpopCV, opt.nsga2.N);
    %[~,opt.nsga2.CD] = apply_nsga3_selection(opt, opt.nsga2.totalpop, opt.nsga2.totalpopObj, opt.nsga2.totalpopCV, opt.nsga2.N);
     
    
    opt.nsga2.pop =  opt.nsga2.totalpop(selectedPopIndex,:);
    opt.nsga2.popObj = opt.nsga2.totalpopObj(selectedPopIndex,:);
    opt.nsga2.popCV = opt.nsga2.totalpopCV(selectedPopIndex,:);
    opt.nsga2.popCons = opt.nsga2.totalpopCons(selectedPopIndex,:);
    opt.nsga2.CD = opt.nsga2.CD(selectedPopIndex,:);
    

end


