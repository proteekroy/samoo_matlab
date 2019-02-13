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


function opt = trusted_nsga3_selection(opt)
    
    if opt.trust_region_option_active==1
        
         switch(opt.trust_region_update_option)
            
            case 1 % continuous decreasing
                allpop = normalize(opt, opt.archive, opt.bound(1,:), opt.bound(2,:));
                norm_pop = normalize(opt, opt.nsga2.totalpop, opt.bound(1,:), opt.bound(2,:)); 
                div = pdist2(norm_pop, allpop);%distance in normalized space 
                div = min(div,[],2);
                I1 = find(div<=opt.TrustDistVar);
                I2 = find(div>opt.TrustDistVar);
            case 2 % non-linear/diffusion trust region
                leader = normalize(opt, opt.archive(opt.ParetoIndex,:), opt.bound(1,:), opt.bound(2,:));
                if isempty(leader)
                    leader = normalize(opt, opt.archive, opt.bound(1,:), opt.bound(2,:));
                end
                norm_pop = normalize(opt, opt.nsga2.totalpop, opt.bound(1,:), opt.bound(2,:));
                div = trust_distance(opt, norm_pop, leader);
                div = min(div,[],2);
                div = (div - min(div))/(max(div)-min(div));
                I1 = find(div<=opt.delta);
                I2 = find(div>opt.delta);                
            case 3 % adaptive
                n  = size(opt.archive, 1);
                allpop = normalize(opt, opt.archive, opt.bound(1,:), opt.bound(2,:));
                norm_pop = normalize(opt, opt.nsga2.totalpop, opt.bound(1,:), opt.bound(2,:)); 
                div = pdist2(norm_pop, allpop);%distance in normalized space 
                index = any(   div  <=  repmat(opt.TrustRadiusDeltaK(1:n), size(norm_pop,1),1),      2);
                I1 = find(index==1);
                I2 = find(index==0);
                
            otherwise
                input('trust region option invalid');
        end
        
    else %no trust region option
        I1 = (1:1:size(opt.nsga2.totalpop, 1))';%every solution is valid
        I2 = [];
    end  
    
    %----------------------------------------------------------------------
    %---------------Select for Next Generation-----------------------------
    %----------------------------------------------------------------------
    
    if size(I1,1)>opt.nsga2.N % we have more solutions inside trust region 
        
        I1 = vertcat(I1, I2);
        [temp_selectedPopIndex, opt] = apply_nsga3_selection(opt, opt.nsga2.totalpop(I1,:), opt.nsga2.totalpopObj(I1,:), opt.nsga2.totalpopCV(I1,:), opt.nsga2.N);
        selectedPopIndex = I1(temp_selectedPopIndex);
        
    else % we have more solutions inside trust region, so select all from inside and some from outside based on distance
         selectedPopIndex = I1;
         temp_N = opt.nsga2.N - size(I1,1);
         if temp_N>0
            div = min(div,[],2);
            [~, I3] = sort(div);
            I3 = setdiff(I3, I1);
            temp_selectedPopIndex = I3(1:temp_N);
            selectedPopIndex = vertcat(selectedPopIndex, temp_selectedPopIndex);
         end
    end
    
    %[~,opt.nsga2.CD] = apply_nsga2_selection(opt, opt.nsga2.totalpop, opt.nsga2.totalpopObj, opt.nsga2.totalpopCV, opt.nsga2.N);  
    
    opt.nsga2.pop =  opt.nsga2.totalpop(selectedPopIndex,:);
    opt.nsga2.popObj = opt.nsga2.totalpopObj(selectedPopIndex,:);
    opt.nsga2.popCV = opt.nsga2.totalpopCV(selectedPopIndex,:);
    opt.nsga2.popCons = opt.nsga2.totalpopCons(selectedPopIndex,:);
    %opt.nsga2.CD = opt.nsga2.CD(selectedPopIndex,:);
    
    
       
end


%---------------------------END OF FILE------------------------------------