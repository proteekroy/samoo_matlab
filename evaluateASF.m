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


%CLUSTER IS DETERMINED BY M-EMO method

function [asf, cluster] = evaluateASF(opt)

    
    asf = zeros(size(opt.normalizedObj,1),1);
    
    switch(opt.methodology)
        
        %case 6
            
            %cluster = find_euclidean_cluster(opt, opt.archiveObj);
        %    [asf, cluster] = find_memo_cluster(opt, opt.archiveObj);
        
        case {11, 12, 21, 22, 31, 32, 41, 42} %find ASF w.r.t one direction, targeting one optimal point
            
            cluster  = repmat(opt.curcluster,size(opt.normalizedObj,1),1);
            asf = calculate_Asf(opt.normalizedObj, opt.dirs(opt.curcluster,:), opt.utopian);
            %asf = calculate_Asf(opt.archiveObj, opt.dirs(opt.curcluster,:), opt.utopian);
            %asf = calculate_Asf(opt.archiveObj, ref_point{opt.curcluster}, opt.utopian);
        
        case {5}
            
            cluster  = repmat(opt.curcluster,size(opt.normalizedObj,1),1);
            asf = calculate_Asf(opt.normalizedObj, opt.dirs(opt.curcluster,:), opt.utopian);
            
            %-------------COMBINE CV WITH ASF--------------------------
            if opt.C>0 
                index = opt.archiveCV<=0;
                feasibleASF = asf(index,:);
                fmax = max(feasibleASF);
                if isempty(fmax)
                    fmax = 0;%max(asf);
                end
                %asf = asf./fmax;
                infeasible_index = find(opt.archiveCV>0);
                
                if ~isempty(infeasible_index)
                    temp_asf = fmax + opt.archiveCV(infeasible_index);
                    asf(infeasible_index) = temp_asf;
                end
            end
            
            opt.bestASFCV = min(asf);
            
    
        case {6} %MEMO approach, similar to case 2, only code is different, without calculating distance to asf directions
             
            asfcv = cell2mat(opt.archiveASFCVAll);
            [asf, cluster] = min(asfcv,[],2);
            
            %{
            %-------------COMBINE CV WITH ASF------------------------------
            if opt.C>0 
                index = opt.archiveCV<=0;
                feasibleASF = asf(index,:);
                fmax = max(feasibleASF);
                if isempty(fmax)
                    fmax = 0;%max(asf);
                end
                %asf = asf./fmax;
                infeasible_index = find(opt.archiveCV>0);
                
                if ~isempty(infeasible_index)
                    temp_asf = fmax + opt.archiveCV(infeasible_index);
                    asf(infeasible_index) = temp_asf;
                end
            end
            %}
            opt.bestASFCV = min(asf);
            
            
        case {-1} % It uses the euclidean distance, which may not be equal to smallest ASF distance
            
            cluster = zeros(size(opt.normalizedObj,1),1);
            
            for i=1:size(opt.normalizedObj,1)
                min_val = intmax;
                for j=1:size(opt.dirs,1)
                    w = opt.dirs(j,:);
                    obj = opt.normalizedObj(i,:);
                    d = norm(obj-((w*obj')*w)/(norm(w)*norm(w)));
                    
                    if(d<min_val)%find the direction for which it is closest 
                       min_val = d;
                       cluster(i) = j;
                    end
                end

                %find asf value for that closest direction
                asf(i) = calculate_Asf(opt.normalizedObj(i,:), opt.dirs(cluster(i),:), opt.utopian);
            end
            
        otherwise
            input('Inside evaluateASF, ASF calculation function not defined');
    end
        
  
    
    
end

function d = calculateNormDist(opt, obj, t)
    if(opt.M==2)
        d = point_to_line([obj 0],[0,0,0],[t 0]);
    elseif opt.M==3
       d = point_to_line(obj,[0,0,0],t);   
    end
end

function d = point_to_line(pt, v1, v2)

      a = v1 - v2;
      b = pt - v2;
      d = norm(cross(a,b)) / norm(a);
      
end

% function asf = calculate_Asf(obj, dir, ~)%utopian)
% 
%     %w = (dir-utopian)./norm(dir-utopian);%unit direction
%     %w = dir;
%     %asf  = max((obj-utopian)./w);
%     
%     %dir = dir./norm(dir);
%     %asf = max(obj.*dir);
%     %asf  = max(obj./dir); % max(w.*(obj-utopian));  
%     
%     %No Normalization is done
%     %asf = max(obj - repmat(dir-1/size(obj,2),size(obj,1),1),[],2);%+0.01*sum(obj-(dir-0.5));%+0.0001*sum(obj.*dir);
%     %asf = max(obj - repmat(dir + (1/size(obj,2)),size(obj,1),1),[],2);%+0.01*sum(obj-(dir-0.5));%+0.0001*sum(obj.*dir);
%     
%     asf = max( (obj-repmat(dir, size(obj,1),1)),[],2);
% end



%     %--------------FIND CLUSTER--------------------------------------------
%     S = zeros(size(opt.normalizedObj,1),1);
%     S(1:end)=intmax;
% 
%     cluster = ones(size(opt.normalizedObj,1),1);
%     minASFcluster = zeros(opt.numdir,1);
%     minASFcluster(1:end) = intmax;
%     ref_point = cell(1, opt.numdir);
%     
%     for i=1:size(opt.dirs,1)
%         w = opt.dirs(i,:);
%         %w = w./norm(w);
%         %w_org = w;
%         w_org = zeros(1,opt.M);
%         for j=1:opt.M
%             w_org(j) =  opt.min_val(j)+ w(j)*(opt.max_val(j)-opt.min_val(j));
%         end
%         %Z  = max((opt.archiveObj-repmat(opt.min_val,size(opt.archiveObj,1),1))./repmat(w_org,size(opt.archiveObj,1),1),[],2);
%         %Z  = max(opt.normalizedObj./repmat(w_org,size(opt.normalizedObj,1),1),[],2);
%         small_z = (opt.max_val-opt.min_val)/2;
%         ref_point{i} = w_org - small_z;
%         Z = max(      (opt.archiveObj-repmat(ref_point{i}, size(opt.archiveObj,1),1)),    [],2);
%         
%         for j = 1:size(opt.normalizedObj,1)
%             if Z(j)<S(j)
%                 cluster(j) = i;
%                 S(j) = Z(j);
%             end
%             if Z(j)<minASFcluster(i)
%                 minASFcluster(i) = Z(j); 
%             end
%         end
%         %{
%         for j = 1:size(opt.normalizedObj,1)
%             %obj = opt.normalizedObj(j,:);
%             %Z = max((obj)./w) + opt.rho*sum(obj,2);
%             obj = opt.archiveObj(j,:);
%             w_org(1) =  opt.min_val(1)+ w(1)*(opt.max_val(1)-opt.min_val(1));
%             w_org(2) =  opt.min_val(2)+ w(2)*(opt.max_val(2)-opt.min_val(2));
%             Z = max((obj)./w_org) + opt.rho*sum(obj,2);
%             %Z = max((obj-0.05)./w);%+opt.rho*sum(obj-0.05./w,2);%max(obj./w);%max((obj-opt.utopian)./w);
%             %Z = max(obj-(w-0.5)) ;%+ 0.5*norm(obj./norm(obj),0.2); %0.5*(sum(obj-(w-0.5)));
%             if(Z<S(j))
%                 cluster(j) = i;
%                 S(j) = Z;
%             end
%             if Z<minASFcluster(i)
%                 minASFcluster(i) = Z; 
%             end
%         end
%         %}
%         %S = bsxfun(@min,Z,S);
%     end

