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
% Contact: royprote@msu.edu, proteek_buet@yahoo.com


function [selectionFunctionValue] = selectionFunction(opt, pop)

    switch(opt.selection_function_option)

        case 1 %KKTPM
            kktpm  = evaluate_kktpm(opt, pop);
            opt.archiveKKTPM = vertcat(opt.archiveKKTPM, kktpm);%archive kktpm
            selectionFunctionValue = opt.archiveKKTPM;
            
        case 2 %MEMO approach, similar to case 2, only code is different, without calculating distance to asf directions
            asfcv = cell2mat(opt.archiveASFCVAll);
            selectionFunctionValue = min(asfcv,[],2);

        case 3 %NSGA-II selection
            
            [R,~] = bos(opt.archiveObj);
            opt.Rank = zeros(size(opt.archiveObj,1),1);
            opt.Epsilon = 1e-14;
            opt.Inf = 1e14;
            for i=1:max(R)
                front = find(R==i)';
                front_cd = crowdingDistance(opt, front, opt.archiveObj(front,:));
                [~,Idx] = sort(front_cd,'descend');
                front = front(Idx);
                norm_cd_rank = ((1:size(front,2))/(size(front,2)));
                norm_cd_rank = norm_cd_rank - 0.000001;
                norm_cd_rank = norm_cd_rank/1000;
                front_rank = i + norm_cd_rank;
                opt.Rank(front) = front_rank;
            end

            selectionFunctionValue = opt.Rank;
            %--------------find cluster------------------------------------
            
            
            %%{
            S = zeros(size(opt.normalizedObj,1),1);
            %Z = zeros(size(opt.normalizedObj,1),1);
            for i=1:size(S,1)
                S(i) = intmax;
            end
            
            cluster = ones(size(opt.normalizedObj,1),1);
                        
            for i=1:size(opt.dirs,1)
                w = opt.dirs(i,:);
                w = w./norm(w);
                
                
                for j = 1:size(opt.normalizedObj,1)
                    obj = opt.normalizedObj(j,:);
                    Z = max((obj)./w);% + opt.rho*sum(obj,2);
                    if(Z<S(j))
                        cluster(j) = i;
                        S(j) = Z;
                    end
                end
            end
            
        
        otherwise
            input('Inside selectionFunction, sf not defined');
    end
    
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

