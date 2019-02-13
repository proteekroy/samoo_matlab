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

function NormalizedObj =  nsga3_normalization(opt, popObj, ideal_point, nadir_point)

    max_popObj = max(popObj);
    
    if ~isempty(opt.FeasibleIndex)
        combinedObj =  opt.archiveObj(opt.FeasibleIndex, :);
    else
        combinedObj = popObj;%vertcat();
    end
    %combinedObj = popObj;
    
    TranslatedObj = combinedObj - repmat(ideal_point, size(combinedObj, 1), 1);%recalculate objective according to z
            
    ASFLines = (1e-16)*ones(opt.M, opt.M)+ eye(opt.M);%ASF direction inclined along each objective
    S = zeros(opt.M, opt.M);%To collect the intercept
            
    for i=1:opt.M
        w = ASFLines(i,:);
        w = w./norm(w);
        [~, index] = min(max(TranslatedObj./repmat(w,size(TranslatedObj,1),1),[],2));%finding ASF values with a direction inclined in an objective
        S(i,:) = TranslatedObj(index,:); %choise the element with smallest asf w.r.t. ASFLines(i,:)
    end

    %----Check if M points doesn't make M-simplex----------------------
    %----It can also happen when size(lastPopIndex,2)<opt.M------------
        
    if det(S)<1e-16
        A = nadir_point; 
    else
        b = ones(opt.M,1);
        A = linsolve(S,b);
        A = 1./A;
        A = A';
        A(A<1e-16) = opt.nadir_point(A<1e-16);
            
        if ~isempty(max_popObj)
            A(A<1e-16) = max_popObj(A<1e-6);
        end
    end

        
    %------------------NORMALIZE WITH INTERCEPT------------------------
                
    NormalizedObj = TranslatedObj./repmat(A, size(TranslatedObj, 1),1);
           



end



% function NormalizedObj =  nsga3_normalization(opt, popObj)
% 
%     if ~isempty(opt.FeasibleIndex)
%         combinedObj =  opt.archiveObj(opt.FeasibleIndex, :);
%     else
%         combinedObj = popObj;%vertcat();
%     end
%     z = min(combinedObj);%find the population minimum
%     TranslatedObj = combinedObj - repmat(z, size(combinedObj, 1), 1);%recalculate objective according to z
%             
%     ASFLines = eye(opt.M);%ASF direction inclined along each objective
%     S = zeros(opt.M, opt.M);%To collect the intercept
%             
%     for i=1:opt.M
%         w = ASFLines(i,:);
%         w(1:opt.M~=i) = 1e-16;
%         w = w./norm(w);
%         [~, index] = min(max(TranslatedObj./repmat(w,size(TranslatedObj,1),1),[],2));%finding ASF values with a direction inclined in an objective
%         S(i,:) = TranslatedObj(index,:); %choise the element with smallest asf w.r.t. ASFLines(i,:)
%     end
% 
%     %----Check if M points doesn't make M-simplex----------------------
%     %----It can also happen when size(lastPopIndex,2)<opt.M------------
%         
%     if det(S)<1e-16   
%         A = max(TranslatedObj,[], 1);
%     else
%         b = ones(opt.M,1);
%         A = linsolve(S,b);
%         A = 1./A;
%         A = A';
%     end
% 
%         
%     %------------------NORMALIZE WITH INTERCEPT------------------------
%                 
%     NormalizedObj = popObj./repmat(A, size(popObj, 1),1);
%            
% 
% 
% 
% end