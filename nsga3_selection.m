function opt = nsga3_selection(opt)

    [n, ~] = size(opt.totalpopObj);
    opt.R = zeros(size(opt.totalpopObj,1),1);   
     
    selectedPopIndex = [];
    index = opt.totalpopCV<=0;
    FeasiblePopIndex = find(index == 1);    
    InfeasiblePopIndex = find(index == 0);

    %----------------------------------------------------------------------
    %---------Find Non-dominated Sorting of Feasible Solutions If exists---
    %----------------------------------------------------------------------
    
    
    F = cell(n,1);
    opt.R = zeros(n,1);
    if ~isempty(FeasiblePopIndex)
                
        R = bos(opt.totalpopObj(FeasiblePopIndex,:));
        
        for i=1:size(FeasiblePopIndex,1)
            F{R(i)} = horzcat(F{R(i)}, FeasiblePopIndex(i));
        end
    end 
    
    %----------------------------------------------------------------------       
    %----------------------CALCULATE NADIR POINT---------------------------
    %----------------------------------------------------------------------
    
    opt.nadir_point =  max(opt.totalpopObj(F{1},:));
    
    %----------------------------------------------------------------------
    %--If Less #Feasible <N, then we need sort infeasibles to select-------
    %----------------------Don't do niching here---------------------------
    %----------------------------------------------------------------------
    
    % Haitham - Start
    opt.pop2Dir = zeros(opt.N, 1);
    opt.pop2DirDistances = zeros(opt.N, 1);
    opt.associationsReady = false;
    idx = 1;
    % Haitham
    
    if size(FeasiblePopIndex, 1) < opt.N
        
        %-------------Pick all feasible solutions--------------------------
        selectedPopIndex = FeasiblePopIndex';
        
        %--------------Rank the Infeasible Solutions-----------------------
    
        if ~isempty(InfeasiblePopIndex)

            CV = opt.totalpopCV(InfeasiblePopIndex);
            [~,index] = sort(CV,'ascend');
            
            for i = 1: size(index,1)
                selectedPopIndex = horzcat(selectedPopIndex, InfeasiblePopIndex(index(i)));
                if size(selectedPopIndex, 2) >= opt.N
                    break; 
                end
            end

        end
    
    %----------------------------------------------------------------------
    %--If there are more #Feasibles than we need to select, Do Niching-----
    %----------------------------------------------------------------------
    
    else
        
        % Haitham - Start
        opt.associationsReady = true;
        % Haitham - End
        
        %------------------Find Last Front---------------------------------
        
        count = zeros(n,1);
        for i=1:n
            count(i) = size(F{i},2);
        end

        cd = cumsum(count);
        p1 = find(opt.N<cd);
        lastfront = p1(1);
        
        %-----------Select upto the front before last front----------------
        
        for i=1:lastfront-1
            selectedPopIndex = horzcat(selectedPopIndex, F{i});
        end
        
        %do a check in to see if N reached
        
        lastPopIndex = F{lastfront};
        
        combinedObj = vertcat(opt.totalpopObj(selectedPopIndex,:), opt.totalpopObj(lastPopIndex,:));
        combinePopIndex = horzcat(selectedPopIndex, lastPopIndex);
        
        %------------------------------------------------------------------       
        %-----FIND MAXIUM VALUES ACCROSS ALL POP INCLUDING LAST FRONT------
        %------------------------------------------------------------------
        
        opt.max_popObj = max(opt.totalpopObj(selectedPopIndex,:));
        
        %------------------------------------------------------------------       
        %-----------------------UPDATE IDEAL POINT-------------------------
        %------------------------------------------------------------------
        
        temp_pop = vertcat(combinedObj, opt.ideal_point);
        opt.ideal_point = min(temp_pop);
        
        %------------------------------------------------------------------       
        %-----------------TRANSLATE TO IDEAL POINT-------------------------
        %------------------------------------------------------------------
        
        %z = min(combinedObj);%find the population minimum
        temp_combinedObj = vertcat(combinedObj, opt.extreme_point);
        temp_TranslatedObj = temp_combinedObj - repmat(opt.ideal_point, size(temp_combinedObj, 1), 1);%recalculate objective according to z
        TranslatedObj = combinedObj - repmat(opt.ideal_point, size(combinedObj, 1), 1);
        
        ASFLines = (1e-16)*ones(opt.M, opt.M)+ eye(opt.M);%eye(opt.M);%ASF direction inclined along each objective
        S = zeros(opt.M, opt.M);%To collect the intercept
            
        %------------------------------------------------------------------       
        %---------------------UPDATE EXTREME POINTS------------------------
        %------------------------------------------------------------------
        %------------------S is a matrix of extreme points-----------------
        
        for i=1:opt.M
            w = ASFLines(i,:);
            % w(1:opt.M~=i) = 1e-16;
            w = w./norm(w);
            [~, index] = min(max(temp_TranslatedObj./repmat(w, size(temp_TranslatedObj,1),1),[],2));%finding ASF values with a direction inclined in an objective
            S(i,:) = temp_TranslatedObj(index,:); %choose the element with smallest asf w.r.t. ASFLines(i,:)
        end
        
        
        
        %------------------------------------------------------------------       
        %-----------------------FIND INTERCEPT-----------------------------
        %------------------------------------------------------------------
        
        %----Check if M points doesn't make M-simplex----------------------
        %----It can also happen when size(lastPopIndex,2)<opt.M------------
        
        if det(S)<1e-16   
            %--------------TO-DO----------------
            A = opt.nadir_point; %max(TranslatedObj,[], 1);%take max of first front only as intercept
        else
            b = ones(opt.M,1);
            A = linsolve(S,b);
            A = 1./A;
            A = A';
            
            %if any(A<1e-6)
            A(A<0) = opt.nadir_point(A<0);
            
            if ~isempty(opt.max_popObj)
                A(A<1e-6) = opt.max_popObj(A<1e-6);
            end
            %end
            %--------------TO-DO----------------
            %if some intercept is negative replace with corresponding nadir
            %point value
            %if intercept is zero (<10^-6) then replace this with max of
            %population (translated obj)
            %--------------TO-DO----------------
%             if opt.M==3 % check intercept
%                 hold all;
%                 xlabel('x')
%                 ylabel('y')
%                 zlabel('z')
%                 box on
%                 grid on
%                 hold on
%                 Intercept = [A(1) 0 0;
%                             0 A(2) 0;
%                             0 0 A(3)];
%                 fill3(Intercept(:,1),Intercept(:,2),Intercept(:,3),'c');
%                 scatter3(S(:,1), S(:,2), S(:,3),'*');
%             end
        end
        
        %---------------Plot to see if intercept is correct----------------
       
        
        
        %------------------NORMALIZE WITH INTERCEPT------------------------
            
        NormalizedObj = TranslatedObj./repmat(A,size(TranslatedObj,1),1);
                
        %------------------------------------------------------------------    
        %-------------------ASSOCIATION------------------------------------
        %------------------------------------------------------------------
        
        %------------------------------------------------------------------
        %-----For Each Direction select nearest solution for last front----
        %-----For Two cases: Selected fronts, and Last front---------------
        %------------------------------------------------------------------
        %should have been same for all cases
        
        Count = zeros(opt.numdir, 1);
               
        %----------------for selected indices------------------------------
        
        tempPopDir = zeros(size(selectedPopIndex, 2), 1);%zeros(size(opt.totalpopObj, 1), 1);%for each solution save direction
        distance = zeros(size(selectedPopIndex, 2), 1);%zeros(size(opt.totalpopObj, 1), 1);%for each solution save distance of that direction
        distance(1:end) = inf;%intmax;
        
        for i = 1:size(selectedPopIndex, 2)
                    
            obj = NormalizedObj(i,:);%NormalizedObj(find(combinePopIndex==s),:);
            
            for j = 1:opt.numdir
                w = opt.dirs(j,:)./norm(opt.dirs(j,:));                                             
                d = sqrt(power(norm(obj), 2)-power(dot(obj,w),2));%norm(obj-((w*obj')*w)/(norm(w)*norm(w)));

                if d < distance(i)
                    tempPopDir(i) = j;
                    distance(i) = d;
                end
            end
            
            %Count(int8(tempPopDir(int8(s)))) = Count(int8(tempPopDir(int8(s))))+1;
            Count(tempPopDir(i)) = Count(tempPopDir(i))+1;

            % Haitham - Start
            opt.pop2Dir(idx) = tempPopDir(i);%tempPopDir(s);
            opt.pop2DirDistances(idx) = distance(i);%distance(s);
            idx = idx + 1;
            % Haitham - End

        end
        
        %------------------------------------------------------------------
%         figure(opt.fig2);
%         clf;
%         hold all;
%         
%         for i=1:opt.numdir
%             w = opt.dirs(i,:);
%             plot([0 w(1)], [0 w(2)], '-k', 'LineWidth', 1);
%         end
%         
%         plot(NormalizedObj (:,1), NormalizedObj (:,2), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
        
        %----------------for last front indices----------------------------
        DirLast = cell(opt.numdir, 1);%collect the solutions those prefer this direction
        DirLastDist = cell(opt.numdir, 1);%collect the distances
        tempPopDir = zeros(size(lastPopIndex, 2), 1);%save direction for 
        distance = zeros(size(lastPopIndex, 2), 1);%for each solution save distance of that direction
        distance(1:end) = inf;%intmax;
        
        for i = 1:size(lastPopIndex, 2)
            s = lastPopIndex(i);
            obj = NormalizedObj(i, :);%NormalizedObj(find(combinePopIndex==s),:);%take normalized obj
            
            for j = 1:opt.numdir
                w = opt.dirs(j,:)./norm(opt.dirs(j,:));
                d = sqrt(power(norm(obj), 2)-power(dot(obj, w),2));%norm(obj-((w*obj')*w)/(norm(w)*norm(w)));
                
                if d < distance(i)
                    tempPopDir(i) = j;
                    distance(i) = d;
                end
            end
            DirLast{tempPopDir(i)} = horzcat(DirLast{tempPopDir(i)}, s);
            DirLastDist{tempPopDir(i)} = horzcat(DirLastDist{tempPopDir(i)}, distance(i));
        end
        
        %-------------------NICHING PRESERVATION---------------------------
        
        for i = 1:opt.numdir
            [DirLastDist{i}, I] = sort(DirLastDist{i},'ascend');
            DirLast{i} = DirLast{i}(I);
        end
                                               
        while size(selectedPopIndex, 2)<opt.N

            [p_j, j] = min(Count);
            
            if isempty(DirLast{j}) 
                %continue;
                %Count(j)=intmax; %excluded for further consideration    
                Count(j)= Inf;
            elseif p_j==0                
                selectedPopIndex = horzcat(selectedPopIndex, DirLast{j}(1));%choose the first
                
                % Haitham - Start
                opt.pop2Dir(idx) = DirLast{j}(1);
                opt.pop2DirDistances(idx) = DirLastDist{j}(1);
                idx = idx + 1;
                % Haitham - End

                Count(j) = Count(j)+1;
                DirLast{j}(1) = [];
                DirLastDist{j}(1) = [];

            else % when p_j>=1                
                
                % NSGA-III
                index = randi(size(DirLast{j},2)); %chose randomly
                selectedPopIndex = horzcat(selectedPopIndex, DirLast{j}(index));

                % Haitham - Start
                opt.pop2Dir(idx) = DirLast{j}(index);
                opt.pop2DirDistances(idx) = DirLastDist{j}(index);
                idx = idx + 1;
                % Haitham - End

                DirLast{j}(index) = [];
                DirLastDist{j}(1) = [];
                Count(j) = Count(j)+1;     
                

            end                                   
        end
        
        %I = ismember(combinePopIndex, selectedPopIndex);
        %plot(NormalizedObj (I,1), NormalizedObj (I,2), 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
        
    end%if size(FeasiblePopIndex,1)<opt.N           
    
    %---------------Select for Next Generation-----------------------------
    opt.pop =  opt.totalpop(selectedPopIndex,:);
    opt.popObj = opt.totalpopObj(selectedPopIndex,:);
    opt.popCV = opt.totalpopCV(selectedPopIndex,:);
    opt.popCons = opt.totalpopCons(selectedPopIndex,:);
            
end

% % Copyright [2018] [Proteek Chandan Roy]
% % 
% % Licensed under the Apache License, Version 2.0 (the "License");
% % you may not use this file except in compliance with the License.
% % You may obtain a copy of the License at
% % 
% %     http://www.apache.org/licenses/LICENSE-2.0
% % 
% % Unless required by applicable law or agreed to in writing, software
% % distributed under the License is distributed on an "AS IS" BASIS,
% % WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% % See the License for the specific language governing permissions and
% % limitations under the License.
% 
% function opt = nsga3_selection(opt)
% 
%     [n, ~] = size(opt.totalpopObj);
%     opt.R = zeros(size(opt.totalpopObj,1),1);   
%      
%     selectedPopIndex = [];
%     index = opt.totalpopCV<=0;
%     FeasiblePopIndex = find(index == 1);    
%     InfeasiblePopIndex = find(index == 0);
% 
%     %----------------------------------------------------------------------
%     %---------Find Non-dominated Sorting of Feasible Solutions If exists---
%     %----------------------------------------------------------------------
%     
%     
%     F = cell(n,1);
%     opt.R = zeros(n,1);
%     if ~isempty(FeasiblePopIndex)
%                 
%         [R,~] = bos(opt.totalpopObj(FeasiblePopIndex,:));
%         
%         for i=1:size(FeasiblePopIndex,1)
%             F{R(i)} = horzcat(F{R(i)}, FeasiblePopIndex(i));
%         end
%     end 
%     
%     
%     %----------------------------------------------------------------------
%     %--If Less #Feasible <N, then we need sort infeasibles to select-------
%     %----------------------Don't do niching here---------------------------
%     %----------------------------------------------------------------------
%     
%     %% Haitham - Start
%     opt.pop2Dir = zeros(opt.N, 1);
%     opt.pop2DirDistances = zeros(opt.N, 1);
%     opt.associationsReady = false;
%     idx = 1;
%     % Haitham
%     
%     if size(FeasiblePopIndex,1)<opt.N
%         
%         %-------------Pick all feasible solutions--------------------------
%         selectedPopIndex = FeasiblePopIndex';
%         
%         %--------------Rank the Infeasible Solutions-----------------------
%     
%         if ~isempty(InfeasiblePopIndex)
% 
%             CV = opt.totalpopCV(InfeasiblePopIndex);
%             [~,index] = sort(CV,'ascend');
%             
%             for i = 1: size(index,1)
%                 selectedPopIndex = horzcat(selectedPopIndex, InfeasiblePopIndex(index(i)));
%                 if size(selectedPopIndex, 2) >= opt.N
%                     break; 
%                 end
%             end
% 
%         end
%     
%     %----------------------------------------------------------------------
%     %--If there are more #Feasibles than we need to select, Do Niching-----
%     %----------------------------------------------------------------------
%     
%     else
%         
%         %% Haitham - Start
%         opt.associationsReady = true;
%         % Haitham - End
%         
%         %------------------Find Last Front---------------------------------
%         
%         count = zeros(n,1);
%         for i=1:n
%             count(i) = size(F{i},2);
%         end
% 
%         cd = cumsum(count);
%         p1 = find(opt.N<cd);
%         lastfront = p1(1);
%         
%         %-----------Select upto the front before last front----------------
%         
%         for i=1:lastfront-1
%             selectedPopIndex = horzcat(selectedPopIndex, F{i});
%         end
%         
%         %do a check in to see if N reached
%         
%         lastPopIndex = F{lastfront};
%         
%         combinedObj = vertcat(opt.totalpopObj(selectedPopIndex,:), opt.totalpopObj(lastPopIndex,:));
%         combinePopIndex = horzcat(selectedPopIndex, lastPopIndex);
%                
%         %------------------------------------------------------------------       
%         %---------------------------FIND INTERCEPT-------------------------
%         %------------------------------------------------------------------
%             
%         z = min(combinedObj);%find the population minimum
%         TranslatedObj = combinedObj - repmat(z, size(combinedObj, 1), 1);%recalculate objective according to z
%             
%         ASFLines = eye(opt.M);%ASF direction inclined along each objective
%         S = zeros(opt.M, opt.M);%To collect the intercept
%             
%         for i=1:opt.M
%             w = ASFLines(i,:);
%             w(1:opt.M~=i) = 1e-16;
%             w = w./norm(w);
%             [~, index] = min(max(TranslatedObj./repmat(w,size(TranslatedObj,1),1),[],2));%finding ASF values with a direction inclined in an objective
%             S(i,:) = TranslatedObj(index,:); %choise the element with smallest asf w.r.t. ASFLines(i,:)
%         end
% 
%         %----Check if M points doesn't make M-simplex----------------------
%         %----It can also happen when size(lastPopIndex,2)<opt.M------------
%         
%         if det(S)<1e-16   
%             A = max(TranslatedObj,[], 1);
%         else
%             b = ones(opt.M,1);
%             A = linsolve(S,b);
%             A = 1./A;
%             A = A';
%             
% %             if opt.M==3 % check intercept
% %                 hold all;
% %                 xlabel('x')
% %                 ylabel('y')
% %                 zlabel('z')
% %                 box on
% %                 grid on
% %                 hold on
% %                 Intercept = [A(1) 0 0;
% %                             0 A(2) 0;
% %                             0 0 A(3)];
% %                 fill3(Intercept(:,1),Intercept(:,2),Intercept(:,3),'c');
% %                 scatter3(S(:,1), S(:,2), S(:,3),'*');
% %             end
%         end
%         
%         %---------------Plot to see if intercept is correct----------------
%        
%         
%         
%         %------------------NORMALIZE WITH INTERCEPT------------------------
%             
%         NormalizedObj = TranslatedObj./repmat(A,size(TranslatedObj,1),1);
%                 
%         %------------------------------------------------------------------    
%         %-------------------ASSOCIATION------------------------------------
%         %------------------------------------------------------------------
%         
%         %------------------------------------------------------------------
%         %-----For Each Direction select nearest solution for last front----
%         %-----For Two cases: Selected fronts, and Last front---------------
%         %------------------------------------------------------------------
%         %should have been same for all cases
%         
%         Count = zeros(opt.numdir, 1);
%                
%         %----------------for selected indices------------------------------
%         
%         tempPopDir = zeros(size(opt.totalpopObj, 1), 1);%for each solution save direction
%         distance = zeros(size(opt.totalpopObj, 1), 1);%for each solution save distance of that direction
%         distance(1:end) = intmax;
%         
%         for i = 1:size(selectedPopIndex, 2)
%             s = selectedPopIndex(i);        
%             obj = NormalizedObj(find(combinePopIndex==s),:);
%             
%             for j = 1:opt.numdir
%                 w = opt.dirs(j,:);                                              
%                 d = norm(obj-((w*obj')*w)/(norm(w)*norm(w)));
%                 
%                 if d < distance(s)
%                     tempPopDir(s) = j;
%                     distance(s) = d;
%                 end
%             end
%             Count(tempPopDir(s)) = Count(tempPopDir(s))+1;
% 
%             %% Haitham - Start
%             opt.pop2Dir(idx) = tempPopDir(s);
%             opt.pop2DirDistances(idx) = distance(s);
%             idx = idx + 1;
%             % Haitham - End
% 
%         end
%         
%         %----------------for last front indices----------------------------
%         DirLast = cell(opt.numdir,1);%collect the solutions those prefer this direction
%         DirLastDist = cell(opt.numdir,1);%collect the distances
%         tempPopDir = zeros(size(opt.totalpopObj, 1), 1);%for each solution save direction
%         distance = zeros(size(opt.totalpopObj, 1), 1);%for each solution save distance of that direction
%         distance(1:end) = intmax;
%         
%         for i = 1:size(lastPopIndex, 2)
%             s = lastPopIndex(i);
%             obj = NormalizedObj(find(combinePopIndex==s),:);%take normalized obj
%             
%             for j = 1:opt.numdir
%                 w = opt.dirs(j,:);
%                 d = norm(obj-((w*obj')*w)/(norm(w)*norm(w)));
%                 if d < distance(s)
%                     tempPopDir(s) = j;
%                     distance(s) = d;
%                 end
%             end
%             DirLast{tempPopDir(s)} = horzcat(DirLast{tempPopDir(s)}, s);
%             DirLastDist{tempPopDir(s)} = horzcat(DirLastDist{tempPopDir(s)}, distance(s));
%         end
%         
%         %-------------------NICHING PRESERVATION---------------------------
%         
%         for i = 1:opt.numdir
%             [~, I] = sort(DirLastDist{i},'ascend');
%             DirLast{i} = DirLast{i}(I);
%         end
%                                                
%         while size(selectedPopIndex, 2)<opt.N
% 
%             [p_j, j]=min(Count);
%             
%             if isempty(DirLast{j})
%                 Count(j)=intmax; %excluded for further consideration              
%             elseif p_j==0                
%                 selectedPopIndex = horzcat(selectedPopIndex, DirLast{j}(1));%choose the first
%                 
%                 %% Haitham - Start
%                 opt.pop2Dir(idx) = DirLast{j}(1);
%                 opt.pop2DirDistances(idx) = DirLastDist{j}(1);
%                 idx = idx + 1;
%                 % Haitham - End
% 
%                 Count(j) = Count(j)+1;
%                 DirLast{j}(1) = [];
% 
%             else % when p_j>=1                
%                 index = randi(size(DirLast{j},2)); %chose randomly
%                 selectedPopIndex = horzcat(selectedPopIndex, DirLast{j}(index));
% 
%                 %% Haitham - Start
%                 opt.pop2Dir(idx) = DirLast{j}(index);
%                 opt.pop2DirDistances(idx) = DirLastDist{j}(index);
%                 idx = idx + 1;
%                 % Haitham - End
% 
%                 DirLast{j}(index) = [];
%                 Count(j) = Count(j)+1;                 
% 
%             end                                   
%         end
%         
%         
%     end%if size(FeasiblePopIndex,1)<opt.N           
%     
%     %---------------Select for Next Generation-----------------------------
%     opt.pop =  opt.totalpop(selectedPopIndex,:);
%     opt.popObj = opt.totalpopObj(selectedPopIndex,:);
%     opt.popCV = opt.totalpopCV(selectedPopIndex,:);
%     opt.popCons = opt.totalpopCons(selectedPopIndex,:);
%             
% end
