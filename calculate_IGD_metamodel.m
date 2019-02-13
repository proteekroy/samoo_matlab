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


clear all;


problem = 'tnk';
method = 7;

r = 1;
%for r = 1:10
    h = figure;
    hold all;
    %datafilename = strcat('IGD Calculation/methodology',num2str(method),'_',problem,'_obj_',num2str(r),'.txt');
    %cvfilename = strcat('IGD Calculation/methodology',num2str(method),'_',problem,'_cv_',num2str(r),'.txt');
    datafilename = strcat('methodology',num2str(method),'_',problem,'_obj_',num2str(r),'.txt');
    cvfilename = strcat('methodology',num2str(method),'_',problem,'_cv_',num2str(r),'.txt');
    data = load(datafilename);
    cv =  load(cvfilename);

    index = cv<=0.0;
    data = data(index==1,:);

    index = paretofront(data);
    data = data(index==1,:);

    if size(data,2)==2
        str1 = strcat('../Metamodel Classification/Pareto Front Reference Direction/',upper(problem),'.2D.pf');
        PFStar = load(str1);
        if(size(PFStar,1))>21
           disp('Wrong Pareto front'); 
        end
    elseif size(data,2)==3
        str1 = strcat('../Metamodel Classification/Pareto Front Reference Direction/',upper(problem),'.3D.pf');
        %str2 = 'CAR-PF.txt';
        PFStar = load(str1);
        if(size(PFStar,1))>91
           disp('Wrong Pareto front'); 
        end
    end
    
    
    igd = IGD(PFStar', data');
    disp(igd);
    gd = IGD(data', PFStar');
    disp(gd);
    PFStar = load(str1);
    if size(data,2)==2
        plot(PFStar(:,1),PFStar(:,2),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',5); 
        plot(data(:,1),data(:,2),'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',5); 
    else
        plot3(PFStar(:,1),PFStar(:,2),PFStar(:,3),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',1);
        plot3(data(:,1),data(:,2),data(:,3),'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',5);
    end

    memo_leg = strcat('M',num2str(method),', IGD = ',num2str(igd));
    title(strcat('M',num2str(method),' run-',num2str(r)));
    lg = legend('Pareto Front',memo_leg);
    name = strcat('methodology',num2str(method), '_',problem,'_',num2str(r),'.fig');
    %saveas(h,name) 
%end