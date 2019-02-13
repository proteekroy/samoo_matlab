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
% Contact: proteek_buet@yahoo.com

%==========================================================================
%======================BUILD ALL MODELS====================================
%==========================================================================

function opt = collect_models(opt, TrainIndex)

    x = opt.archive(TrainIndex,:);
    f = opt.archiveObj(TrainIndex, :);
    normalized_f = opt.normalizedObj(TrainIndex, :);
    cons = opt.archiveCons(TrainIndex,:);
    cv = opt.archiveCV(TrainIndex,:); %archive constrain violation
    acv = opt.archiveACV(TrainIndex,:);
    y_asf = cell2mat(opt.archiveASFAll);
    y_asfcv = cell2mat(opt.archiveASFCVAll);
    y_asf = y_asf(TrainIndex,:);
    y_asfcv = y_asfcv(TrainIndex,:);
    
    [~,ia,~] = unique(x,'rows');
    x = x(ia,:);
    cons = cons(ia,:);
    f = f(ia,:);
    SF_VAL = opt.archiveSF_VAL(TrainIndex,:);
    SF_VAL = SF_VAL(ia,:);
    normalized_f = normalized_f(ia, :);
    cv = cv(ia,:);
    acv = acv(ia, :);
    y_asf = y_asf(ia,:);
    y_asfcv = y_asfcv(ia,:);
    
    
    n = size(x,1);
   
    %----------------MODEL OBJECTIVES--------------------------------------
    
    [opt.dmodel_obj, ~] = dacefit(x, f, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model obj
    [opt.dmodel_normalized_obj, ~] = dacefit(x, normalized_f, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model obj
    
    %----------------MODEL CONSTRAINTS AND CV------------------------------
    
    if opt.C>0
        [opt.dmodel_cons, ~] = dacefit(x, cons, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model cv
        [opt.dmodel_cv, ~] = dacefit(x, acv, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model cv
    end
    
    %---------------MODEL ASF FOR EACH DIRECTIONS--------------------------
    
    [opt.dmodel_asf, ~] = dacefit(x, y_asf, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model asf       
    [opt.dmodel_asfcv, ~] = dacefit(x, y_asfcv, @regpoly2, @corrgauss, opt.theta, opt.lob, opt.upb);%model asfcv
    
    
    %---------------MODEL SELECTION FUNCTION-------------------------------
    opt.net = train_neural_network(x, SF_VAL);
    
end

