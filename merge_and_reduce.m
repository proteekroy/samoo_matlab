function [opt] = merge_and_reduce(opt)

% [X, Y, TotConsViol, MinASFs, clusterID]

k = opt.V;
M = opt.M;
p = opt.numdir-1;
ncon = opt.C;
Np = opt.numdir;
%pop3 = [pop1(:, 1:k+M+1); pop2(:, 1:k+M+1)];
pop3 = horzcat(opt.nsga2.totalpop, opt.nsga2.totalpopObj, opt.nsga2.totalpopCV);

Nfit = opt.nsga2.N;

%======================FIND UNIQUE=========================================
[~, index] = unique(pop3(:, k+1:k+M+1), 'rows'); % Compare in true objective space
new_index = sort(index);
pop3 = pop3( new_index, : );
opt.nsga2.totalpop =  opt.nsga2.totalpop(index,:);
opt.nsga2.totalpopObj = opt.nsga2.totalpopObj(index,:);
opt.nsga2.totalpopCV = opt.nsga2.totalpopCV(index,:);
opt.nsga2.totalpopCons = opt.nsga2.totalpopCons(index,:);
[~,opt.nsga2.CD] = apply_nsga2_selection(opt, opt.nsga2.totalpop, opt.nsga2.totalpopObj, opt.nsga2.totalpopCV, opt.nsga2.N);

%=====================FIND CLUSTER AND ASF IN COMBINED POP=================
[opt, pop3, cluster_sizes] = memo_clustering(opt, pop3);

%pop3 = pop3(1:end-2,:);

cluster_sizes2 = cluster_sizes;
zeroIndex = cluster_sizes2 == 0;
numberOfZeros = sum(zeroIndex);
cluster_sizes2(zeroIndex) = Inf;
[~, sortedIndex] = sort(cluster_sizes2);
sortedIndex = sortedIndex(1 : Np-numberOfZeros);
C = cell(Np-numberOfZeros, 1);
for i=1:Np-numberOfZeros
    if ncon == 0
        clusterMinASFs = find( pop3(:,end) == sortedIndex(i) );
        [~, temp_idx] = sort( pop3(clusterMinASFs, end-1) );
        C{i} = clusterMinASFs(temp_idx)';
    else % ncon > 0
        members = find( pop3(:,end) == sortedIndex(i) );
        localSort = sortClusterMembers( pop3(members,:), k, M );
        C{i} = members(localSort)';
    end
end


%---------PROTEEK----------------------------------------------------------
% allpop = linear_normalize(opt, opt.hifipop(:,1:opt.k), opt.Xmin, opt.Xmax);
% norm_pop = linear_normalize(opt, pop3(:,1:opt.k), opt.Xmin, opt.Xmax); 
% div = pdist2(norm_pop, allpop);%distance in normalized space
% div = min(div,[],2);
% I1 = find(div<=opt.TrustDistVar);

% if size(I1,1)>0 && size(I1,1)<=Nfit
%     count = size(I1,1);
%     if count>0
%         %disp('dd');
%         %pop(1:count,:) = pop3(I1, :);
%         %pop3 = pop3(setdiff(1:size(pop3,1), I1), :);
%         selectedPopIndex = I1';
%         count=count+1;
%     end
%     
% else
    count = 1;
    selectedPopIndex = [];
%end
%---------PROTEEK----------------------------------------------------------
%count = 1; 
j = 0;
while count <= Nfit
    j = j + 1;
    for i=1:Np-numberOfZeros
        while j > size( C{i},2 )
            i = i+1;
            if i > Np-numberOfZeros
                break
            end
        end
        %pop(count,:) = pop3( C{i}(1,j), : );
        if ~ismember(C{i}(1,j),selectedPopIndex)
            selectedPopIndex = horzcat(selectedPopIndex, C{i}(1,j));
            count = count + 1;
        end
        if count > Nfit
            break
        end
    end
end


%---------PROTEEK----------------------------------------------------------
% if max(selectedPopIndex)>size(opt.nsga2.totalpop,1)
%    disp('Wrong'); 
% end
opt.nsga2.pop =  opt.nsga2.totalpop(selectedPopIndex,:);
opt.nsga2.popObj = opt.nsga2.totalpopObj(selectedPopIndex,:);
opt.nsga2.popCV = opt.nsga2.totalpopCV(selectedPopIndex,:);
opt.nsga2.popCons = opt.nsga2.totalpopCons(selectedPopIndex,:);
opt.nsga2.CD = opt.nsga2.CD(selectedPopIndex,:);
    
%---------PROTEEK----------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function localSort = sortClusterMembers(pop, k, M)

N_total = size(pop,1);
localSort = [];

feasID = find(pop(:, k+M+1) == 0); N_feas = size(feasID,1);
if N_feas > 0
    [~, idx1] = sort( pop( feasID, k+M+2) ); % Sort with respect to MinASF value
    localSort = [localSort; feasID(idx1)];
end
if N_total > N_feas
    infeasID = find(pop(:, k+M+1) ~= 0);
    [~, idx2] = sort( pop( infeasID, k+M+2) ); % Sort with respect to CV value
    localSort = [localSort; infeasID(idx2)];
end




