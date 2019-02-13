function [opt] = unique_population(opt)

    
    %-------------FIND UNIQUE SOLUTIONS--------------------
%     [~,ia,~] = unique(opt.archive,'rows');
%     opt.archive = opt.archive(ia,:);
%     opt.archiveObj = opt.archiveObj(ia,:);
%     opt.archiveASF = opt.archiveASF(ia,:);
%     opt.archiveCV = opt.archiveCV(ia,:); 
%     opt.archiveCons = opt.archiveCons(ia,:);
%     opt.archiveCluster = opt.archiveCluster(ia,:);
%     opt.normalizedObj = opt.normalizedObj(ia,:);
%     opt.archiveKKTPM = opt.archiveKKTPM(ia,:);
    
    
    %-------------ACTIVE ARCHIVE SET---------------------------------------
    
    [~,ia,~] = unique(opt.activeArchive,'rows');
    opt.activeArchive = opt.activeArchive(ia,:);
    opt.activeArchiveObj = opt.activeArchiveObj(ia,:);
    opt.activeArchiveASF = opt.activeArchiveASF(ia,:);
    opt.activeArchiveCons = opt.activeArchiveCons(ia,:);
    opt.activeArchiveCV = opt.activeArchiveCV(ia,:);
    
end