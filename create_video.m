
shuttleVideo = VideoReader('shuttle.avi');
workingDir = pwd;
imageNames = dir(fullfile(workingDir,'video-m1-1','*.jpg'));
%imageNames = {imageNames.name}';

imageNamesNew = cell(1,length(imageNames));
% 
% for K = 1 : length(imageNames)
%    filename = imageNames(K).name;
%    c = strsplit(imageNames(K).name,'.');
%    d = strsplit(c{1},'_');
%    imageNamesNew{str2num(d{3})} = filename;
% end



imageNames = dir(fullfile(workingDir,'video-m1-1','*.jpg'));
imageNames = {imageNames.name}';

for i=1:159
    imageNamesNew{i} =  imageNames{219+i};
end
for i=160:378
    imageNamesNew{i} =  imageNames{i-159};
end
outputVideo = VideoWriter(fullfile(workingDir,'methodology_m11.avi'));
outputVideo.FrameRate = shuttleVideo.FrameRate;
open(outputVideo)


for ii = 1:length(imageNamesNew)
   img = imread(fullfile(workingDir,'video',imageNamesNew{ii}));
   writeVideo(outputVideo,img)
end

close(outputVideo);