% Generate image data info. Note that image files with sample code not 
% equal to '01' (means solid tumor) are excluded. And we also excluded 
% images with 5x magnification. One patient may have multiple slides (1 to 3).

% ==image data info
% imageInfo.file: list(i).name is 'folder/imageFilename'. Files are ranked
%   according to the order of the first 15 letters of the barcode
% imageInfo.pid: cell array of strings containing patient ID of each image 
%   file, the first 15 letters of barcode.
% imageInfo.mppX: 
% imageInfo.mppY:
% imageInfo.width:
% imageInfo.height:
% imageInfo.mag: manification

clear
tic

%% generate imageInfo.file, imageInfo.pid
dirData = 'G:/LUSC/WSI-jiankang/';
list = dir(dirData);
list = list(3:end);

nList = numel(list);
imageInfo.file = cell(nList, 1);
imageInfo.pid = cell(nList, 1);

tissue = cell(nList, 1);
for i = 1:nList
    sublist = dir([dirData, list(i).name, '/*.svs']);
    tissue{i} = sublist.name(18:19);            %�ĳ�ʵ�ʵ�λ��
    imageInfo.file{i} = [list(i).name, '/', sublist.name];
    imageInfo.pid{i} = sublist.name(1:15);      %�ĳ�ʵ�ʵ�λ��
    fprintf('%d/%d finished\n', i, nList);
end

% remove non-solid tumor files (sample code ~= 01) �Ƴ���ʵ�������ļ�(��������~= 01)
% indTumor = strcmp(tissue, '01');        %ֻҪһ��������ͼƬ
% imageInfo.file = imageInfo.file(indTumor);
% imageInfo.pid = imageInfo.pid(indTumor);

% sort
[imageInfo.pid, ind] = sort(imageInfo.pid);
imageInfo.file = imageInfo.file(ind);
% fprintf('%s finished\n', imageInfo.file{330});

%% generate imageInfo.mppX, imageInfo.mppY, imageInfo.width
%  imageInfo.height, and imageInfo.mag

% Load openslide library
openslide_load_library();

nFiles = numel(imageInfo.file);
disp(nFiles)
imageInfo.mppX = zeros(nFiles, 1);
imageInfo.mppY = zeros(nFiles, 1);
imageInfo.width = zeros(nFiles, 1);
imageInfo.height = zeros(nFiles, 1);
imageInfo.mag = zeros(nFiles, 1);
for i = 1:nFiles    
    fprintf('%d %s examed\n',i,imageInfo.file{i} );
    % Open whole-slide image
    slidePtr = openslide_open([dirData, imageInfo.file{i}]);  %��ͼƬ

    % Get whole-slide image properties
    [imageInfo.mppX(i, 1),...
        imageInfo.mppY(i, 1),...
        imageInfo.width(i, 1),...
        imageInfo.height(i, 1),...
        numberOfLevels,...
        downsampleFactors,...
        imageInfo.mag(i, 1)] = openslide_get_slide_properties(slidePtr);
    fprintf('%d/%d finished\n', i, nFiles);
end

toc

% Close whole-slide image, note that the slidePtr must be removed manually �ر����Żõ�Ƭͼ��ע������ֶ��Ƴ��õ�Ƭ
openslide_close(slidePtr)
clear slidePtr

% Unload library
openslide_unload_library

imageInfo = struct2table(imageInfo);
% remove images with 5x magnification �Ƴ��Ŵ�5����ͼ��  �����Ϊʲô
indOk = imageInfo.mag ~= 5;
imageInfo = imageInfo(indOk, :);

save imageInfo-jiankang imageInfo
