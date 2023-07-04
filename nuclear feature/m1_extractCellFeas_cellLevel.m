% Extract 10 cell features for each image in 'imageInfo'
clear
tic

% nucleus area in real size    %��ʱ�����ʲô��˼
mpp_40 = 0.2525; % 40X
realMin = 80*mpp_40^2;
realMax = 1100*mpp_40^2;
realPSize = 2000*mpp_40;


% Load openslide library
openslide_load_library();

strc = load('imageInfo-jiankang.mat');
imInfo = strc.imageInfo;
dirData = 'G:/LUSC/WSI-jiankang/';

% process each image
nFile = size(imInfo, 1);
for i = 1:nFile
    t1 = tic;  %���浱ǰʱ��
    file = imInfo.file{i};
    parts = regexp(file, '/', 'split');
    filename = parts{2}(1:end-4);
    
    % test if the image is already processed
    if exist(['extractCellFeas_cellLevel_jiankang/', filename, '.mat'], 'file');        
        fprintf('image %d/%d already processed\n', i, nFile);
        continue;
    end
%     fprintf('%s,%s',dirData,file);
    slidePtr = openslide_open([dirData, file]);
    mpp = imInfo.mppX(i);
    width = imInfo.width(i);
    height = imInfo.height(i);
    ps = round(realPSize/mpp);
    
    %% select patches using thumnail  ʹ������ͼѡ�񲹶�
    % Get thumbnail
    thumnail = openslide_read_associated_image(slidePtr, 'thumbnail');
    thumnail = thumnail(:, :, 2:4);
    [ht, wt, ~] = size(thumnail);
    ratio = width/wt;         %ԭͼ������ͼ֮��
    
    % patch size in thumbnail  ����ͼ�еĲ�����С
    pst = round(ps/ratio);     %�õ�����ͼ��С

    % select patches of interest in thumbnail image  ������ͼ��ѡ�����Ȥ�Ĳ���
    [X, Y] = meshgrid(1:pst:wt-pst+1, 1:pst:ht-pst+1);
    %meshgrid�������ɵ�X��Y�Ǵ�С��ȵľ���xgv��ygv����������ʸ����xgv��ygv������������
    %X��ͨ����xgv����length(ygv)�У��ϸ���������length(ygv)-1�У��õ�
    %Y�����ȶ�ygv����ת�õõ�ygv'����ygv'���ƣ�length(xgv)-1���εõ���
    xy = [X(:), Y(:)];   %���ÿһ������ϵ����
    d1 = size(xy, 1);    %�ܹ����ٸ��� Ϊÿһ������ϵ�
    indOk = zeros(d1, 1);
    for j = 1:d1
        r = xy(j, 2);
        c = xy(j, 1);
        rows = r:r+pst-1;
        cols = c:c+pst-1;
        tile = thumnail(rows, cols, :);   %��ȡ����һ��
        tMean = mean(tile, 3);            %3��ͨ��ȡƽ��
        if sum(tMean(:)>210) < pst*pst/2  %ѡȡ��һ�� ԭ����
            indOk(j) = 1;
        end
    end
    xy = xy(indOk==1, :)-1; % move upper-left point to (0, 0) �����ϵ��ƶ���(0��0)
    xy = round(xy*ratio);   %������ͼ��ԭ
    
    % remvoe coordinates exceeding width or height �Ƴ����곬����Ȼ�߶�
    indOk = xy(:, 1)<=width-ps & xy(:, 2)<=height-ps;
    xy = xy(indOk, :);
    
    %% extract cell features for each tile  ��ȡÿ����Ƭ�ĺ˼�����
    d1xy = size(xy, 1);            %���ս�ȡС�������
    mypitchCount=d1xy;
    cellFeas = cell(d1xy, 1);
    for ixy = 1:d1xy
        x = xy(ixy, 1);
        y = xy(ixy, 2);
%         tile = openslide_read_region(slidePtr, x, y, ps, ps, 0);    %��ȡһС��  RGBA
        try
            tile = openslide_read_region(slidePtr, x, y, ps, ps, 0);    %��ȡһС��  RGBA
        catch
            mypitchCount=mypitchCount-1;
            fprintf('pitch index %d,and the total count %d',ixy,mypitchCount);
            continue;
        end
        tile = tile(:, :, 2:4);  %RGB
        areaMin = round(realMin/mpp^2);
        bw = hmt(tile, areaMin);    %�ָ�ϸ����

        statsR = regionprops('table', bw, tile(:, :, 1), 'area', 'MeanIntensity',...
            'MajorAxisLength', 'MinorAxisLength', 'centroid');  
        statsR.Area = statsR.Area*mpp^2;   %ϸ���˴�С
        statsR.MajorAxisLength = statsR.MajorAxisLength*mpp;   %���᳤��
        statsR.MinorAxisLength = statsR.MinorAxisLength*mpp;   %���᳤��
        
        indOk = statsR.Area <= realMax;
        statsR = statsR(indOk, :);
        if size(statsR, 1) < 80    %�������80��ϸ���˾���һ��
            continue;
        end
        
        statsG = regionprops('table', bw, tile(:, :, 2), 'MeanIntensity');
        statsB = regionprops('table', bw, tile(:, :, 3), 'MeanIntensity');  %Bͨ����ƽ����ɫ����
        
        statsG = statsG(indOk, :);
        statsB = statsB(indOk, :);

        feas1_7 = [statsR.Area, statsR.MajorAxisLength, statsR.MinorAxisLength,...
            statsR.MajorAxisLength ./ statsR.MinorAxisLength,...
            statsR.MeanIntensity, statsG.MeanIntensity, statsB.MeanIntensity];


        % last 3 features derived from Delaunay graph ��Delaunayͼ���������3������
        centroids = statsR.Centroid;   %Centroid����
        feas8_10 = zeros(size(statsR, 1), 3);
        DT = delaunayTriangulation(double(centroids));
        E = edges(DT);
        for k = 1:size(centroids, 1)  %����ÿһ�����ģ�ϸ���ˣ�
            edgesCell = E(sum(E==k, 2)~=0, :);   %�õ����ϸ���˵����б�
            dist = zeros(size(edgesCell, 1), 1);  %��������
            for m = 1:numel(dist)                %����ÿһ����
                p1 = centroids(edgesCell(m, 1), :);  %һ�˵������
                p2 = centroids(edgesCell(m, 2), :);  %��һ�˵������
                dist(m) = norm(p1-p2);             %�������
            end
            feas8_10(k, :) = [mean(dist), max(dist), min(dist)];
        end

        cellFeas{ixy, 1} = [feas1_7, feas8_10*mpp];        
    end
        
    cellFeas(cellfun(@isempty, cellFeas)) = []; % remove empty elements
    save(['extractCellFeas_cellLevel_jiankang', '/', filename, '.mat'], 'cellFeas');
    fprintf('image %d/%d processed, time %f\n', i, nFile, toc(t1));
end


toc

% Close whole-slide image, note that the slidePtr must be removed manually
openslide_close(slidePtr)
clear slidePtr

% Unload library
openslide_unload_library