% Extract 10 cell features for each image in 'imageInfo'
clear
tic

% nucleus area in real size    %暂时不清楚什么意思
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
    t1 = tic;  %保存当前时间
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
    
    %% select patches using thumnail  使用缩略图选择补丁
    % Get thumbnail
    thumnail = openslide_read_associated_image(slidePtr, 'thumbnail');
    thumnail = thumnail(:, :, 2:4);
    [ht, wt, ~] = size(thumnail);
    ratio = width/wt;         %原图与略缩图之比
    
    % patch size in thumbnail  缩略图中的补丁大小
    pst = round(ps/ratio);     %得到略缩图大小

    % select patches of interest in thumbnail image  在缩略图中选择感兴趣的补丁
    [X, Y] = meshgrid(1:pst:wt-pst+1, 1:pst:ht-pst+1);
    %meshgrid函数生成的X，Y是大小相等的矩阵，xgv，ygv是两个网格矢量，xgv，ygv都是行向量。
    %X：通过将xgv复制length(ygv)行（严格意义上是length(ygv)-1行）得到
    %Y：首先对ygv进行转置得到ygv'，将ygv'复制（length(xgv)-1）次得到。
    xy = [X(:), Y(:)];   %存放每一块的左上点矩阵
    d1 = size(xy, 1);    %总共多少个点 为每一块的左上点
    indOk = zeros(d1, 1);
    for j = 1:d1
        r = xy(j, 2);
        c = xy(j, 1);
        rows = r:r+pst-1;
        cols = c:c+pst-1;
        tile = thumnail(rows, cols, :);   %截取其中一块
        tMean = mean(tile, 3);            %3个通道取平均
        if sum(tMean(:)>210) < pst*pst/2  %选取这一块 原因不清
            indOk(j) = 1;
        end
    end
    xy = xy(indOk==1, :)-1; % move upper-left point to (0, 0) 将左上点移动到(0，0)
    xy = round(xy*ratio);   %将略缩图还原
    
    % remvoe coordinates exceeding width or height 移除坐标超过宽度或高度
    indOk = xy(:, 1)<=width-ps & xy(:, 2)<=height-ps;
    xy = xy(indOk, :);
    
    %% extract cell features for each tile  提取每个切片的核级特征
    d1xy = size(xy, 1);            %最终截取小块的数量
    mypitchCount=d1xy;
    cellFeas = cell(d1xy, 1);
    for ixy = 1:d1xy
        x = xy(ixy, 1);
        y = xy(ixy, 2);
%         tile = openslide_read_region(slidePtr, x, y, ps, ps, 0);    %读取一小块  RGBA
        try
            tile = openslide_read_region(slidePtr, x, y, ps, ps, 0);    %读取一小块  RGBA
        catch
            mypitchCount=mypitchCount-1;
            fprintf('pitch index %d,and the total count %d',ixy,mypitchCount);
            continue;
        end
        tile = tile(:, :, 2:4);  %RGB
        areaMin = round(realMin/mpp^2);
        bw = hmt(tile, areaMin);    %分割细胞核

        statsR = regionprops('table', bw, tile(:, :, 1), 'area', 'MeanIntensity',...
            'MajorAxisLength', 'MinorAxisLength', 'centroid');  
        statsR.Area = statsR.Area*mpp^2;   %细胞核大小
        statsR.MajorAxisLength = statsR.MajorAxisLength*mpp;   %长轴长度
        statsR.MinorAxisLength = statsR.MinorAxisLength*mpp;   %短轴长度
        
        indOk = statsR.Area <= realMax;
        statsR = statsR(indOk, :);
        if size(statsR, 1) < 80    %如果少于80个细胞核就下一个
            continue;
        end
        
        statsG = regionprops('table', bw, tile(:, :, 2), 'MeanIntensity');
        statsB = regionprops('table', bw, tile(:, :, 3), 'MeanIntensity');  %B通道的平均颜色像素
        
        statsG = statsG(indOk, :);
        statsB = statsB(indOk, :);

        feas1_7 = [statsR.Area, statsR.MajorAxisLength, statsR.MinorAxisLength,...
            statsR.MajorAxisLength ./ statsR.MinorAxisLength,...
            statsR.MeanIntensity, statsG.MeanIntensity, statsB.MeanIntensity];


        % last 3 features derived from Delaunay graph 从Delaunay图导出的最后3个特征
        centroids = statsR.Centroid;   %Centroid质心
        feas8_10 = zeros(size(statsR, 1), 3);
        DT = delaunayTriangulation(double(centroids));
        E = edges(DT);
        for k = 1:size(centroids, 1)  %对于每一个质心（细胞核）
            edgesCell = E(sum(E==k, 2)~=0, :);   %得到这个细胞核的所有边
            dist = zeros(size(edgesCell, 1), 1);  %计算数量
            for m = 1:numel(dist)                %对于每一条边
                p1 = centroids(edgesCell(m, 1), :);  %一端点的坐标
                p2 = centroids(edgesCell(m, 2), :);  %另一端点的坐标
                dist(m) = norm(p1-p2);             %计算距离
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