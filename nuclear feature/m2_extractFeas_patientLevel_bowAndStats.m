% Each cell has 10 features. Below we get patient-level feaures of each 
% patient.

% ---BoW histogram feature--- 
% For each type of cell feature, a ten-dimensional BoW histogram feature
% is generated for each patient. So a patient will be described by 10*10
% features in total.

% ---statistics of distribution---
% For each type of feature, five statistics are computed. So there will be 
% 50 features in total.

clear

dir1 = 'extractCellFeas_cellLevel_jiankang/';
list = dir([dir1, '*.mat']);
dSize = 10; %聚类为10个簇

%% clustering each type of feature 聚类每种类型的要素
% If dictionary exists, load it
if exist(['dic', num2str(dSize), 'jiankang.mat'], 'file')    %如果之前有做过这个图像的聚类则读取
    strc = load(['dic', num2str(dSize), 'jiankang.mat']);
    dic = strc.dic;
else % else do clustering
    n = 136;     %N=300  之前 一定要改回来！！！
    rng('default');    %控制随机数的生成
    indices = randperm(numel(list), n);  %返回从1到numel(list)的随机n个数，这n个数不相同

    feas = cell(n, 1);
    for i = 1:numel(indices)
        strc = load([dir1, list(indices(i)).name]);
        feas{i} = cell2mat(strc.cellFeas);       %随机获得每个病例的细胞级特征
        fprintf('%d/%d\n', i, numel(indices));
    end
    feas = cell2mat(feas);
    s1 = size(feas, 1);
    feas = feas(randperm(s1, round(s1/10)), :);   %从其中再随机获得 round(s1/10)个细胞级特征
    
    nf = size(feas, 2);    %细胞级特征数量
    dic = cell(nf, 1);
    for i = 1:nf
        t1 = tic;
        [~, dic{i}] = kmeans(feas(:, i), dSize, 'maxIter', 500, 'replicates', 5);
        dic{i} = sort(dic{i});
        fprintf('clustering %d/%d finisehd, time %f\n', i, nf, toc(t1));
    end
    save(['dic', num2str(dSize), 'jiankang.mat'], 'dic');
end


%% generate BoW features for each patients  为每个患者生成图像级特征
tab = struct2table(list);   %得到每个病例的细胞级特征
disp(tab.name);
s1 = size(tab, 1);
pids = cell(s1, 1);
for i = 1:s1
    mua=tab.name{i}(1:15);     %一个的时候不用加{} 为tab.name(1:15) 之后改回来  pids{i} = tab.name{i}(1:15);
    pids{i} = mua;
end

[upid, ip, iu] = unique(pids); %[a,b,c]=unique() a相当于set()  b为第一次出现的位置 c为对应位置的原数字相应的排名 越大排名越大
nupid = numel(upid);           %set后的病例数量
feasBow = zeros(nupid, dSize*10);
feasStats = zeros(nupid, 50);
for i = 1:nupid
    inds = find(iu==i);       %找到所有排名低的病例
    feas  = [];
    for j = 1:numel(inds)
        strc = load([dir1, tab.name{inds(j)}]);     %与上一个问题一样 tab.name 改回来 tab.name{inds(j)}
        feas = [feas; cell2mat(strc.cellFeas)];   %获得相同病例的所有细胞级特征
    end
    
    ent  = zeros(1, 10);
    for j = 1:size(feas, 2)   %10个细胞级特征
        idx = knnsearch(dic{j}, feas(:, j));  %knnsearch(a,b) 找到每一个b距离a最短的点的位置  idx为a的行坐标
        counts = hist(idx, 1:dSize);          %十箱直方图
        counts = counts/sum(counts);          %得到频率
        feasBow(i, (1:dSize)+(j-1)*10) = counts;  %储存
        
        a = counts.*log2(counts);
        a(isnan(a))=0;
        ent(j) = -sum(a);        %熵 
    end
    
    feasStats(i, :) = [mean(feas), std(feas), skewness(feas),...
        kurtosis(feas), ent];
    fprintf('patient no. %d/%d finished\n', i, nupid);
end

pid = upid;
imFeas = [feasBow, feasStats];
imData = table(pid, imFeas);
save('imData_jiankang', 'imData');

