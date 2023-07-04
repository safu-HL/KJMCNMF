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
dSize = 10; %����Ϊ10����

%% clustering each type of feature ����ÿ�����͵�Ҫ��
% If dictionary exists, load it
if exist(['dic', num2str(dSize), 'jiankang.mat'], 'file')    %���֮ǰ���������ͼ��ľ������ȡ
    strc = load(['dic', num2str(dSize), 'jiankang.mat']);
    dic = strc.dic;
else % else do clustering
    n = 136;     %N=300  ֮ǰ һ��Ҫ�Ļ���������
    rng('default');    %���������������
    indices = randperm(numel(list), n);  %���ش�1��numel(list)�����n��������n��������ͬ

    feas = cell(n, 1);
    for i = 1:numel(indices)
        strc = load([dir1, list(indices(i)).name]);
        feas{i} = cell2mat(strc.cellFeas);       %������ÿ��������ϸ��������
        fprintf('%d/%d\n', i, numel(indices));
    end
    feas = cell2mat(feas);
    s1 = size(feas, 1);
    feas = feas(randperm(s1, round(s1/10)), :);   %�������������� round(s1/10)��ϸ��������
    
    nf = size(feas, 2);    %ϸ������������
    dic = cell(nf, 1);
    for i = 1:nf
        t1 = tic;
        [~, dic{i}] = kmeans(feas(:, i), dSize, 'maxIter', 500, 'replicates', 5);
        dic{i} = sort(dic{i});
        fprintf('clustering %d/%d finisehd, time %f\n', i, nf, toc(t1));
    end
    save(['dic', num2str(dSize), 'jiankang.mat'], 'dic');
end


%% generate BoW features for each patients  Ϊÿ����������ͼ������
tab = struct2table(list);   %�õ�ÿ��������ϸ��������
disp(tab.name);
s1 = size(tab, 1);
pids = cell(s1, 1);
for i = 1:s1
    mua=tab.name{i}(1:15);     %һ����ʱ���ü�{} Ϊtab.name(1:15) ֮��Ļ���  pids{i} = tab.name{i}(1:15);
    pids{i} = mua;
end

[upid, ip, iu] = unique(pids); %[a,b,c]=unique() a�൱��set()  bΪ��һ�γ��ֵ�λ�� cΪ��Ӧλ�õ�ԭ������Ӧ������ Խ������Խ��
nupid = numel(upid);           %set��Ĳ�������
feasBow = zeros(nupid, dSize*10);
feasStats = zeros(nupid, 50);
for i = 1:nupid
    inds = find(iu==i);       %�ҵ����������͵Ĳ���
    feas  = [];
    for j = 1:numel(inds)
        strc = load([dir1, tab.name{inds(j)}]);     %����һ������һ�� tab.name �Ļ��� tab.name{inds(j)}
        feas = [feas; cell2mat(strc.cellFeas)];   %�����ͬ����������ϸ��������
    end
    
    ent  = zeros(1, 10);
    for j = 1:size(feas, 2)   %10��ϸ��������
        idx = knnsearch(dic{j}, feas(:, j));  %knnsearch(a,b) �ҵ�ÿһ��b����a��̵ĵ��λ��  idxΪa��������
        counts = hist(idx, 1:dSize);          %ʮ��ֱ��ͼ
        counts = counts/sum(counts);          %�õ�Ƶ��
        feasBow(i, (1:dSize)+(j-1)*10) = counts;  %����
        
        a = counts.*log2(counts);
        a(isnan(a))=0;
        ent(j) = -sum(a);        %�� 
    end
    
    feasStats(i, :) = [mean(feas), std(feas), skewness(feas),...
        kurtosis(feas), ent];
    fprintf('patient no. %d/%d finished\n', i, nupid);
end

pid = upid;
imFeas = [feasBow, feasStats];
imData = table(pid, imFeas);
save('imData_jiankang', 'imData');

