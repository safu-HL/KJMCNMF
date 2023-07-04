% Whole solid image
X1 = csvread('WSI_150_411.csv');
X1=abs(X1');
X1=log2(X1+1);
% mRNA
X2 = csvread('mRNA_3347_411.csv');
X2=X2';
X2=log2(X2+1);
% miRNA
X3 = csvread('miRNA_168_411.csv');
X3=X3';
X3=log2(X3+1);
%lncRNA
X4=csvread('lncRNA_2282_411.csv');
X4=X4';
X4=log2(X4+1);
