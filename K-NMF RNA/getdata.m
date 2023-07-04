% mRNA
X1 = csvread('mRNA_3347_411_2.csv');
X1=X1';
% miRNA
X2 = csvread('miRNA_168_411_2.csv');
X2=X2';
%lncRNA
X3=csvread('lncRNA_2282_411_2.csv');
X3=X3';
%Correlation miRNA&mRNA miRNA&lncRNA
A = csvread('mRNA_miRNA_3347_168.csv');
B = csvread('lncRNA_miRNA_2282_168.csv');