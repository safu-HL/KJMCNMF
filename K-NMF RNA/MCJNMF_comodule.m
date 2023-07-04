function [Q1,Q2,Q3,H1,H2,H3] = MCJNMF_comodule(X1,X2,X3,A,B,a1,r1,r2,L1,L2,K,w,myfid)
%
% INPUT:
% X1 (N,M1): input profile matrix
% X2 (N,M2): input profile matrix
% X3 (N,M3): input profile matrix
% A (M1,M2): input adjacent matrix
% B (M1,M3): input adjacent matrix
% r1       : limit the growth of W
% r2       : limit the growth of H1, H2, H3
% L1       : weigh the must link constraints in A
% L2       : weigh the must link constraints in B
% K        : Number of components

% avoid this kind of column or row: sum == 0
index = find(sum(X1,1) == 0);
X1(:,index) = X1(:,index) + eps;

index = find(sum(X2,1) == 0);
X2(:,index) = X2(:,index) + eps;

index = find(sum(X3,1) == 0);
X3(:,index) = X3(:,index) + eps;

index = find(sum(A,1) == 0);
A(:,index) = A(:,index) + eps;

index = find(sum(B,1) == 0);
B(:,index) = B(:,index) + eps;  


% set the iteration number and initiate the output matrices
nloop = 1; 
verbose=1;
load('k11_2_kernal.mat');
load('k22_2_kernal.mat');
load('k33_2_kernal.mat');
[n,m1] = size(X1);
[n,m2] = size(X2);
[n,m3] = size(X3);

bestQ1=zeros(m1,K);
bestQ2=zeros(m2,K);
bestQ3=zeros(m3,K);
bestH1=zeros(K,m1);
bestH2=zeros(K,m2);
bestH3=zeros(K,m3);

bestobj1=1000000000;
bestobj2=1000000000;
bestobj3=1000000000;
fid = fopen([ 'guocheng-K/Record_K=' int2str(K) '_L1=' num2str(L1) '_L2=' num2str(L2) '_a1=' num2str(a1) 'r1=' num2str(r1) '_r2=' num2str(r2) '_w=' num2str(w) '.txt'],'wt+');
for iloop=1:nloop
    if verbose 
        fprintf(fid,' iteration %d\n',iloop); 
    end 
    
    maxiter=500; 
    speak=1; 
    [Q1,Q2,Q3,H1,H2,H3]=MCJNMF(X1, X2, X3, A, B, a1, r1, r2, L1, L2, K, maxiter, speak, fid,w,myfid);
    % compute residue
    newobj1 = trace((k11-k11*Q1*H1-H1'*Q1'*k11+H1'*Q1'*k11*Q1*H1));
    newobj2 =trace((k22-k22*Q2*H2-H2'*Q2'*k22+H2'*Q2'*k22*Q2*H2));
    newobj3 = trace((k33-k33*Q3*H3-H3'*Q3'*k33+H3'*Q3'*k33*Q3*H3));
    
    if (newobj1<bestobj1)||(newobj2<bestobj2)||(newobj3<bestobj3)
        bestobj1 = newobj1;
        bestobj2 = newobj2;
        bestobj3 = newobj3;
        bestQ1 = Q1;
        bestQ2 = Q2;
        bestQ3 = Q3;
        bestH1 = H1;
        bestH2 = H2;
        bestH3 = H3;
    end
end
fclose(fid);
%%  compute the modules according to bestW, bestH1 bestH2 and bestH3
Q1 = bestQ1;Q2 = bestQ2;Q3 = bestQ3;H1 = bestH1; H2 = bestH2; H3=bestH3;
% save major_Q1_best.mat Q1;
% save major_Q2_best.mat Q2;
% save major_Q3_best.mat Q3;
% save major_H1_best.mat H1;
% save major_H2_best.mat H2;
% save major_H3_best.mat H3;
