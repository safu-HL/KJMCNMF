function [W,H1,H2,H3,H4] = JNMF_comodule(X1,X2,X3,X4,K)
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

index = find(sum(X4,1) == 0);
X4(:,index) = X4(:,index) + eps;

% set the iteration number and initiate the output matrices
nloop = 1; 
verbose=1;

[n,m1] = size(X1);
[n,m2] = size(X2);
[n,m3] = size(X3);
[n,m4] = size(X4);

bestW=zeros(n,K);
bestH1=zeros(K,m1);
bestH2=zeros(K,m2);
bestH3=zeros(K,m3);
bestH4=zeros(K,m4);

bestobj1=1000000000;
bestobj2=1000000000;
bestobj3=1000000000;
bestobj4=1000000000;
fid = fopen([ 'KRecord_K=' int2str(K)  '.txt'],'wt+');
for iloop=1:nloop
    if verbose 
        fprintf(fid,' iteration %d\n',iloop); 
    end 
    
    maxiter=1000; 
    speak=1; 
    [W,H1,H2,H3,H4]=JNMF(X1, X2, X3, X4,K, maxiter, speak, fid);
    % compute residue
    newobj1 = sum(sum((X1-W*H1).^2));
    newobj2 = sum(sum((X2-W*H2).^2));
    newobj3 = sum(sum((X3-W*H3).^2));
    newobj4 = sum(sum((X4-W*H4).^2));
    
    if (newobj1<bestobj1)||(newobj2<bestobj2)||(newobj3<bestobj3)||(newobj4<bestobj4)        
        bestobj1 = newobj1;
        bestobj2 = newobj2;
        bestobj3 = newobj3;
        bestobj4 = newobj4;
        bestW = W;
        bestH1 = H1;
        bestH2 = H2;
        bestH3 = H3;
        bestH4 = H4;
    end
end
fclose(fid);
%%  compute the modules according to bestW, bestH1 bestH2 and bestH3
W = bestW; H1 = bestH1; H2 = bestH2; H3=bestH3;H4=bestH4;
save KW_best.mat W;
% save KH1_best.mat H1;
save KH2_best.mat H2;
save KH3_best.mat H3;
save KH4_best.mat H4;
