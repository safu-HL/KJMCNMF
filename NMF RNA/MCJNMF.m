function [W,H1,H2,H3,H4] = MCJNMF(X1, X2, X3, X4, A, B, C, a1, r1, r2, L1, L2,L3, K, maxiter, speak, fid)
%
% Multiple NMF using euclidean distance update equations:
%
% Lee, D..D., and Seung, H.S., (2001), 'Algorithms for Non-negative Matrix
% Factorization', Adv. Neural Info. Proc. Syst. 13, 556-562.
%
% INPUT:
% X1 (N,M1): N (dimensionallity) x M1 (feature 1) non negative input matrix
% X2 (N,M2): N (dimensionallity) x M2 (feature 2) non negative input matrix
%%%  X3(N,M3): N (dimensionallity) x M3 (feature 3) non negative input matrix
% A (M2,M2): M2 x M2, non negative input matrix
% B (M1,M2): M1 x M2, non negative input matrix 
% C (M3,M2):M3 x M2, non negative input matrix 

% r1       : limit the growth of W
% r2       : limit the growth of H1 H2 and H3
% L1       : weigh the must link constraints in A
% L2       : weigh the must link constraints in B
% K        : Number of components
% maxiter  : Maximum number of iterations to run
% speak    : prints iteration count and changes in connectivity matrix 
%            elements unless speak is 0
%fid       : file identifier which is used to store the changes record

% OUTPUT:
% W        : N x K matrix
% H1       : K x M1 matrix
% H2       : K x M2 matrix
% H3       : K x M3 matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User adjustable parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print_iter = 20; % iterations between print on screen and convergence test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for negative values in X1 and X2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (min(min(X1)) < 0) || (min(min(X2)) < 0) || (min(min(X3)) < 0) || (min(min(X4)) < 0)
    error('Input matrix elements can not be negative');     %不能为负
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for same rows in X1 X2 and X3 x4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n1,m1] = size(X1);
[n2,m2] = size(X2);
[n3,m3] = size(X3);
[n4,m4] = size(X4);

if (n1 ~= n2) || (n1 ~= n3) || (n1~=n4)|| (n2 ~= n3) || (n2~=n4)|| (n3~=n4)
    error('Input matrices should have the same rows');
    return
end
n = n1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize random W, H1 H2 H3 H4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% W = rand(n,K);
% H1 = rand(K,m1);
% H2 = rand(K,m2);
% H3 = rand(K,m3);
% H4 = rand(K,m4);
load('W_original.mat');
load('H1_original.mat');
load('H2_original.mat');
load('H3_original.mat');
load('H4_original.mat');

% use W*H to test for convergence
Xr_old1 = W*H1;
Xr_old2 = W*H2;
Xr_old3 = W*H3;
Xr_old4 = W*H4;

E1 = ones(K,m1);
E2 = ones(K,m2);
E3 = ones(K,m3);
E4 = ones(K,m4);

for iter = 1 : maxiter
    % Euclidean multiplicative method
    W = W.*([H1 H2 H3 H4]*[X1 X2 X3 X4]')'./(W*([H1 H2 H3 H4]*[H1 H2 H3 H4]'+r1*eye(K))+eps);    
    HH1 = H1.*(W'*X1 + 2*a1*H1 + L1/2*H2*C')./((W'*W + 2*a1*(H1*H1'))*H1 + r2/2*E1+eps);     
    HH2 = H2.*(W'*X2 + 2*a1*H2 + L1/2*H1*C + L2/2*H3*A)./((W'*W + 2*a1*(H2*H2'))*H2 + r2/2*E2+eps);
    HH3 = H3.*(W'*X3 + 2*a1*H3 + L2/2*H2*A' + L3/2*H4*B')./((W'*W + 2*a1*(H3*H3'))*H3 + r2/2*E3+eps);
    HH4 = H4.*(W'*X4 + 2*a1*H4 + L3/2*H3*B)./((W'*W + 2*a1*(H4*H4'))*H4 + r2/2*E4+eps);
 
    H1 = HH1;
    H2 = HH2;
    H3 = HH3;
    H4 = HH4;

   % iter
    % print to screen
    if (rem(iter,print_iter) == 0) && speak 
        Xr1 = W*H1+eps;            
        Xr2 = W*H2+eps;
        Xr3 = W*H3+eps;
        Xr4 = W*H4+eps;
        diff_step = sum(sum(abs(Xr_old1-Xr1)))+sum(sum(abs(Xr_old2-Xr2)))+sum(sum(abs(Xr_old3-Xr3)))+sum(sum(abs(Xr_old4-Xr4)));
        
 	    diff2 = -L1*trace(H1*C*H2');
        diff3 = -L2*trace(H3*A*H2');
        diff4 = -L3*trace(H3*B*H4');
        diff5 = r1*sum(sum(W.*W));     
        diff6 = r2*(sum(sum(H1).^2)+sum(sum(H2).^2)+sum(sum(H3).^2)+sum(sum(H4).^2));
        diff7 = a1*(trace(H1'*H1*H1'*H1)+trace(H2'*H2*H2'*H2)+trace(H3'*H3*H3'*H3)+trace(H4'*H4*H4'*H4));
        diff8 = r2*(sum(sum(H1).^2)+sum(sum(H2).^2)+sum(sum(H3).^2)+sum(sum(H4).^2))+a1*(trace(H1'*H1*H1'*H1)+trace(H2'*H2*H2'*H2)+trace(H3'*H3*H3'*H3)+trace(H4'*H4*H4'*H4));    %sum求列和,两个sum相当于求所有元素之和
	    Xr_old1 = Xr1;
        Xr_old2 = Xr2;
        Xr_old3 = Xr3;
        Xr_old4 = Xr4;
        eucl_dist1 = nmf_euclidean_dist(X1,W*H1);
        eucl_dist2 = nmf_euclidean_dist(X2,W*H2);
        eucl_dist3 = nmf_euclidean_dist(X3,W*H3);
        eucl_dist4 = nmf_euclidean_dist(X4,W*H4);
        diff1 = eucl_dist1 + eucl_dist2+ eucl_dist3+eucl_dist4;
 		diff = diff1 + diff2 + diff3 + diff4 + diff5 +diff6+diff7;
        errorx1 = mean(mean(abs(X1-W*H1)))/mean(mean(X1));
        errorx2 = mean(mean(abs(X2-W*H2)))/mean(mean(X2));
        errorx3 = mean(mean(abs(X3-W*H3)))/mean(mean(X3));
        errorx4 = mean(mean(abs(X4-W*H4)))/mean(mean(X4));
        errorx = errorx1 + errorx2 + errorx3+errorx4;
          disp(['SVD_Iter = ',int2str(iter),...
            ', relative error = ',num2str(errorx),...
            ', diff_step = ', num2str(diff_step),...
            ', eucl dist ' num2str(diff1),...
            ', NR1 =', num2str(diff2),...
            ', NR2 =', num2str(diff3),...
            ', NR3 =', num2str(diff4),...
            ', SC1 =', num2str(diff5),...
            ', SC2 =', num2str(diff6),...
            ', OC =',  num2str(diff7),...
            ', errorx1=', num2str(errorx1),...
            ', errorx2=', num2str(errorx2),...
            ', errorx3=', num2str(errorx3),...
            ', errorx4=', num2str(errorx4)])
         fprintf(fid,'%s\n',[sprintf('Iter = \t'),int2str(iter),...
            sprintf('\t relative error = \t'),num2str(errorx),...
            sprintf('\t diff_step = \t'),num2str(diff_step),...
            sprintf('\t diff = \t'), num2str(diff),...
            sprintf('\t diff1 = \t'), num2str(diff1),...
            sprintf('\t diff2 = \t'), num2str(diff2),...
            sprintf('\t diff3 = \t'), num2str(diff3),...
            sprintf('\t diff4 = \t'), num2str(diff4),...
		    sprintf('\t diff5 = \t'), num2str(diff5),...
		    sprintf('\t diff6 = \t'), num2str(diff6),...
            sprintf('\t diff7 = \t'), num2str(diff7),...
            sprintf('\t diff8 = \t'),num2str(diff8)]);
%         if iter==maxiter
%             fprintf(myfid,'%s\n',[sprintf('\t relative error= '),num2str(errorx),...
%             sprintf('\t eucl dist= '),num2str(diff1)]);
%         end
        if errorx < 10^(-5), break, end
    end
end

function err = nmf_euclidean_dist(X,Y)
err = sum(sum((X-Y).^2));
