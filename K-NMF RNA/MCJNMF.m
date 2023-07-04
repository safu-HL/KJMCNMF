function [Q1,Q2,Q3,H1,H2,H3] = MCJNMF(X1, X2, X3, A, B, a1, r1, r2, L1, L2, K, maxiter, speak, fid,w,myfid)
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

% if (min(min(X1)) < 0) || (min(min(X2)) < 0) || (min(min(X3)) < 0) || (min(min(X4)) < 0)
%     error('Input matrix elements can not be negative');     %不能为负
%     return
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for same rows in X1 X2 and X3 x4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n1,m1] = size(X1);
[n2,m2] = size(X2);
[n3,m3] = size(X3);

if (n1 ~= n2) || (n1 ~= n3) || (n2 ~= n3)
    error('Input matrices should have the same rows');
    return
end
n = n1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize random W, H1 H2 H3 H4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%W = rand(n,K);
Q1=rand(m1,K);
Q2=rand(m2,K);
Q3=rand(m3,K);
H1 = rand(K,m1);
H2 = rand(K,m2);
H3 = rand(K,m3);
% load('Q1_original.mat');
% load('Q2_original.mat');
% load('Q3_original.mat');
% load('H1_original.mat');
% load('H2_original.mat');
% load('H3_original.mat');
load('k11_2_kernal.mat');
load('k12_2_kernal.mat');
load('k13_2_kernal.mat');
load('k21_2_kernal.mat');
load('k22_2_kernal.mat');
load('k23_2_kernal.mat');
load('k31_2_kernal.mat');
load('k32_2_kernal.mat');
load('k33_2_kernal.mat');
% use W*H to test for convergence


Xr_old1 = Q1*H1+eps;
Xr_old2 = Q2*H2+eps;
Xr_old3 = Q3*H3+eps;

E1 = ones(K,m1);
E2 = ones(K,m2);
E3 = ones(K,m3);

for iter = 1 : maxiter
    % Euclidean multiplicative method
    QQ1=Q1.*(k11*H1'+w*k12*Q2+w*k13*Q3)./(k11*Q1*H1*H1'+r1*Q1+2*w*k11*Q1+eps);
    QQ2=Q2.*(k22*H2'+w*k21*Q1+w*k23*Q3)./(k22*Q2*H2*H2'+r1*Q2+2*w*k22*Q2+eps);
    QQ3=Q3.*(k33*H3'+w*k31*Q1+w*k32*Q2)./(k33*Q3*H3*H3'+r1*Q3+2*w*k33*Q3+eps);
    HH1 = H1.*(Q1'*k11 + 2*a1*H1 + L1/2*H2*A)./(Q1'*k11*Q1*H1+2*a1*H1*H1'*H1+r2/2*E1+eps);     
    HH2 = H2.*(Q2'*k22 + 2*a1*H2 + L1/2*H1*A'+L2/2*H3*B')./(Q2'*k22*Q2*H2+2*a1*H2*H2'*H2+r2/2*E2+eps);     
    HH3 = H3.*(Q3'*k33 + 2*a1*H3 + L2/2*H2*B)./(Q3'*k33*Q3*H3+2*a1*H3*H3'*H3+r2/2*E3+eps);     

    Q1 = QQ1;
    Q2 = QQ2;
    Q3 = QQ3;
%     H1 = log2(HH1+1);
%     H2 = log2(HH2+1);
%     H3 = log2(HH3+1);
    H1 = HH1;
    H2 = HH2;
    H3 = HH3;
   % iter
    % print to screen
    if (rem(iter,print_iter) == 0) && speak 
        Xr1 = Q1*H1+eps;            
        Xr2 = Q2*H2+eps;
        Xr3 = Q3*H3+eps;
        diff_step = sum(sum(abs(Xr_old1-Xr1)))+sum(sum(abs(Xr_old2-Xr2)))+sum(sum(abs(Xr_old3-Xr3)));
        
 	    diff2 = -L1*trace(H2*A*H1');
        diff3 = -L2*trace(H2*B*H3');
        diff4 = r1*(sum(sum(Q1.*Q1))+sum(sum(Q2.*Q2))+sum(sum(Q3.*Q3)));     
        diff5 = r2*(sum(sum(H1).^2)+sum(sum(H2).^2)+sum(sum(H3).^2));
        diff6 = a1*(trace(H1'*H1*H1'*H1)+trace(H2'*H2*H2'*H2)+trace(H3'*H3*H3'*H3));
        diff7 = r2*(sum(sum(H1).^2)+sum(sum(H2).^2)+sum(sum(H3).^2))+a1*(trace(H1'*H1*H1'*H1)+trace(H2'*H2*H2'*H2)+trace(H3'*H3*H3'*H3));    %sum求列和,两个sum相当于求所有元素之和
	    Xr_old1 = Xr1;
        Xr_old2 = Xr2;
        Xr_old3 = Xr3;
        eucl_dist1 = trace((k11-k11*Q1*H1-H1'*Q1'*k11+H1'*Q1'*k11*Q1*H1));
        eucl_dist2 = trace((k22-k22*Q2*H2-H2'*Q2'*k22+H2'*Q2'*k22*Q2*H2));
        eucl_dist3 = trace((k33-k33*Q3*H3-H3'*Q3'*k33+H3'*Q3'*k33*Q3*H3));
        diff1 = eucl_dist1 + eucl_dist2+ eucl_dist3;
 		diff = diff1 + diff2 + diff3 + diff4 + diff5 +diff6;
        errorx1 = mean(mean(abs(k11-k11*Q1*H1)))/mean(mean(k11));
        errorx2 = mean(mean(abs(k22-k22*Q2*H2)))/mean(mean(k22));
        errorx3 = mean(mean(abs(k33-k33*Q3*H3)))/mean(mean(k33));
        errorx = errorx1 + errorx2 + errorx3;
          disp(['SVD_Iter = ',int2str(iter),...
            ', relative error = ',num2str(errorx),...
            ', diff_step = ', num2str(diff_step),...
            ', eucl dist ' num2str(diff1),...
            ', errorx1 ' num2str(errorx1),...
            ', errorx2 ' num2str(errorx2),...
            ', errorx3 ' num2str(errorx3),...
            ', NR1 =', num2str(diff2),...
            ', NR2 =', num2str(diff3),...
            ', SC1 =', num2str(diff4),...
            ', SC2 =', num2str(diff5),...
            ', OC =',  num2str(diff6)])
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
            sprintf('\t diff7 = \t'), num2str(diff7)]);
        if iter==maxiter
            fprintf(myfid,'%s\n',[sprintf('\t relative error= '),num2str(errorx),...
            sprintf('\t eucl dist= '),num2str(diff1)]);
        end
        if errorx < 10^(-8), break, end
    end
end

% function err = nmf_euclidean_dist(X,Y)
% err = sum(sum((X-Y).^2));
