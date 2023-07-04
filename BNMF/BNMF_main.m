%% load the real input data
getdata;
%% SVD iniialization
X = [X1 X2 X3 X4];
K = 42;
[U,S,V]=svd(X);
s_sum = sum(S);   %所有列相加
W(:,1) = sqrt(S(1,1))*U(:,1);
H(1,:) = sqrt(S(1,1))*V(:,1)';
for j = 2 : K
    x = U(:,j);
    y = V(:,j);
    xp = (x>=0).*x;
    xn = (x<0).*(-x);
    yp = (y>=0).*y;
    yn = (y<0).*(-y);
    xpnrm = norm(xp);
    ypnrm = norm(yp);
    mp = xpnrm * ypnrm;
    xnnrm = norm(xn);
    ynnrm = norm(yn);
    mn = xnnrm * ynnrm;
    if mp > mn
        u = xp/xpnrm;
        v = yp/ypnrm;
        sigma = mp;
    else
        u = xn/xnnrm;
        v = yn/ynnrm;
        sigma = mn;
    end
    W(:,j) = sqrt(S(j,j)*sigma)*u;
    H(j,:) = sqrt(S(j,j)*sigma)*v';
end
W = abs(W);
save W_original.mat W;
H = abs(H);
H1 = H(:,1:150);
H2 = H(:,151:3497);
H3 = H(:,3498:3665);
H4=H(:,3666:5947);
save H1_original.mat H1;
save H2_original.mat H2;
save H3_original.mat H3;
save H4_original.mat H4;
%% 计算
tic  %计时开始
L1=std2(X1);
L2=std2(X2);
L3=std2(X3);
L4=std2(X4);
[W,H1,H2,H3,H4] = BNMF_comodule(X1,X2,X3,X4,L1,L2,L3,L4,K);
toc  %计时结束










