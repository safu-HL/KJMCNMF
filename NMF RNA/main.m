%% load the real input data
getdata;
%% preprocess the real input data
A = sparse(abs(A')); %miRNA&mRNA   稀疏性矩阵   格式为 （坐标）  数据
B = sparse(abs(B')); %miRNA&lncRNA
C = sparse(abs(C')); %WSI&mRNA
%% SVD iniialization
X = [X1 X2 X3 X4];
K = 42;
[U,S,V]=svd(X);
s_sum = sum(S);   %所有列相加
W(:,1) = sqrt(S(1,1))*U(:,1);
H(1,:) = sqrt(S(1,1))*V(:,1)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% myfid=fopen( 'canshuxuanze4.txt','wt+');
% weici=1;
% for i1=2:70
%     tic  %计时开始
%     K=i1;
%     H=zeros(K,5947);
%     W=zeros(411,K);
%     W(:,1) = sqrt(S(1,1))*U(:,1);
%     H(1,:) = sqrt(S(1,1))*V(:,1)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%% applying MDJNMF
canshu=[0.0001,0.001,0.01,0.1,1,10];
L1 = 0.0001; L2 = 1;L3=0.0001; r1 = 0.0001; r2 =0.0001;  a = 0.0001; K=15;
%% L1,L2,L3,a参数选择
% myfid=fopen( 'canshuxuanze1.txt','wt+');
% weici=0;
% for i1=1:6
%     for i2=1:6
%         for i3=1:6
%             for i4=1:6
%                 L1=canshu(i1);
%                 L2=canshu(i2);
%                 L3=canshu(i3);
%                 a=canshu(i4);
%                 tic  %计时开始
%                 weici=weici+1;
%                 fprintf(myfid,'%s\n',['Record_K=' int2str(K) '_L1=' num2str(L1) '_L2=' num2str(L2) '_L3=' num2str(L3) '_a1=' num2str(a) 'r1=' num2str(r1) '_r2=' num2str(r2)]);
%                 fprintf(myfid,'%s',[sprintf('Iter= '),int2str(weici)]);
%                 MCJNMF_comodule(X1,X2,X3,X4,A,B,C,a,r1,r2,L1,L2,L3,K,myfid);
%                 toc  %计时结束
%             end
%         end
%     end
% end
% fclose(myfid);
%% r1,r2参数选择
% myfid=fopen( 'canshuxuanze2.txt','wt+');
% weici=0;
% for i1=1:6
%     for i2=1:6
%                 r1=canshu(i1);
%                 r2=canshu(i2);
%                 tic  %计时开始
%                 weici=weici+1;
%                 fprintf(myfid,'%s\n',['Record_K=' int2str(K) '_L1=' num2str(L1) '_L2=' num2str(L2) '_L3=' num2str(L3) '_a1=' num2str(a) 'r1=' num2str(r1) '_r2=' num2str(r2)]);
%                 fprintf(myfid,'%s',[sprintf('Iter= '),int2str(weici)]);
%                 MCJNMF_comodule(X1,X2,X3,X4,A,B,C,a,r1,r2,L1,L2,L3,K,myfid);
%                 toc  %计时结束
%     end
% end
% fclose(myfid);
%% K值选择
%                 
%                 weici=weici+1;
%                 fprintf(myfid,'%s\n',['Record_K=' int2str(K) '_L1=' num2str(L1) '_L2=' num2str(L2) '_L3=' num2str(L3) '_a1=' num2str(a) 'r1=' num2str(r1) '_r2=' num2str(r2)]);
%                 fprintf(myfid,'%s',[sprintf('Iter= '),int2str(weici)]);
%                 [W,H1,H2,H3,H4] =MCJNMF_comodule(X1,X2,X3,X4,A,B,C,a,r1,r2,L1,L2,L3,K,myfid);
%                 tt = 2;
%                 [Co_module] =Comodule_selection(W, H1, H2, H3, H4, tt);
%                 newX1=W*H1;
%                 newX2=W*H2;
%                 newX3=W*H3;
%                 newX4=W*H4;
%                 sumcorr1=0;
%                 sumcorr2=0;
%                 sumcorr3=0;
%                 sumcorr4=0;
%                 for i2=1:411
%                     sumcorr1=sumcorr1+corrcoef(X1(i2,:),newX1(i2,:));
%                     sumcorr2=sumcorr2+corrcoef(X2(i2,:),newX2(i2,:));
%                     sumcorr3=sumcorr3+corrcoef(X3(i2,:),newX3(i2,:));
%                     sumcorr4=sumcorr4+corrcoef(X4(i2,:),newX4(i2,:));
%                 end
%                 sumcorr1=sumcorr1/411;
%                 sumcorr2=sumcorr2/411;
%                 sumcorr3=sumcorr3/411;
%                 sumcorr4=sumcorr4/411;
%                 minH1=5000;
%                 minH2=5000;
%                 minH3=5000;
%                 minH4=5000;
%                 CB1=zeros(1,150);
%                 CB2=zeros(1,3347);
%                 CB3=zeros(1,168);
%                 CB4=zeros(1,2282);
%                 for i2=1:K
%                     %H1
%                     CH1=cell2mat(Co_module(i2,1));
%                     LH=length(CH1);
%                     if LH<minH1
%                         minH1=length(CH1);
%                     end
%                     for i3=1:LH
%                         CB1(1,CH1(1,i3))=CB1(1,CH1(1,i3))+1;
%                     end
%                     
%                     %H2
%                     CH2=cell2mat(Co_module(i2,2));
%                     LH=length(CH2);
%                     if LH<minH2
%                         minH2=length(CH2);
%                     end
%                     for i3=1:LH
%                         CB2(1,CH2(1,i3))=CB2(1,CH2(1,i3))+1;
%                     end
%                     
%                     %H3
%                     CH3=cell2mat(Co_module(i2,3));
%                     LH=length(CH3);
%                     if LH<minH3
%                         minH3=length(CH3);
%                     end
%                     for i3=1:LH
%                         CB3(1,CH3(1,i3))=CB3(1,CH3(1,i3))+1;
%                     end
%                     
%                     %H4
%                     CH4=cell2mat(Co_module(i2,4));
%                     LH=length(CH4);
%                     if LH<minH4
%                         minH4=length(CH4);
%                     end
%                     for i3=1:LH
%                         CB4(1,CH4(1,i3))=CB4(1,CH4(1,i3))+1;
%                     end
%                 end
%                 if (K/4)>=1
%                     yuzhi=K/4;
%                 else
%                     yuzhi=1;
%                 end
%                 CC1=(length(find(CB1(1,:)>yuzhi)))/length(find(CB1(1,:)>0));
%                 CC2=(length(find(CB2(1,:)>yuzhi)))/length(find(CB2(1,:)>0));
%                 CC3=(length(find(CB3(1,:)>yuzhi)))/length(find(CB3(1,:)>0));
%                 CC4=(length(find(CB4(1,:)>yuzhi)))/length(find(CB4(1,:)>0));
%                 fprintf(myfid,'%s\n\n',[sprintf('corr1= '),num2str(sumcorr1(1,2)),...
%                     sprintf('\t corr2= '),num2str(sumcorr2(1,2)),...
%                     sprintf('\t corr3= '),num2str(sumcorr3(1,2)),...
%                     sprintf('\t corr4= '),num2str(sumcorr4(1,2)),...
%                     sprintf('\t chong1= '),num2str(CC1),...
%                     sprintf('\t chong2= '),num2str(CC2),...
%                     sprintf('\t chong3= '),num2str(CC3),...
%                     sprintf('\t chong4= '),num2str(CC4)]);
%                 toc  %计时结束
% end
% fclose(myfid);
% for i1=1:21
tic  %计时开始
[W,H1,H2,H3,H4] = MCJNMF_comodule(X1,X2,X3,X4,A,B,C,a,r1,r2,L1,L2,L3,K);
toc  %计时结束
 %% get comodules(index)
% tt = 2;
% [Co_module] =Comodule_selection(W, H1, H2, H3, H4, tt);
% save tt/Comodule_15_tt_2.mat Co_module;
% % filename_chong=['tt/Comodule_tt_2_',num2str(i1),'.mat'];
% % save(filename_chong,'Co_module');
% % 
% % end
% 
% %% module elements extraction
% B1=zeros(15,1500);
% B2=zeros(15,1500);
% B3=zeros(15,1500);
% B4=zeros(15,1500);
% B1(B1==0)=[];
% B2(B2==0)=[];
% B3(B3==0)=[];
% B4(B4==0)=[];
% for i=1:15
%     A1=Co_module(i,1);
%     A2=Co_module(i,2);
%     A3=Co_module(i,3);
%     A4=Co_module(i,4);
%     A11=cell2mat(A1);
%     A22=cell2mat(A2);
%     A33=cell2mat(A3);
%     A44=cell2mat(A4);
%     k1=length(A11);
%     k2=length(A22);
%     k3=length(A33);
%     k4=length(A44);
%     for j1=1:k1
%         B1(i,j1)=A11(1,j1);
%     end
%     for j2=1:k2
%         B2(i,j2)=A22(1,j2);
%     end
%     for j3=1:k3
%         B3(i,j3)=A33(1,j3);
%     end
%     for j4=1:k4
%         B4(i,j4)=A44(1,j4);
%     end
% end
% B1(B1==0)=NaN;
% B2(B2==0)=NaN;
% B3(B3==0)=NaN;
% B4(B4==0)=NaN;
% 
% xlswrite('tt/Co_module_tt_2_K_15_svd_WSI_2.xlsx',B1');
% xlswrite('tt/Co_module_tt_2_K_15_svd_mRNA_2.xlsx',B2');
% xlswrite('tt/Co_module_tt_2_K_15_svd_miRNA_2.xlsx',B3');
% xlswrite('tt/Co_module_tt_2_K_15_svd_lncRNA_2.xlsx',B4');