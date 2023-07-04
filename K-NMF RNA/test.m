% b=[0.5,0.2,0.1;
%    0.1,0.2,0.3];
% 
% c=[0.2,1.2,-3.1;
%    -0.3,0.4,0.2];
% 
% k12=kij(b,c,1);
% for i=1:2
%     disp(k12)
% end
% getdata;
% o=100;
% tic  %计时开始
% k11=kij(X1,X1,o);
% k12=kij(X1,X2,o);
% k13=kij(X1,X3,o);
% k21=k12';
% k22=kij(X2,X2,o);
% k23=kij(X2,X3,o);
% k31=k13';
% k32=k23';
% k33=kij(X3,X3,o);
% save k11_2_kernal.mat k11;
% save k12_2_kernal.mat k12;
% save k13_2_kernal.mat k13;
% save k21_2_kernal.mat k21;
% save k22_2_kernal.mat k22;
% save k23_2_kernal.mat k23;
% save k31_2_kernal.mat k31;
% save k32_2_kernal.mat k32;
% save k33_2_kernal.mat k33;
% toc  %计时结束
% eucl_dist1 = trace((k11-k11*Q1*H1-H1'*Q1'*k11+H1'*Q1'*k11*Q1*H1));
% eucl_dist2 = trace((k22-k22*Q2*H2-H2'*Q2'*k22+H2'*Q2'*k22*Q2*H2));
% eucl_dist3 = trace((k33-k33*Q3*H3-H3'*Q3'*k33+H3'*Q3'*k33*Q3*H3));
% diff1 = eucl_dist1 + eucl_dist2+ eucl_dist3;
% disp(diff1);
% yue=k22*Q2*H2;
% yue2=H2'*Q2'*k22;
% yue3=H2'*Q2'*k22*Q2*H2;
% yue4=k22-yue-yue2+yue3;
% yue5=trace(yue4);
% disp(yue5);
% disp(max(max(Q1)))
% disp(max(max(Q2)))
% disp(max(max(Q3)))
% disp(min(min(Q1)))
% disp(min(min(Q2)))
% disp(min(min(Q3)))
% disp(log2(max(max(Q3))+1))
% MH1 = mean(H1,2);   MH2 = mean(H2,2);  MH3 = mean(H3,2);
% VH1 = std(H1,0,2);  VH2 = std(H2,0,2); VH3 = std(H3,0,2);%为0时，为n-1  1时为n
% save tt/Comodule3_tt_2.mat Co_module;

% c=[1,2,3,4,5];
% d=[-1,3,4,5,6];
% e=corr2(c,d);
% f=corrcoef(c,d);
newX1=k11*Q1*H1;
newX2=k22*Q2*H2;
newX3=k33*Q3*H3;
sumcorr2=0;
disp(corrcoef(k33(6,:),newX3(6,:)))
% for i2=1:2282
%     sumcorr2=sumcorr2+corrcoef(k33(i2,:),newX3(i2,:));
% end
% sumcorr2=sumcorr2/2282;
% disp(sumcorr2);
% sumcorr1=corr2(k11,newX1);
% sumcorr2=corr2(k22,newX2);
% sumcorr3=corr2(k33,newX3);
