function [ fX ] = kij( X,Y,o)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    [~,m1] = size(X);
    [n,m2] = size(Y);
    fX=zeros(m1,m2);
    X=X';
    max=0;
    min=10000000;
    for i=1:m1
        for j=1:m2
            sum=0;
            for k=1:n
                sum=sum+(X(i,k)-Y(k,j))^2;
            end
            if sum>max
                max=sum;
            end
            if sum<min
                min=sum;
            end
            fX(i,j)=exp(-(sum)/(2*o));
        end
    end
    disp('-----------------------------')
    disp(m1)
    disp(max)
    disp(min)
    disp('-----------------------------')
end


