SD=[3,4,5,1,1,8];
c1 = sum(SD(SD(1,:) > 1));
for i=1:3
    filename=['testsave_',num2str(i),'.mat'];
    save(filename,'SD');
end

disp(c1)