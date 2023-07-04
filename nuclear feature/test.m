jiankang_shuzu=[];
for i=1:186
    disp(imData{i,2});
    jiankang_shuzu=[jiankang_shuzu;imData{i,2}];
end
xlswrite('WSI_150_186',jiankang_shuzu');

     