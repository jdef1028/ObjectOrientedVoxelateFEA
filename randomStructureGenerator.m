function randomStructureGenerator(size, VF)
Bimg_coarse = zeros(size, size, size);
currentVF = 0;
while currentVF < VF
    xs = randi(size);
    ys = randi(size);
    zs = randi(size);
    Bimg_coarse(xs,ys,zs) = 1;
    currentVF = sum(Bimg_coarse(:))/size^3;
end
save('fine.mat','Bimg_coarse')