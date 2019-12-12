function vox = ras2vox(ras,vox2ras0)

if size(ras,2)==3
    ras2=[ras'; ones(1,size(ras,1))];
else
    ras2=[ras; ones(1,size(ras,2))];
end

vox = inv(vox2ras0) * ras2;

vox=vox([2 1 3],:)-1;  % matlab flip and 1-based indexing


