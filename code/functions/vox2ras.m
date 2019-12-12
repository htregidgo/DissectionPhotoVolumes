function ras = vox2ras(vox,vox2ras0)

if size(vox,2)==3
    vox2=[vox'-1; ones(1,size(vox,1))]; % matlab 1-based indexing
else
    vox2=[vox-1; ones(1,size(vox,2))];
end

ras = vox2ras0 * vox2([2 1 3 4],:); % matlab flip

ras = ras(1:3,:);




