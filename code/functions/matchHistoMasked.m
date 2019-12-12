% Matches the histogram of an image to either: (a) a reference histogram; or
% (b) the histogram of other image. It ignores pixels/voxels equal to zero.
function Imatched = matchHistoMasked(I, REF)

x1=I(I>0);
h1 = imhist(x1); %// Compute histograms
c1 = cumsum(h1) / numel(x1); %// Compute CDFs

if abs(sum(REF(:))-1)<1e-9 % if it's already a histogram
    c2=cumsum(REF(:));
else % if it's an image, compute histogram (and CDF)
    x2=REF(REF>0);
    h2 = imhist(x2);
    c2 = cumsum(h2) / numel(x2);
end

M = zeros(256,1,'uint8');
for idx = 1 : 256
    [~,ind] = min(abs(c1(idx) - c2));
    M(idx) = ind-1;
end
out = M(double(x1)+1);
Imatched=zeros(size(I),'uint8');
Imatched(I>0)=out;
