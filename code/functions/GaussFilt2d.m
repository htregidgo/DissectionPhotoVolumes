% Y=GaussFilt2d(X,sigma,pixdim)
% Gaussian blurring. Pixel dimension is considered (1mm isotropic by default)
% pixdim and sigma can be a vector of length 1 (if isotropic) or 2 (anisotropic)
function Y=GaussFilt2d(X,sigma,pixdim)
if exist('pixdim','var')==0
    pixdim=[1 1];
end
if numel(pixdim)==1
    pixdim=pixdim*ones(1,2);
end
if numel(sigma)==1
    sigma=sigma*ones(1,2);
end
sigma=sigma./pixdim; % in voxels
w=ceil(2.5*sigma);

Y=X;
if sigma(1)>0, v=-w(1):w(1); g=exp(-0.5*v.*v/(sigma(1)*sigma(1))); g=g/sum(g); Y=imfilter(Y,reshape(g,[length(v) 1]),'replicate'); end
if sigma(2)>0, v=-w(2):w(2); g=exp(-0.5*v.*v/(sigma(2)*sigma(2))); g=g/sum(g); Y=imfilter(Y,reshape(g,[1 length(v)]),'replicate'); end




