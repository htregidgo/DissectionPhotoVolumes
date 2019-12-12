function [cost,grad, warpedPhotos, warpedMasks, REFvox2ras0New] = ...
    costFun(params,cogREF,REFmri,Imri,Mmri,IIph,JJph,KKph,...
    REL_NCC_INTRA_WEIGHT,REL_DICE_INTRA_WEIGHT,REL_DICE_INTER_WEIGHT,...
    REL_DETERMINANT_COST, mode, phvals, mvals, refvals, ...
    DiceInterAccumNum, DiceInterAccumDen)

% HARD CODED CONSTANTS
%%%%
% IF YOU CHANGE THIS ONE, THEN ALSO CHANGE IT IN ReconPhotoVolume_joint.m and ReconPhotoVolume_joint_multires.m
FACTOR_AFFINE_MAT = 20; 
FACTOR_SCALING = 20;

Nims = size(Imri.vol,3);

% "Unstack parameters"
if mode == 1 % rigid + similarity
    
    % These are in pixel space, 2D
    theta=params(1:3:end-7)/180*pi;
    tr=params(2:3:end-7);
    tc=params(3:3:end-7);
    
    % For the reference volume
    s=exp(params(end-6)/FACTOR_SCALING);
    rotx=params(end-5)/180*pi;
    roty=params(end-4)/180*pi;
    rotz=params(end-3)/180*pi;
    tx=params(end-2);
    ty=params(end-1);
    tz=params(end);
    
elseif mode == 2  % rigid + affine
    
    % These are in pixel space, 2D
    theta=params(1:3:end-12)/180*pi;
    tr=params(2:3:end-12);
    tc=params(3:3:end-12);
    
    % For the reference volume
    M=reshape(params(end-11:end-3),[3 3])/FACTOR_AFFINE_MAT;
    tx=params(end-2);
    ty=params(end-1);
    tz=params(end);
    
else %  affine + image translation + affine
    
    % These are in pixel space, 2D
    Mph=reshape(params(1:4*Nims),[2 2 Nims])/FACTOR_AFFINE_MAT;
    tr=params(4*Nims+1:2:end-12);
    tc=params(4*Nims+2:2:end-12);
    
    % For the reference volume
    M=reshape(params(end-11:end-3),[3 3])/FACTOR_AFFINE_MAT;
    tx=params(end-2);
    ty=params(end-1);
    tz=params(end);
end

% Resample atlas, unless already provided!
if exist('refvals','var')==0 || isempty(refvals)
    if mode == 1 % similarity
        
        T1 = [1 0 0 -cogREF(1); 0 1 0 -cogREF(2); 0 0 1 -cogREF(3); 0 0 0 1];
        T2=[s 0 0 0;  0 s 0 0; 0 0 s 0;  0 0 0 1];
        T3=[1 0 0 0; 0 cos(rotx) -sin(rotx) 0; 0 sin(rotx) cos(rotx) 0; 0 0 0 1];
        T4=[cos(roty) 0 sin(roty) 0; 0 1 0 0; -sin(roty) 0 cos(roty) 0; 0 0 0 1];
        T5=[cos(rotz) -sin(rotz) 0 0; sin(rotz) cos(rotz) 0 0; 0 0 1 0; 0 0 0 1];
        T6=[1 0 0 cogREF(1); 0 1 0 cogREF(2); 0 0 1 cogREF(3); 0 0 0 1];
        T7=[1 0 0 tx; 0 1 0 ty; 0 0 1 tz; 0 0 0 1];
        T=T7*T6*T5*T4*T3*T2*T1;
        
    else % affine
        T1 = [1 0 0 -cogREF(1); 0 1 0 -cogREF(2); 0 0 1 -cogREF(3); 0 0 0 1];
        T2 = [[M zeros(3,1)]; 0 0 0 1];
        T3=[1 0 0 cogREF(1); 0 1 0 cogREF(2); 0 0 1 cogREF(3); 0 0 0 1];
        T4=[1 0 0 tx; 0 1 0 ty; 0 0 1 tz; 0 0 0 1];
        T=T4*T3*T2*T1;
    end
    
    voxref = (inv(T * REFmri.vox2ras0) * Imri.vox2ras0) * ...
        [JJph(:)'-1;  IIph(:)'-1;  KKph(:)'-1; ones(1,numel(IIph))];
    voxref=voxref([2 1 3],:)+1;
    refvals=interpn(REFmri.vol,voxref(1,:),voxref(2,:),voxref(3,:));
    refvals=reshape(refvals,size(IIph));
    refvals(isnan(refvals))=0;
end

% Now resample photos, unless provided already!
if exist('phvals','var')==0 || isempty(phvals) ...
        || exist('mvals','var')==0 || isempty(mvals)
    phvals = zeros([size(IIph), 3]);
    mvals = zeros(size(IIph));
    imid=size(IIph,1)/2;
    jmid=size(IIph,2)/2;
    I = IIph(:,:,1)-imid;
    J = JJph(:,:,1)-jmid;
    for z=1:size(IIph,3)
        if mode==1 || mode == 2 % rigid
            I2 = cos(theta(z))*I(:) - sin(theta(z)) * J(:) + tr(z) + imid;
            J2 = sin(theta(z))*I(:) + cos(theta(z)) * J(:) + tc(z) + jmid;
        else
            I2 = Mph(1,1,z)*I(:) + Mph(1,2,z) * J(:) + tr(z) + imid;
            J2 = Mph(2,1,z)*I(:) + Mph(2,2,z) * J(:) + tc(z) + jmid;
        end
        for c = 1:3
            vals = interpn(Imri.vol(:,:,z,c),I2(:),J2(:));
            vals = reshape(vals,size(I));
            phvals(:,:,z,c)=vals;
        end
        vals = interpn(Mmri.vol(:,:,z),I2(:),J2(:));
        vals = reshape(vals,size(I));
        mvals(:,:,z)=vals;
    end
    phvals(isnan(phvals))=0;
    mvals(isnan(mvals))=0;
    
end

% Now compute cost
% For the Dice score, we need a little trick, in case that we're only
% computing a mini-volume (see computation of gradient below)  
DiceInterNums=squeeze(sum(sum(mvals.*refvals,1),2));
DiceInterDens=squeeze(sum(sum(mvals.*mvals,1),2)+sum(sum(refvals.*refvals,1),2));
if exist('DiceInterAccumDen','var')==0 || isempty(DiceInterAccumDen)
   DiceInterAccumNum=0;
   DiceInterAccumDen=0;
end
diceInter = 2 * (sum(DiceInterNums) + DiceInterAccumNum) ...
    / (sum(DiceInterDens) + DiceInterAccumDen);

X = phvals/255;
Y = mvals;
rhos = zeros([size(X,3)-1,3]);
diceIntras = zeros([size(X,3)-1,1]);
% el = strel('disk',3);
for z=2:size(X,3)
    m = mvals(:,:,z)>0 | mvals(:,:,z-1)>0;
    % m = imdilate(m, el);
    for c=1:3
        aux = X(:,:,z,c);
        a = aux(m);
        aux=X(:,:,z-1,c);
        b = aux(m);
        rhos(z-1,c)=corr(a,b);
    end
    aux=Y(:,:,z);
    a=aux(m);
    aux=Y(:,:,z-1);
    b=aux(m);
    diceIntras(z-1) = 2 * sum(a.*b) / (sum(a.*a)+sum(b.*b));
end

cost = -  REL_NCC_INTRA_WEIGHT * sum(rhos(:)) ...
       - REL_DICE_INTRA_WEIGHT * sum(diceIntras) ...
       - REL_DICE_INTER_WEIGHT * diceInter;

if mode==3
    dets=zeros(1,size(Mph,3));
    for z=1:size(Mph,3)
        dets=det(Mph(:,:,z));
    end
    cost = cost + REL_DETERMINANT_COST * sum(abs(log(abs(dets))));
end

if nargout >=3
    warpedPhotos = phvals;
    grad=[];
end
if nargout >=4
    warpedMasks = mvals;
end
if nargout >=5
    REFvox2ras0New = T * REFmri.vox2ras0;
end

%%%%%%%%%%%%  GRADIENT %%%%%%%%%%%%%

if nargout == 2
    
    EPS=0.05;
    grad = zeros(size(params));
    
    % First, gradient of reference volume
    % One parameter at the time; we don't need to recompute deformed slices, so we provide those).
    if mode == 1 % similarity
        start = length(params)-6;
    else  % affine
        start = length(params)-11;
    end
    for j=start:length(params)
        paramsG=params;
        paramsG(j)=params(j)+EPS;
        costGplus=costFun(paramsG,cogREF,REFmri,Imri,Mmri,IIph,JJph,KKph,...
            REL_NCC_INTRA_WEIGHT,REL_DICE_INTRA_WEIGHT,REL_DICE_INTER_WEIGHT,...
            REL_DETERMINANT_COST, mode, phvals, mvals);
        paramsG=params;
        paramsG(j)=params(j)-EPS;
        costGminus=costFun(paramsG,cogREF,REFmri,Imri,Mmri,IIph,JJph,KKph,...
            REL_NCC_INTRA_WEIGHT,REL_DICE_INTRA_WEIGHT,REL_DICE_INTER_WEIGHT,...
            REL_DETERMINANT_COST, mode, phvals, mvals);        
        grad(j) = (costGplus-costGminus)/(2*EPS);
    end
    
    
    % Next, gradient of photos
    % Here, what we do is to provide a mini-problem with only 2-3 slices
    for n=1:Nims
        
        % Extract mini-problem
        n1=max(1,n-1);
        n2=min(Nims,n+1);
        
        ImriSl=Imri;
        ImriSl.vol=ImriSl.vol(:,:,n1:n2,:);
        ImriSl.vox2ras0(1:3,4)=ImriSl.vox2ras0(1:3,4)+ImriSl.vox2ras0(1:3,1:3)*[0;0;n1-1];
        
        MmriSl=ImriSl;
        MmriSl.vol=Mmri.vol(:,:,n1:n2);
        
        IIphSl=IIph(:,:,n1:n2);
        JJphSl=JJph(:,:,n1:n2);
        KKphSl=KKph(:,:,n1:n2);

        DiceInterAccumNumSl = sum(DiceInterNums([1:n1-1 n2+1:end]));
        DiceInterAccumDenSl = sum(DiceInterDens([1:n1-1 n2+1:end]));
        
        refvalsSl=refvals(:,:,n1:n2);
        
        % parameters to probe, in full size and in slice subproblem
        if mode==1 % rigid + similarity
            idx=3*n-2:3*n;
            idxSl=[3*n1-2:3*n2 length(params)-6:length(params)];
        elseif mode==2 % rigid + affine
            idx=3*n-2:3*n;
            idxSl=[3*n1-2:3*n2 length(params)-11:length(params)];
        else % affine + affine
            idx=[4*n-3:4*n 4*Nims+2*n-1 4*Nims+2*n];
            idxSl=[4*n1-3:4*n2 4*Nims+2*n1-1:4*Nims+2*n2 length(params)-11:length(params)];
        end
        paramsSl=params(idxSl);
        
%         costSl=costFun(paramsSl,cogREF,REFmri,ImriSl,MmriSl,IIphSl,JJphSl,KKphSl,...
%             REL_NCC_INTRA_WEIGHT,REL_DICE_INTRA_WEIGHT,REL_DICE_INTER_WEIGHT,...
%             REL_DETERMINANT_COST, mode, [], [], refvalsSl,DiceInterAccumNumSl,...
%             DiceInterAccumDenSl);
        
        % Compute gradients
        for i=idx
            j=find(idxSl==i);
            
            paramsSlG=paramsSl;
            paramsSlG(j)=paramsSlG(j)+EPS;
            costSlGplus=costFun(paramsSlG,cogREF,REFmri,ImriSl,MmriSl,IIphSl,JJphSl,KKphSl,...
                REL_NCC_INTRA_WEIGHT,REL_DICE_INTRA_WEIGHT,REL_DICE_INTER_WEIGHT,...
                REL_DETERMINANT_COST, mode, [], [], refvalsSl,DiceInterAccumNumSl,...
                DiceInterAccumDenSl);
            
            paramsSlG=paramsSl;
            paramsSlG(j)=paramsSlG(j)-EPS;
            costSlGminus=costFun(paramsSlG,cogREF,REFmri,ImriSl,MmriSl,IIphSl,JJphSl,KKphSl,...
                REL_NCC_INTRA_WEIGHT,REL_DICE_INTRA_WEIGHT,REL_DICE_INTER_WEIGHT,...
                REL_DETERMINANT_COST, mode, [], [], refvalsSl,DiceInterAccumNumSl,...
                DiceInterAccumDenSl);
            
            grad(i) = (costSlGplus-costSlGminus)/(2*EPS);
            
        end
    end
    
end



