 function hfigout = showimage3(X, nbrk, lgsc,dyy,dxx)
% showimage3(X, nbrk, lgsc,dxx,dyy)
%    X = image
%    nbrk = number of panels for long images
%    lgsc = log scale (min to maax) in dB, negative numbers for linear 
%    dxx,dyy = spacing of samples (default = 1) or vector of positions

if nargin < 2
	help showimage3
	return
elseif nargin > 2
    if (lgsc >= 0)
        % take log of abs now
        xx = abs(X);
        mx = max(xx(:));
        X = 20*log10(xx/mx);
        X(find(X < -lgsc)) = -lgsc;  % handles -Inf as well
    end
end

sty = 0;
if nargin > 4
    if length(dyy) > 1
        dy = dyy(2) - dyy(1);
        sty = dyy(1) - dy;
    else
        dy = dyy;
    end
else
    dy = 1;
end
stx = 0;
if nargin > 3
    if length(dxx) > 1
        dx = dxx(2) - dxx(1);
        stx = dxx(1) - dx;
    else
        dx = dxx;
    end
else
    dx = 1;
end

nrow = nbrk;
ncol = 1;

[nnr nnc] = size(X);
subnc = ceil(nnc/ncol);
nextrac = subnc*ncol - nnc;
subnr = ceil(nnr/nrow);
nextrar = subnr*nrow - nnr;
warning('off')
padX = [X zeros([nnr nextrac]); zeros([nextrar nnc]) zeros([nextrac nextrar])];
warning('on')
mn = min(padX(:));
mx = max(padX(:));

for lpr = 1:nrow
    for lpc = 1:ncol
        subplot(ncol,nrow,lpr+(lpc-1)*ncol)
        vc = [1:subnc]+(lpc-1)*subnc;
        vr = [1:subnr]+(lpr-1)*subnr;
                imagesc(vc*dy + sty,vr*dx + stx,padX(vr,vc),[mn mx])
    end
end
colormap gray

