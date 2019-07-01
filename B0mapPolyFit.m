function [B0param HRmap] = B0mapPolyFit(LRmap,LRmapmag,imagesize)
% function [B0param HRmap] = B0mapPolyFit(LRmap,LRmapmag,imagesize)
%
% fit a high resolution map from the low resolution map and get the linear
% parameter
%
% (c) Xue Feng 2012
% University of Virginia

ncoef = 21;
if (size(LRmap,1) ~= size(LRmap,2))
    error('only square map supported\n');
end
mapsize = size(LRmap,1);
f1 = imagesize/mapsize;
f2 = f1*f1;
f3 = f2*f1;
f4 = f3*f1;
f5 = f4*f1;

x = repmat([-mapsize/2:1:mapsize/2-1],mapsize,1); % x points right, y points down, x corresponds to PE, y corresponds to RO; this is important
y = repmat([-mapsize/2:1:mapsize/2-1]',1,mapsize);
F(:,1) = LRmapmag(:);
F(:,2) = LRmapmag(:).*x(:);
F(:,3) = LRmapmag(:).*y(:);
F(:,4) = LRmapmag(:).*x(:).*x(:);
F(:,5) = LRmapmag(:).*x(:).*y(:);
F(:,6) = LRmapmag(:).*y(:).*y(:);
F(:,7) = LRmapmag(:).*x(:).*x(:).*x(:);
F(:,8) = LRmapmag(:).*x(:).*x(:).*y(:);
F(:,9) = LRmapmag(:).*x(:).*y(:).*y(:);
F(:,10) = LRmapmag(:).*y(:).*y(:).*y(:);
F(:,11) = LRmapmag(:).*x(:).*x(:).*x(:).*x(:);
F(:,12) = LRmapmag(:).*x(:).*x(:).*x(:).*y(:);
F(:,13) = LRmapmag(:).*x(:).*x(:).*y(:).*y(:);
F(:,14) = LRmapmag(:).*x(:).*y(:).*y(:).*y(:);
F(:,15) = LRmapmag(:).*y(:).*y(:).*y(:).*y(:);
F(:,16) = LRmapmag(:).*x(:).*x(:).*x(:).*x(:).*x(:);
F(:,17) = LRmapmag(:).*x(:).*x(:).*x(:).*x(:).*y(:);
F(:,18) = LRmapmag(:).*x(:).*x(:).*x(:).*y(:).*y(:);
F(:,19) = LRmapmag(:).*x(:).*x(:).*y(:).*y(:).*y(:);
F(:,20) = LRmapmag(:).*x(:).*y(:).*y(:).*y(:).*y(:);
F(:,21) = LRmapmag(:).*y(:).*y(:).*y(:).*y(:).*y(:);
dmap = LRmapmag(:).*LRmap(:);

d = F\dmap;

HRx = repmat([-imagesize/2:1:imagesize/2-1],imagesize,1); % x points right, y points down, x corresponds to PE, y corresponds to RO; this is important
HRy = repmat([-imagesize/2:1:imagesize/2-1]',1,imagesize);
x1 = HRx;
x2 = x1.*x1;
x3 = x2.*x1;
x4 = x3.*x1;
x5 = x4.*x1;
y1 = HRy;
y2 = y1.*y1;
y3 = y2.*y1;
y4 = y3.*y1;
y5 = y4.*y1;

HRmap1 = (x1*d(2)+y1*d(3))/f1;
HRmap2 = (x2*d(4)+x1.*y1*d(5)+y2*d(6))/f2;
HRmap3 = (x3*d(7)+x2.*y1*d(8)+x1.*y2*d(9)+y3*d(10))/f3;
HRmap4 = (x4*d(11)+x3.*y1*d(12)+x2.*y2*d(13)+x1.*y3*d(14)+y4*d(15))/f4;
HRmap5 = (x5*d(16)+x4.*y1*d(17)+x3.*y2*d(18)+x2.*y3*d(19)+x1.*y4*d(20)+y5*d(21))/f5;

HRmap = d(1)+HRmap1+HRmap2+HRmap3+HRmap4+HRmap5;

HRmapmags = ones(imagesize);
S = sum(sum(HRmapmags));
Sx = sum(sum(HRx.*HRmapmags));
Sy = sum(sum(HRy.*HRmapmags));
Sf = sum(sum(HRmap.*HRmapmags));
Sxx = sum(sum(HRx.^2.*HRmapmags));
Syy = sum(sum(HRy.^2.*HRmapmags));
Sxy = sum(sum(HRx.*HRy.*HRmapmags));
Sxf = sum(sum(HRx.*HRmap.*HRmapmags));
Syf = sum(sum(HRy.*HRmap.*HRmapmags));
delta = det([S Sx Sy; Sx Sxx Sxy; Sy Sxy Syy]);
deltaf = det([Sf Sx Sy; Sxf Sxx Sxy; Syf Sxy Syy]);
deltax = det([S Sf Sy; Sx Sxf Sxy; Sy Syf Syy]);
deltay = det([S Sx Sf; Sx Sxx Sxf; Sy Sxy Syf]);

f0 = deltaf/delta;
alpha = deltax/delta;
beta = deltay/delta;

B0param = [f0 alpha*imagesize beta*imagesize]; % convert to Hz/FOV