function [B0param B0map] = LRB0map_KM(mapimage,imagesize,mapdelay)

    map = angle(conj(squeeze(mapimage(:,:,2))).*squeeze(mapimage(:,:,1)))/2/pi/mapdelay;
    mapmag = abs(mapimage(:,:,2));

% % calculate the B0param with the low resolution map; replaced by using the high resolution map with polynomial fit
% mapmags = mapmag.^2;
% x = repmat([-mapsize/2:1:mapsize/2-1],mapsize,1); % x points right, y points down, x corresponds to PE, y corresponds to RO
% y = repmat([-mapsize/2:1:mapsize/2-1]',1,mapsize);
% S = sum(sum(mapmags));
% Sx = sum(sum(x.*mapmags));
% Sy = sum(sum(y.*mapmags));
% Sf = sum(sum(map.*mapmags));
% Sxx = sum(sum(x.^2.*mapmags));
% Syy = sum(sum(y.^2.*mapmags));
% Sxy = sum(sum(x.*y.*mapmags));
% Sxf = sum(sum(x.*map.*mapmags));
% Syf = sum(sum(y.*map.*mapmags));
% delta = det([S Sx Sy; Sx Sxx Sxy; Sy Sxy Syy]);
% deltaf = det([Sf Sx Sy; Sxf Sxx Sxy; Syf Sxy Syy]);
% deltax = det([S Sf Sy; Sx Sxf Sxy; Sy Syf Syy]);
% deltay = det([S Sx Sf; Sx Sxx Sxf; Sy Sxy Syf]);
% 
% f0 = deltaf/delta;
% alpha = deltax/delta;
% beta = deltay/delta;
% 
% LRB0param = [f0 alpha*mapsize beta*mapsize]; % convert to Hz/FOV

[B0param B0map] = B0mapPolyFit(map,mapmag,imagesize);

