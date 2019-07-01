function [T2Sval]=correction_T2star(kdat,dat,dwelltime)

    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);                    % Returns Zero-Crossing Indices Of Argument Vector

    kdat2=kdat-0.0000001;
    zx = zci(kdat2(:,1));     % Approximate Zero-Crossing Indices
    zx=zx(3:end-2);
    
    xdat=abs(dat(:,1,1));
    wndw=20;
    
    
    for cpt=1:1:length(zx)
        tmp=zeros(size(xdat));
        tmp(zx(cpt)-wndw:zx(cpt)+wndw)=xdat(zx(cpt)-wndw:zx(cpt)+wndw);
        [echox(cpt) indx(cpt)]=max(tmp);
    end

    
   % echox(:)=echox(:)./echox(1);
    
    
   lb = [1 1]; % S0, f, D : Min   
   ub = [20000 2000]; % Max
   val0 = [4000 200]; % init du fit 0
   options = optimset('Display','off');

   fun_relax_T2s = @(x,xdata)exp(-xdata/x(1));
   [val,resnorm,residual,exitflag] = lsqcurvefit(fun_relax_T2s,val0, indx*dwelltime, echox/max(echox), lb, ub, options);  % exp(-T2/TE)  

   T2Sval = val(1);  % us
  % T2S0 = val(2);
%     zy = zci(kdat(:,2)); 
%     zy=zy(3:end-2);
%     
%     ydat=abs(dat(:,3));
%     for cpt=1:1:length(zy)
%         ydat2=ydat(zy(cpt)-wndw:zy(cpt)+wndw);
%         d(cpt) = -sign(kdat(zx(cpt)-wndw,2))*finddelay(ydat2,ydatn2);
%     end
% 
%     Delayy=round(mean(d));

end