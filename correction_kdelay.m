
function [Delay]=correction_kdelay(gdat,kdat,dat)

  


    Delay=0;
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);                    % Returns Zero-Crossing Indices Of Argument Vector


    zx = zci(gdat);     % Approximate Zero-Crossing Indices
   % zx=zx(3:end-1);     % remove the first Echos
    
   
   vect_size=max(diff(zx));
   wndw=max(vect_size);
   
   %zx=zx(zx>wndw);
   
   VectX=nan(wndw,length(zx));
   VectY=nan(wndw,length(zx));
   kdat2=kdat;

    for cpt=2:1:length(zx)
        vect_ind=zx(cpt-1):1:zx(cpt);
        VectX(vect_ind,cpt)=kdat(vect_ind);
        VectY(vect_ind,cpt)=dat(vect_ind);
    end
    %figure(1),plot(VectX(:,1:20),VectY(:,1:20))
    r2_min=9999999;
    Pos_Echo=[];
    k=1;
    MaxE=min(length(zx),20);
    for cpt_D=-20:1:20
        
       VectX=nan(wndw,length(zx));
       VectY=nan(wndw,length(zx));
       if cpt_D<0
            kdat2(1:end-abs(cpt_D)+1)=kdat(abs(cpt_D):end);
       elseif cpt_D>0
            kdat2(abs(cpt_D):end)=kdat(1:end-abs(cpt_D)+1);
       else
            kdat2=kdat;
       end 
       
       for cpt=2:1:length(zx)
        vect_ind=zx(cpt-1):1:zx(cpt);
        VectX(vect_ind,cpt)=kdat2(vect_ind);
        VectY(vect_ind,cpt)=dat(vect_ind);
       end 
        
       figure(1)
       plot(VectX(:,1:MaxE),VectY(:,1:MaxE))
       xx=xlim;
       yy=ylim;
       text(3*xx(1)/4,9.5*yy(2)/10,['Current Delay=' num2str(cpt_D) 'us; Ref Delay=' num2str(Delay) 'us']);
       [ValM,Indx]=max(VectY(:,1:MaxE));
       for cpt=1:1:MaxE
        Pos_Echo(cpt)=(VectX(Indx(cpt),cpt));
       end
       

       if(sum(abs(Pos_Echo))<r2_min)
          r2_min=sum(abs(Pos_Echo));
          Delay=cpt_D;
       end
       k=k+1;
  
    end
end
