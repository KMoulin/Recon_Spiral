function [dKshift] = correction_traj2(RAW,G,K,RECON)


K2=K;
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);                    % Returns Zero-Crossing Indices Of Argument Vector
Gamma=42575.6;

    zx = zci(G);     % Approximate Zero-Crossing Indices
   % zx=zx(3:end-1);     % remove the first Echo
    
   
   vect_size=max(diff(zx));
   wndw=max(vect_size);
   
   %zx=zx(zx>wndw);
   
   VectX=nan(wndw,length(zx)-1);
   VectY=nan(wndw,length(zx)-1);
%%
   RAW2=RAW;
   r2_min=999999;
   dKshift=0;
   for cpt_dK=0:1:100
     dK=-0.0001*cpt_dK;
     K2=K+sign(G)*dK.*RECON.Time';
     
    for cpt=2:1:length(zx)
              
        vect_ind=zx(cpt-1):1:zx(cpt);       
        VectX(1:length(vect_ind),cpt-1)=K2(vect_ind);      
        VectY(1:length(vect_ind),cpt-1)=abs(RAW2(vect_ind)); 
        
        [ValM,IndEcho1]=max(abs(VectY(:,cpt-1)));
        Echo2(cpt-1,1)=VectX(IndEcho1,cpt-1);
        Echo2(cpt-1,2)=VectY(IndEcho1,cpt-1);
        Echo2(cpt-1,3)=zx(cpt-1)+IndEcho1;
        
    end
    
%     figure
%     subplot(3,1,1)
%     plot(Echo2(:,1),Echo2(:,2),'o')
%     subplot(3,1,2)
%     plot(G2),hold on,plot(Echo2(:,3),zeros(length(Echo2(:,3)),1),'o')
%     subplot(3,1,3)
%     plot(K2),hold on,plot(Echo2(:,3),zeros(length(Echo2(:,3)),1),'o')
    figure(3)
    plot(Echo2(:,1),Echo2(:,2),'o')
%     hold on
%     plot(RECON.Time'/2,RECON.Time'.*0.006.*sign(G),'o')
%     plot(K2,abs(RAW2))
    xx=xlim;
    yy=ylim;
    text(3*xx(1)/4,9.5*yy(2)/10,['Current Delay=' num2str(dK) 'us; Ref Delay=' num2str(dKshift) 'us']);
    axis([-2000 2000 0 yy(2)])
    % pause(0.5)
     if(sum(abs( Echo2(:,1)))<r2_min)
      r2_min=sum(abs( Echo2(:,1)));
      dKshift=dK;
     end  
   end
   
    K2=K+sign(G)*dKshift.*RECON.Time';
    figure,plot(K2,abs(RAW))
%     
%    
%     
%     figure
%     plot((Echo2(1:2:end,3)),abs(Echo2(1:2:end,1)),'o')
%     figure
%     plot((Echo2(2:2:end,3)),abs(Echo2(2:2:end,1)),'o')
    
    
end