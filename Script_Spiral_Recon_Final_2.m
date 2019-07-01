close all;
clear all;

file=uigetfile('.dat');
twix_obj = mapVBVD(file);

%% %%%%%%%%%%%%%%  EXTRACT RAW DATA INFORMATION %%%%%%%%%%%%%%

ICE=[];
ICE.BW=twix_obj{1, 2}.hdr.MeasYaps.sRXSPEC.alDwellTime{1};
ICE.DwellTime=ICE.BW/1000; % us
ICE.GradRasterTime=10; % us
ICE.dT=1e-6; % s
ICE.Gamma=42575.6; % Hz/mT from 42.576 MHz/T
ICE.dTwoPiPrecis = 6.2831853;
ICE.dGolden=2.3998277;%1.9416; % GoldenAngle =137.5
ICE.diff_nb_bval=twix_obj{1, 2}.hdr.MeasYaps.sDiffusion.lDiffWeightings;
ICE.diff_bval=cell2mat(twix_obj{1, 2}.hdr.MeasYaps.sDiffusion.alBValue);
ICE.diff_dir=twix_obj{1, 2}.hdr.MeasYaps.sDiffusion.lDiffDirections;
ICE.diff_avg=cell2mat(twix_obj{1, 2}.hdr.MeasYaps.sDiffusion.alAverages);

ICE.Grad_Delay=twix_obj{1,2}.image.iceParam(17,1); % us
ICE.Ramp_Down=twix_obj{1,2}.image.iceParam(18,1);  % us
ICE.Oversamp=twix_obj{1,2}.image.iceParam(19,1);
ICE.OversampSign=twix_obj{1,2}.image.iceParam(20,1); % 0 is negative
ICE.Nbleaves=twix_obj{1,2}.image.iceParam(21,1);
ICE.Grad_Ampl=twix_obj{1,2}.image.iceParam(22,1)/100; % mT/m
ICE.RiseTime=twix_obj{1,2}.image.iceParam(23,1)/100;  % us
ICE.RiseTimeFactor=twix_obj{1, 2}.hdr.MeasYaps.sWipMemBlock.adFree{2};
ICE.FOV=twix_obj{1, 2}.hdr.Config.RoFOV/1000; % mm -> m
ICE.Matrix=twix_obj{1, 2}.hdr.Config.NImageCols;
ICE.TE1= twix_obj{1, 2}.hdr.MeasYaps.alTE{1}/1000;
ICE.TE2= twix_obj{1, 2}.hdr.MeasYaps.sWipMemBlock.adFree{1};


if(ICE.OversampSign)ICE.FOV1=(ICE.Oversamp/100)*ICE.FOV;
else ICE.FOV1=(-ICE.Oversamp/100)*ICE.FOV;
end
ICE.Kmax=ICE.Matrix/(2*ICE.FOV);	
ICE.SlewRate= 1/(ICE.RiseTime*ICE.dT); % mT/(m.s)
ICE.Additional_Delay=0;
ICE.TE_Delay_KPoint=round(ICE.Grad_Delay/ICE.DwellTime)+ICE.Additional_Delay; % Points
ICE.Ramp_Down_KPoint=round(ICE.Ramp_Down*ICE.GradRasterTime/ICE.DwellTime);  % Points
ICE.Seg=twix_obj{1, 2}.image.NSeg;

%% %%%%%%%%%%%%%% GENERATE K-SPACE TRAJECTORY FROM PROTOCOL %%%%%%%%%%%%%%
[K,G,K_ADC,G_ADC]=calc_vds(ICE); % K: K-space trajectory every 1 us, K_ADC: K-space trajectory "interpolated" to dwell Times, G_ADC: Gradients interpolated to dwell Times, G: Gradients at every grad raster Times
ICE.ADCPoint=size(K_ADC,1);
K_ADC2= cumsum(G_ADC,1);

%% %%%%%%%%%%%%%% RECON PARAMETER %%%%%%%%%%%%%%

RECON=[];
RECON.MatrixOverSampling=1;
RECON.Delay=[1,1];  % Gradient delay units in points
RECON.Kshift=[0,0];
RECON.Matrix=ICE.Matrix*1;
RECON.Debug=0; % 1 for display the plot
RECON.Nbleaves=twix_obj{1, 2}.hdr.MeasYaps.sDiffusion.alAverages{1};
RECON.Time=0:ICE.DwellTime:(ICE.ADCPoint-ICE.TE_Delay_KPoint+ICE.Additional_Delay)*ICE.DwellTime;
RECON.VectorDat=ICE.TE_Delay_KPoint-ICE.Additional_Delay:ICE.ADCPoint;
RECON.B0_dTE=(ICE.TE1-ICE.TE2)*1e-3;

%% %%%%%%%%%%%%%% Reorder the RAW data %%%%%%%%%%%%%%
% RAW_diff[points, coils, b-values, directions, averages, shots];
%

RAW_tmp=squeeze(twix_obj{1,2}.image.unsorted);
RAW=[];
for cpt_seg=1:1:ICE.Seg
    RAW=[RAW; RAW_tmp(:,:,cpt_seg:ICE.Seg:end)]; 
end

% RAWdata 1 to 4 are pre scans:
% 1 & 2 are for estimating the gradient delay in X and Y
% 3 & 4 are a single shot B0 map with short and long TE
RAW_PreScan=RAW(:,:,1:4);  

cpt2=5; % Start RAWdata counter after the pre scans
cpt_shot=1;
cpt_avg=1;
for cpt_global_avg=1:1:min(ICE.diff_avg) % For now toss low TE data
    if cpt_shot>ICE.Nbleaves
        cpt_shot=1;
        cpt_avg=cpt_avg+1;
    end
        
    for cpt_b=1:1:ICE.diff_nb_bval
            if (cpt_b>1)
                for cpt_dir=1:1:ICE.diff_dir
                    RAW_diff(:,:,cpt_b,cpt_dir,cpt_avg,cpt_shot)=RAW(RECON.VectorDat,:,cpt2);
                    cpt2=cpt2+1;
                end
            else
                 RAW_diff(:,:,cpt_b,1,cpt_avg,cpt_shot)=RAW(RECON.VectorDat,:,cpt2);       
                 cpt2=cpt2+1;                 
            end
    end
    cpt_shot=cpt_shot+1;
end
clear cpt_avg cpt_b cpt_dir cpt_shot cpt_global_avg cpt2;

%% %%%%%%%%%%%%%% ESTIMATE GRAD DELAY CORRECTION AND APPLY IT (K_ADC2 -> K_cor) %%%%%%%%%%%%%%

[RECON.Delay(1)]=correction_kdelay(G_ADC(RECON.VectorDat,1),K_ADC2(RECON.VectorDat,1),abs(RAW_PreScan(RECON.VectorDat,1,1)));
[RECON.Delay(2)]=correction_kdelay(G_ADC(RECON.VectorDat,2),K_ADC2(RECON.VectorDat,2),abs(RAW_PreScan(RECON.VectorDat,1,2)));

[G_cor, K_cor] = correction_traj_delay(G_ADC(RECON.VectorDat,:),K_ADC2(RECON.VectorDat,:),RECON);

%% %%%%%%%%%%%%%% (WIP) ESTIMATE TRAJECTORY DISTORTION BUT DO NOT APPLY IT (K_cor -> K_cor2) %%%%%%%%%%%%%%

[RECON.Kshift(1)] = correction_traj2(RAW_PreScan(RECON.VectorDat,1,1),G_cor(:,1),K_cor(:,1),RECON);
[RECON.Kshift(2)] = correction_traj2(RAW_PreScan(RECON.VectorDat,1,2),G_cor(:,2),K_cor(:,2),RECON);

[K_cor2] = correction_traj_shift2(K_cor,G_cor,RECON);

if RECON.Debug>0 % Display Correction (for debug)
    figure,plot(K_cor(:,1),abs(RAW_PreScan(RECON.VectorDat,1,1))),hold on,%plot(-K_cor(:,1),abs(RAW_PreScan(ICE.TE_Delay_KPoint:ICE.ADCPoint,1,2)))
    axis([ -1600 1600 0 0.035]);
    figure,plot(K_cor(:,2),abs(RAW_PreScan(RECON.VectorDat,1,2))),hold on,%plot(-K_cor(:,2),abs(RAW_PreScan(ICE.TE_Delay_KPoint:ICE.ADCPoint,1,4))),
    axis([ -1600 1600 0 0.035]);

    figure,plot(K_cor2(:,1),abs(RAW_PreScan(RECON.VectorDat,1,1))),hold on,%plot(-K_cor(:,1),abs(RAW_PreScan(ICE.TE_Delay_KPoint:ICE.ADCPoint,1,2)))
    axis([ -1600 1600 0 0.035]);
    figure,plot(K_cor2(:,2),abs(RAW_PreScan(RECON.VectorDat,1,2))),hold on,%plot(-K_cor(:,2),abs(RAW_PreScan(ICE.TE_Delay_KPoint:ICE.ADCPoint,1,4))),
    axis([ -1600 1600 0 0.035]);
end
K_cor2=K_cor; % As it's still WIP we do not apply it for now. 
 
%% %%%%%%%%%%%%%% DENSITY COMPENSATION CORRECTION (Work well for linear density but not too well for advance density sampling) %%%%%%%%%%%%%%

RECON.Dcf=[];
RECON.Dcf_s=[0; abs(diff(complex(K_cor2(:,1),K_cor2(:,2))))/max(abs(diff(complex(K_cor2(:,1),K_cor2(:,2)))))];

%% %%%%%%%%%%%%%% (WIP) PHASE CORRECTION (Doesn't work with NUFFT) %%%%%%%%%%%%%%


for cpt_diff=1:1:size(RAW_diff,3) % For now take only the first direction and first average of diffusion image and recon it
        [RAW_tmp]= phase_correction(squeeze(RAW_diff(:,:,cpt_diff,1,:)),K_cor2,RECON);   
        RAW_diff2(:,:,cpt_diff,1,:)=RAW_tmp;
end
RAW_diff2=RAW_diff; % As it's still WIP we do not apply it for now. 

%% %%%%%%%%%%%%%% (WIP) T2 STAR CORRECTION %%%%%%%%%%%%%%

[RECON.T2Sval]=correction_T2star(K_cor2,abs(RAW_PreScan(RECON.VectorDat,1,1:2)),ICE.DwellTime); % us
RECON.T2cf_s=exp(-RECON.Time./RECON.T2Sval)';
%RECON.Dcf_s=RECON.Dcf_s./RECON.T2cf_s;  % As it's still WIP we do not apply it for now. 

%% %%%%%%%%%%%%%% GENERATE K-SPACE trajectory for every shot (From K_cor -> K_rot) %%%%%%%%%%%%%%

for cpt=1:1:RECON.Nbleaves
    %ICE.dRotAngle(cpt) = ICE.dTwoPiPrecis *(cpt-1)/(ICE.Nbleaves); % Fix Angle version
    ICE.dRotAngle(cpt) = ICE.dGolden *(cpt-1);                    % Golden Angle version
    ICE.RotMat(:,:,cpt) = [cos(ICE.dRotAngle(cpt)) -sin(ICE.dRotAngle(cpt));sin(ICE.dRotAngle(cpt)) cos(ICE.dRotAngle(cpt))];
    
    for cpt2=1:1:length(K_cor2(:,1))
        K_rot(cpt2,cpt,:)=K_cor2(cpt2,:)*squeeze(ICE.RotMat(:,:,cpt));
    end
    
    for cpt2=1:1:length(G_cor)
        G_rot(cpt2,cpt,:)=G_cor(cpt2,:)*squeeze(ICE.RotMat(:,:,cpt));
    end
    
end

%% %%%%%%%%%%%%%% GLOBAL TRAJECTORY AND NORMALIZATION (From K_rot -> K1) %%%%%%%%%%%%%%

K1=complex(squeeze(K_rot(:,:,1)),squeeze(K_rot(:,:,2)));
K1=RECON.MatrixOverSampling*K1./(2*max(max(abs(K1)))); % Normalization

%% %%%%%%%%%%%%%% RECONSTRUCT PRESCAN 3 & 4 TO GENERATE B0 MAP (From K1 & RAW_PreScan -> Grid_K) %%%%%%%%%%%%%%

RECON.Dcf=RECON.Dcf_s;
Final_eMap=[];
Final_ePhase=[];
[Final_eMap(:,:,1) S]= nufft_n_combine(RAW_PreScan(RECON.VectorDat,:,3),K1(:,1)*4,RECON,RECON.Matrix/4,1);
[Final_eMap(:,:,2) S]= nufft_n_combine(RAW_PreScan(RECON.VectorDat,:,4),K1(:,1)*4,RECON,RECON.Matrix/4,1);
Final_eMap(isnan(Final_eMap))=0;
Final_ePhase(:,:,1)=phase_unwrap(angle(Final_eMap(:,:,1)));
Final_ePhase(:,:,2)=phase_unwrap(angle(Final_eMap(:,:,2)));

[B0param B0map] = LRB0map_KM(Final_eMap,RECON.Matrix,RECON.B0_dTE);

[RAW_diff3, K3] = correctionB0(RAW_diff2,K1,RECON,B0param);

%% %%%%%%%%%%%%%% REGRIDDING IMG (From K1 & Raw_diff -> Grid_K) %%%%%%%%%%%%%%

for cpt_diff=1:1:size(RAW_diff3,3) % Take only the first direction and first average of diffusion image and recon it
    RAW_tmp=[];
    K_tmp=[];
    RECON.Dcf=[];
    disp('Recon Img');
    for cpt_lv=1:1:RECON.Nbleaves 
%          RAW_tmp=[];
%         K_tmp=[];
%         RECON.Dcf=[];
        RAW_tmp=[RAW_tmp ;squeeze(RAW_diff3(:,:,cpt_diff,1,cpt_lv))];    
        RECON.Dcf=[RECON.Dcf ;RECON.Dcf_s];
        K_tmp=[K_tmp ;K3(:,cpt_lv)];  
%         [Final_Img(:,:,cpt_diff,cpt_lv) S]= nufft_n_combine(RAW_tmp,K_tmp,RECON,RECON.Matrix);
%         Final_Phase(:,:,cpt_diff,cpt_lv)=phase_unwrap(angle(Final_Img(:,:,cpt_diff,cpt_lv)));
    end
    %[Final_Img(:,:,cpt_lv,cpt_diff) S]= grid_n_combine(RAW_tmp,K_tmp,RECON,RECON.Matrix);
    [Final_Img(:,:,cpt_diff) S]= nufft_n_combine(RAW_tmp,K_tmp,RECON,RECON.Matrix,2);
      Final_Img(isnan(Final_Img))=0;
      Final_Phase(:,:,cpt_diff)=phase_unwrap(angle(Final_Img(:,:,cpt_diff)));
   
end
clear cpt cpt_coil cpt_diff cpt_lv;

%% %%%%%%%%%%%%%% DISPLAY FINAL IMAGES & ADC CALCULATION %%%%%%%%%%%%%%

figure,imagesc(abs(Final_Img(:,:,1))),colormap('gray');
figure,imagesc(Final_Phase(:,:,1));

figure,imagesc(abs(Final_Img(:,:,2))),colormap('gray');
figure,imagesc(Final_Phase(:,:,2));


ADC=abs(log(abs(Final_Img(:,:,2))./abs(Final_Img(:,:,1)))/-350);
figure,imagesc(ADC,[0 0.005]),colormap('jet') %% Everything is too high by a factor 2....



%% %%%%%%%%%%%%%%  LOCAL FUNCTIONS %%%%%%%%%%%%%%
function S = adaptive_coil(data)
    data_2d(:,:,1,:)=data;
    [Nx,Ny,Nz,Nc] = size(data_2d);
    S = zeros(Nx,Ny,Nz,Nc);
    M = zeros(Nx,Ny,Nz);
    w = 5;
    for i = 1:Nx
        ii = max(i-w,1):min(i+w,Nx);
        for j = 1:Ny
            jj = max(j-w,1):min(j+w,Ny);
            for k = 1:Nz
                kk = max(k-w,1):min(k+w,Nz);
                kernel = reshape(data_2d(ii,jj,kk,:),[],Nc);
                [V,D] = eigs(conj(kernel'*kernel),1);
                S(i,j,k,:) = V*exp(-1j*angle(V(1)));
                M(i,j,k) = sqrt(D);
            end
        end
    end
    S = squeeze(S.*(M>0.01*max(abs(M(:)))));
end

function [Img, S] = grid_n_combine(RAW,K,RECON,Matrix,CoilMode)
    
    Grid_K=[];  
    S=[];
     
     %% Regridding 
     disp('Coil');
     h = waitbar(0,'Coil...');
    for cpt_coil=1:1:size(RAW,2)
        Data=RAW(:,cpt_coil); % Take the acquisition part after grad delay;
        Data=Data(:);
        K=K(:);
        
       [tmp_dat m] = gridmat(K,Data,RECON.Dcf,Matrix);
        Grid_K(:,:,cpt_coil)=tmp_dat;
        waitbar(cpt_coil/size(RAW,2),h);
    end   
    close(h)
    
    %% FFT (Grid_K -> Grid_Img)
    Grid_Img=fftshift(fft2(fftshift(Grid_K)));
    Grid_Img(isnan(Grid_Img))=0;
    
    %% Deconvolution
    Grid_Img=Grid_Img./m; 
    
    %% % COIL COMBINATION (From Grid_Img -> Final_Img): Sum of square of the magnitude
    
    %% % COIL COMBINATION (From Grid_Img -> Final_Img): Sum of square of the magnitude
    if (CoilMode==2) % Adaptive coil combine (eSPIRIT), works well with a pre-coil selection (WIP)
        [val index]=maxk(squeeze(median(median(abs(Grid_Img(:,:,:))))),12);
        S = adaptive_coil(Grid_Img(:,:,index));
        Img=sum(Grid_Img(:,:,index).*conj(S),3)./sum(S.*conj(S),3);
    else% Sum of square
        Img=sqrt(sum((Grid_Img).^2,3));   
    end   
    
end

function [Img, S] = nufft_n_combine(RAW,K,RECON,Matrix,CoilMode)  
     
     Grid_Img=[];
     S=[];
     
     %% Nufft
     st = nufft_init([2*pi*real(K) 2*pi*imag(K)] , [Matrix Matrix], [5 5], [Matrix*2 Matrix*2],[Matrix/2 Matrix/2], 'minmax:kb');  % Basic parameters
     disp('Coil');
     h = waitbar(0,'Coil...');
    for cpt_coil=1:1:size(RAW,2)
        Data=complex(double(real(RAW(:,cpt_coil))), double(imag(RAW(:,cpt_coil)))).*RECON.Dcf ; % Take the acquisition part after grad delay;    
        Grid_Img(:,:,cpt_coil)  = nufft_adj(Data,st); 
        waitbar(cpt_coil/size(RAW,2),h);
    end   
    close(h)
    
    %% % COIL COMBINATION (From Grid_Img -> Final_Img): Sum of square of the magnitude
    if (CoilMode==2) % Adaptive coil combine (eSPIRIT), works well with a pre-coil selection (WIP)
        [val index]=maxk(squeeze(median(median(abs(Grid_Img(:,:,:))))),12);
        S = adaptive_coil(Grid_Img(:,:,index));
        Img=sum(Grid_Img(:,:,index).*conj(S),3)./sum(S.*conj(S),3);
    else% Sum of square
        Img=sqrt(sum((Grid_Img).^2,3));   
    end
end


function [RAW2, K2] = correctionB0(RAW,K,RECON,B0)   
    RAW2=zeros(size(RAW));
    K2=K;
    Time = repmat(RECON.Time',1,size(RAW,2),size(RAW,3),size(RAW,4),size(RAW,5))*1e-6;
    Time2 = repmat(RECON.Time',1,size(K,2))*1e-6;
    
    %% Off resonance phase correction on the data, shift the data based on the B0 shift (very simple)
    B0Cor=exp(i*2*pi*Time*B0(1));
    RAW2 =RAW.*B0Cor;  
    
    %% Off resonance phase correction on the trajectoy, shift the trajectory based on the B0 map (very simple)
    B0traj = B0(2).*Time2+i*B0(3).*Time2;
    K2 = K+B0traj;
end

function [G2, K2] = correction_traj_delay(G,K,RECON)
    G2=zeros(size(G));
    K2=zeros(size(K));
    for cpt=1:1:2
        if RECON.Delay(cpt)<0
            K2(1:end-abs(RECON.Delay(cpt))+1,cpt)=K(abs(RECON.Delay(cpt)):end,cpt);
            G2(1:end-abs(RECON.Delay(cpt))+1,cpt)=G(abs(RECON.Delay(cpt)):end,cpt);
        elseif RECON.Delay(cpt)>0
            K2(abs(RECON.Delay(cpt)):end,cpt)=K(1:end-abs(RECON.Delay(cpt))+1,cpt);
             G2(abs(RECON.Delay(cpt)):end,cpt)=G(1:end-abs(RECON.Delay(cpt))+1,cpt);
        else
            K2(:,cpt)=K(:,cpt);
            G2(:,cpt)=G(:,cpt);
        end
    end
end

function [RAW2] = correction_traj_shift(RAW,G,RECON)
    
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);                    % Returns Zero-Crossing Indices Of Argument Vector
RAW2=RAW;
for cpt_ax=1:1:size(G,2) 
    zx = zci(G(:,cpt_ax));     % Approximate Zero-Crossing Indices
    for cpt_coil=1:1:size(RAW,2)
        for cpt=2:1:length(zx)
            vect_ind=zx(cpt-1):1:zx(cpt);
            Kshift=vect_ind'*RECON.Kshift(cpt_ax);%
            R_tmp=RAW2(vect_ind,cpt_coil);           
            RAW2(vect_ind,cpt_coil)=ifft(((fft(R_tmp)).*exp(1i*2*pi*Kshift)));
            %RAW2(vect_ind,cpt_coil)=ifft((fft(R_tmp)));
        end
    end
end
end

function [K2] = correction_traj_shift2(K,G,RECON)   

    for cpt_ax=1:1:size(G,2) 
       K2(:,cpt_ax)=K(:,cpt_ax)+sign(G(:,cpt_ax))*RECON.Kshift(cpt_ax).*RECON.Time';
    end
end

function [RAW2]= phase_correction(RAW,K,RECON)

    
    % Normalize the trajectory  
    K(:,1)=K(:,1)./(2*max(abs(K(:,1))));
    K(:,2)=K(:,2)./(2*max(abs(K(:,2))));
    
    % Apply the Dcf on the data 
    RAW1=RAW.*RECON.Dcf_s;
    %RAW1=RAW;
    
    % Generate a gaussian window
    w = window(@gausswin,size(RAW1,1)*2,5);
    w=w(end/2+1:end);
    
    % Init
    RAW2=[];
    RAW_tmp=[];
    st = nufft_init([2*pi*K(:,1) 2*pi*K(:,2)] , [RECON.Matrix RECON.Matrix], [5 5], [RECON.Matrix*2 RECON.Matrix*2],[RECON.Matrix/2 RECON.Matrix/2], 'minmax:kb'); %kb  
    
    % For each data substract the low resolution phase image 
    for cpt_shot=1:1:size(RAW1,3)   
         for cpt_coil=1:1:size(RAW1,2) 
             
            RAW_tmp=double(squeeze(RAW1(:,cpt_coil,cpt_shot)).*w);

            IMG_triangle=nufft_adj(RAW_tmp,st);

            IMG=nufft_adj(double(squeeze(RAW1(:,cpt_coil,cpt_shot))),st);

            IMG2=abs(IMG).*exp(i*( angle(IMG)-angle(IMG_triangle)));
            IMG2_final(:,:,cpt_coil,cpt_shot)=IMG2;
            RAW2(:,cpt_coil,cpt_shot)=nufft(IMG2,st);
         end 
    end
 
end