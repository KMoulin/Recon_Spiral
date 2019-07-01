
function [K,G,K_ADC,G_ADC]=calc_vds(ICE)

	 kr =0;	    % Current value of kr	
	 kr1=0;		% Current value of 1st derivative of kr 
	 kr2=0;		% Current value of 2nd derivative of kr 

	theta =0;		% Current value of theta 
	theta1=0;		% Current value of 1st derivative of theta 
	theta2=0;		% Current value of 2nd derivative of theta 

		
	% Calculate the k-space path 
    k=1;
    p=0;
    kx=[];
    ky=[];
	kx(1)=0;   % x components of current k-location 
    ky(1)=0;   % y components of current k-location 
	while (ICE.Grad_Delay>k)
		kx(k+1)=0;   % x components of current k-location 
        ky(k+1)=0;   % y components of current k-location 
		k=k+1;
		if(rem(k,ICE.GradRasterTime)==0)		
			p=p+1;
        end
    end
  
	m_lKPoint=k;	
	m_lGradPoint=p;
	
	while (kr < ICE.Kmax)
        
		[theta2,kr2] =calcthetadotdot(ICE.SlewRate,ICE.Grad_Ampl,kr,kr1,ICE.GradRasterTime,ICE.dT,ICE.Nbleaves, [ICE.FOV ICE.FOV1] ,2,  k, ICE.Kmax, ICE.Gamma);

		% Integrate to obtain new values of kr, kr1, theta and theta1 */

		theta1 = theta1 + theta2 * ICE.dT;
		theta  = theta  + theta1 * ICE.dT;

		kr1 = kr1 + kr2 * ICE.dT;
		kr  = kr  + kr1 * ICE.dT;
       
        kx(k+1)=(kr * cos(theta));
		ky(k+1)=(kr * sin(theta));

		k=k+1;
		if(rem(k,ICE.GradRasterTime)==0)		
			p=p+1;
        end
    end
    m_lKPoint=k;
	m_lGradPoint=p;
	
    [gx, gy]=grad_fact_spiral(ICE, kx,ky,m_lGradPoint,m_lKPoint);
   
    m_lGradPoint=m_lGradPoint+ICE.Ramp_Down+1;
    t=0:1:m_lKPoint-1; % 1us dwell time;
    t_dwell=0:ICE.DwellTime:(m_lKPoint-1); % ADC dwell time;
    
    t_g=0:ICE.GradRasterTime:(m_lGradPoint-1)*ICE.GradRasterTime; % 10us dwell time;
    t_g_dwell=0:ICE.DwellTime:(m_lGradPoint-1)*ICE.GradRasterTime; % Gradient ADC dwell time;
    
    
    kx_ADC = interp1(t,kx,t_dwell);
    ky_ADC = interp1(t,ky,t_dwell);
    
    
    if(length(t_g)>length(gx))
        gx=[0 gx];
        gy=[0 gy];
    end
    
    gx_ADC = interp1(t_g,gx,t_g_dwell);
    gy_ADC = interp1(t_g,gy,t_g_dwell);
   % G_ADC= 2*[0, 0 ;diff(K_ADC)/((ICE.Gamma*ICE.GradRasterTime*ICE.dT))];

   
%     for k=1:1:m_lKPoint        
%             if(rem((k-1),round(t))==0)
%                 K_ADC(p,:)=K(k,:);
%                 p=p+1;
%                 t=t+ICE.DwellTime;
%             end
%     end

    m_lADCPoint=length(t_dwell);
     
    G=[gx; gy]';
    K=[kx; ky]';
    K_ADC=[kx_ADC; ky_ADC]';
    G_ADC=[gx_ADC; gy_ADC]';
end


function [theta2,kr2] =calcthetadotdot( slewmax, gradmax, kr, kr1,  Tgsample, m_dT, Ninterleaves,  fov, numfov, k, m_dKrmax , m_dGammaH  )

	 dPI = 3.141592653589793;
     dTwoPiPrecis = 6.2831853;	
	
	
	 fovval=0;      % FOV for this value of kr	
	 dfovdrval=0;	% dFOV/dkr for this value of kr	
% 	 gmaxfov;		% FOV-limited Gmax.	
% 	 maxkr1;
% 	 count;
% 
% 	 tpf;     % Used to simplify expressions. */
% 	 tpfsq;	% 	" 		"        */
% 
% 	 qdfA, qdfB, qdfC;	%nQuadratic formula coefficients */
% 	 rootparta,rootpartb;


	for (count=1:1:numfov)
	
		fovval = fovval + fov(count)*(kr/m_dKrmax)^(count-1);
		if (count > 1)
            dfovdrval = dfovdrval + (count-1)*fov(count)*(kr/m_dKrmax)^(count-2)/m_dKrmax;
        end
    end
  

	gmaxfov = 1/( m_dGammaH * fovval * m_dT) ;	 %Tdsample
	if (gradmax > gmaxfov)
        gradmax = gmaxfov;
    end
 
	
	% Maximum dkr/dt, based on gradient amplitude. 

	maxkr1 = sqrt((m_dGammaH*gradmax)^2 / (1+(dTwoPiPrecis*fovval*kr/Ninterleaves)^2));
    % maxkr1 = sqrt(pow(m_dGammaH*gradmax,2) / (1+pow(dTwoPiPrecis*fovval*kr/Ninterleaves,2)));

 
	% These two are just to simplify expressions below 
	tpf = dTwoPiPrecis*fovval/Ninterleaves;
	tpfsq = tpf^2;
    
	if (kr1 > maxkr1)	% Then choose krdotdot so that krdot is in range 	
		kr2 = (maxkr1 - kr1)/(m_dT);
	else			% Choose krdotdot based on max slew rate limit. 
			% Set up for quadratic formula solution. 

		qdfA = 1+tpfsq*kr*kr;
		qdfB = 2*tpfsq*kr*kr1*kr1 + 2*tpfsq/fovval*dfovdrval*kr*kr*kr1*kr1;
		qdfC = (tpfsq*kr*kr1*kr1)^2 + 4*tpfsq*(kr1^4) + (tpf*dfovdrval/fovval*kr*kr1*kr1)^2 + 4*tpfsq*dfovdrval/fovval*kr*kr1^4 - (m_dGammaH*slewmax)^2;


		rootparta = -qdfB/(2*qdfA);
		rootpartb = qdfB*qdfB/(4*qdfA*qdfA) - qdfC/qdfA;
		
		% Safety check - if complex, take real part.
		if (rootpartb < 0)	
            kr2 = rootparta;
        else
            kr2 = rootparta + sqrt(rootpartb);
        end
     end

	% Calculate thetadotdot */
      theta2 = tpf*dfovdrval/fovval*kr1*kr1 + tpf*(kr2);
end



function [m_pGx, m_pGy]=grad_fact_spiral(ICE, kx,ky,lGradPoint,lKPoint)

    % Ramp of the gradients
	lSpiralRampPoint=ICE.Ramp_Down;
	
    % Allocate memory for gradients. */
    % m_pGx= new float[m_lGradPoint+m_lSpiralRampPoint];
    % m_pGy= new float[m_lGradPoint+m_lSpiralRampPoint];

	 MaxAmplX=0;
	 MaxAmplY=0;

	for p=1:1:(lGradPoint+lSpiralRampPoint)	
		m_pGx(p)=0;
		m_pGy(p)=0;
    end
   
	p=2;
	for k=1:1:lKPoint-1
	
		if(rem(k,ICE.GradRasterTime)==0)
			m_pGx(p) = ((kx(k+1)-kx(k-ICE.GradRasterTime+1))/(ICE.Gamma*ICE.GradRasterTime*ICE.dT));
			m_pGy(p) = ((ky(k+1)-ky(k-ICE.GradRasterTime+1))/(ICE.Gamma*ICE.GradRasterTime*ICE.dT));

			if(abs(m_pGx(p))>MaxAmplX) 
                MaxAmplX=abs(m_pGx(p));
            end
			if(abs(m_pGy(p))>MaxAmplY) 
                MaxAmplY=abs(m_pGy(p));
            end
            p=p+1; 	
        end
    end

	% Add a Ramp down
    for k=0:1:(lSpiralRampPoint-1)
        m_pGx(p)=(m_pGx(lGradPoint-1)-m_pGx(lGradPoint-1)*k/lSpiralRampPoint);
        m_pGy(p)=(m_pGy(lGradPoint-1)-m_pGy(lGradPoint-1)*k/lSpiralRampPoint);
        p=p+1;
    end

end






