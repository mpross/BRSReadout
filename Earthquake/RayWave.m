close all

%6.7 Vanuatu 4/6/16, 7.2 New Caledonia 8/12/16, 7.1 Atlantic 8/29/16, 
%7.8 New Zealand 11/13/16, 7.9 Papua New Guinea 12/17/16, 
%7.9 Papua New Guinea 1/22/17
fileName={'GPS1143961145_6_7_Vanuatu.mat','GPS1155000413_7_2_NewCaledonia.mat',...
    'GPS1156480214_7_1_Atlantic.mat','GPS1163070193_7_8_NewZealand.mat','GPS1166007087_7_9_PapuaNewGuinea.mat'...
    ,'GPS1169094640_7_9_PapuaNewGuinea2.mat'};

sampf =8;
noiseThreshold=0.5; % approx 5*BRS noise
phaseThreshold=20; % threshold for phase between seismometer signal and BRS signal

avgAng=[];
avgErrAng=[];

for m=0:5
    Astop1 = 20;
    Apass  = .5;
    Astop2 = 20;
    velA=[];
    velS=[];
    ang=[];
    errAng=[];
    delta_t_X=[];
    delta_t_Y=[]; 
    seriesX=[];
    seriesY=[];
    seriesZ=[];
    seriesRX=[];    
   
    startFreq=0.025;
    freqStep=.005;
    mockFreq=0.06;
    iter=floor((.1-startFreq)/freqStep);
    % Mock Data to produce 8*sqrt(2)/2 km/s with both methods at 45 degrees
    if m==-1
        mockFreq=0.06;
        BRSY= awgn(sqrt(2)*4*10^-8*sin(2*pi*mockFreq*(1:8*4000)/sampf).*exp(-((1:8*4000)/sampf-8*2000/sampf).^2/(50)^2),180);        
        ETMYZ = awgn(3.2*10^-4*sin(2*pi*mockFreq*(1:8*4000)/sampf).*exp(-((1:8*4000)/sampf-8*2000/sampf).^2/(50)^2),180);
        ETMXZ = awgn(3.2*10^-4*sin(2*pi*mockFreq*(1:8*4000)/sampf).*exp(-((1:8*4000)/sampf-8*2000/sampf).^2/(50)^2),180);
        ITMYZ = awgn(3.2*10^-4*sin(2*pi*mockFreq*(1-6:8*4000-6)/sampf).*exp(-((1-6:8*4000-6)/sampf-8*2000/sampf).^2/(50)^2),180);   
    else
        % Data Reading        
        if (exist(fileName{m+1},'file'))
            myfile = load(fileName{m+1});
            mydata = myfile.rawdata8Hz1;
            rawBRSY= mydata(:,4);        
            rawETMXZ = mydata(:,3);
            rawETMYZ = mydata(:,2);
            rawITMYZ = mydata(:,1);       
        end
        Sttime =1000*8;
        Endtime=length(rawBRSY)-9;

        ETMXZ=1e-9 *rawETMXZ(Sttime:Endtime);
        ETMYZ=1e-9 *rawETMYZ(Sttime:Endtime);
        ITMYZ=1e-9 *rawITMYZ(Sttime:Endtime);
        BRSY=1e-9*rawBRSY(Sttime:Endtime);
        %% Uncomment to mock data using real seismometer signal. 
        % Must also change BRSYcal_out to BRSY in line 143 to avoid
        % filtering out BRS oscillation on a STS signal.
        % Expected velocity: 4km/s @ 90 deg
        
%         ETMXZ=1e-9 *rawETMYZ(Sttime:Endtime);
%         ETMYZ=1e-9 *rawETMYZ(Sttime+8:Endtime+8);
%         ITMYZ=1e-9 *rawETMYZ(Sttime:Endtime);
%         BRSY=1/4e3*1e-9*rawETMYZ(Sttime+8:Endtime+8);

    end    
    %% Filters
    % Filter to remove low frequency junk caused by reponse inversion
    [bb,aa] = butter(3,[2*0.01/sampf 2*.200/sampf],'bandpass');
    
    Rot_time = transpose(1/sampf * (0:1:length(BRSY)-1));
    
    %STS response inversion filter
    STSInvertFilt = zpk(-2*pi*[pairQ(8.2e-3,0.7)],-2*pi*[0 0],1);
    STSInvertFilt = 1*STSInvertFilt/abs(freqresp(STSInvertFilt,2*pi*100));
    
    % Apply filterS
    T240cal_vel = lsim(STSInvertFilt,ETMYZ,Rot_time);    
    ETMYZ_out=filter(bb,aa,T240cal_vel);
    
    T240cal_vel = lsim(STSInvertFilt,ETMXZ,Rot_time);
    ETMXZ_out=filter(bb,aa,T240cal_vel);
    
    T240cal_vel = lsim(STSInvertFilt,ITMYZ,Rot_time);
    ITMYZ_out=filter(bb,aa,T240cal_vel);
    
    % BRSY repsonse inversion is done before data is written
    BRSY_out=filter(bb,aa,BRSY);

    time=(1:length(ETMYZ_out))/sampf;
 
    ETMXZ_out=ETMXZ_out(100*sampf:length(ETMXZ_out));
    ETMYZ_out=ETMYZ_out(100*sampf:length(ETMYZ_out));
    ITMYZ_out=ITMYZ_out(100*sampf:length(ITMYZ_out));
    BRSY_out=BRSY_out(100*sampf:length(BRSY_out));  

    figure(27)
    spectrogram(ETMYZ,1024*2,512*2,1024*2,8,'yaxis')
    ax = gca;
    ax.YScale = 'log';
    ylim([5e-3 0.5])
    title('Z Velocity')

    figure(28)
    spectrogram(BRSY,1024*2,512*2,1024*2,8,'yaxis')
    ax = gca;
    ax.YScale = 'log';
    ylim([5e-3 0.5])
    title('BRS')

    %Gives the filter time to ring down
    startTime=300*sampf;
    endTime=length(ETMXZ_out);
    ampThreshold=max(BRSY_out*1e9)/10;
      
    for a=0:iter
        %% Data Crunching  
        % Bandpass filtering to get data into frequency bins
        freq=(startFreq+a*freqStep);
        [bb,aa]= butter(3,[2*((a-1/2)*freqStep+startFreq)/sampf 2*((a+1/2)*freqStep+startFreq)/sampf],'bandpass');
        
        filtData=filter(bb,aa,ETMXZ_out-mean(ETMXZ_out));
        X=filtData(startTime:endTime);        
        filtData=filter(bb,aa,ETMYZ_out-mean(ETMYZ_out));
        Y=filtData(startTime:endTime);
        filtData=filter(bb,aa,ITMYZ_out-mean(ITMYZ_out));
        C=filtData(startTime:endTime);
        seriesZ=Y;
        filtData=filter(bb,aa,BRSY_out-mean(BRSY_out)); 
        seriesRX=filtData(startTime:endTime);
        
        % Fitting of the BRS and EndY seiesmometer signals to a*sin(w*t)+b*cos(w*t)
        FitData=zeros(length(seriesZ),1);

        tempZ=[];
        tempRX=[];
        r2z=[];
        r2rx=[];        
        fitLength=floor(1/(freq/sampf)/4);
        for j=1:floor(length(seriesRX)/fitLength)-2
           if (max(abs(1e9*seriesRX(j*fitLength:(j+2)*fitLength)))>ampThreshold)
               tim=(j*fitLength:(j+1)*fitLength)'./sampf;
               %Fitting the STS_Z. The signal is scaled to allow better fitting
               cut=1e9*seriesZ(j*fitLength:(j+1)*fitLength);
               g = fittype( @(a,b,cen_fr,x) a*sin(2*pi*cen_fr*x)+b*cos(2*pi*cen_fr*x), 'problem', 'cen_fr' );
               [myfit,st] = fit(tim,cut, g,'problem',freq,'StartPoint', [1, 1]);
               r2z=[r2z;st.rsquare];
               a1=abs(coeffvalues(myfit));
               %Fitting the BRS
               cut=1e12*seriesRX(j*fitLength:(j+1)*fitLength);
               g = fittype( @(a,b,cen_fr,x) a*sin(2*pi*cen_fr*x)+b*cos(2*pi*cen_fr*x), 'problem', 'cen_fr' );
               [myfit,st] = fit(tim,cut, g,'problem',freq,'StartPoint', [1, 1]);
               r2rx=[r2rx;st.rsquare];
               a2=abs(coeffvalues(myfit));
               FitData(j*fitLength:(j+1)*fitLength) = a2(1)*sin(2*pi*freq.*tim)+a2(2)*cos(2*pi*freq.*tim);
               % Extracting amplitude as a complex number
               if (abs(r2rx)>0) % Fit quality check
                   if (abs(r2z)>0)
                       if(abs(angle(a1(2)+i*a1(1))-angle((a2(2)+i*a2(1))))<(phaseThreshold*pi/180) && ...
                               abs((a2(2)+a2(1)*i)/1e3)>=noiseThreshold && abs((a2(2)+a2(1)*i)/1e3)>=ampThreshold) % Phase and Amplitude Threshold
                           tempZ=[tempZ abs(a1(2))+abs(a1(1))*i];
                           tempRX=[tempRX (abs(a2(2))+abs(a2(1))*i)/1e3];
                       else
                           tempZ=[tempZ NaN];
                           tempRX=[tempRX NaN];
                       end
                   end
               end
           else
               tempZ=[tempZ NaN];
               tempRX=[tempRX NaN];
           end
        end
        
        % Diagnostic plots
        figure(1)
        plot((1:length(tempZ))*fitLength/sampf,abs(tempZ)/1000,(1:length(tempRX))*fitLength/sampf,abs(tempRX))
        title('Fit Amplitudes')
        ylabel('Amplitude (nrad or um)')
        xlabel('Time (s)')
        legend('STS_Z','BRS')
        set(gca,'FontSize',12)
        grid on
        
        figure(2)
        line1=plot((1:length(FitData))/sampf,seriesRX*1e9,(1:length(FitData))/sampf,FitData/1e3,'--',(1:length(FitData))/sampf,...
            noiseThreshold*(1:length(FitData))./(1:length(FitData)),(1:length(FitData))/sampf,ampThreshold*(1:length(FitData))./(1:length(FitData)));
        ylabel('Angle (nrad)')
        xlabel('Time (s)')
        title('BRS Filtered Data and Fit')
        legend('Data','Fit','Noise Threshold','Amp Threshold')
        set(gca,'FontSize',12)
        grid on

        figure(3)            
        line2=plot((startTime:endTime)/sampf,X,(startTime:endTime)/sampf,Y+1e-6,(startTime:endTime)/sampf,C+2e-6);
        title('Filtered Seismometer Signals')
        legend('EndX','EndY','Corner')
        set(line2,'LineWidth',1.2)
        set(gca,'FontSize',12)
        grid on
        
        % Time of arrival delay calculations using the cross correlation
        % between vertical components of each seismometer        
        delta_t_X=[];
        delta_t_Y=[];
        len=floor(10/freq*sampf);
        peakLength=1/2; %Periods
        for j=1:floor(length(X)/len)-2
            if(max(abs(tempRX(floor(j*len/fitLength):floor((j+1)*len/fitLength))))>=noiseThreshold &&...
                    max(abs(tempRX(floor(j*len/fitLength):floor((j+1)*len/fitLength))))>=ampThreshold) % Amplitude noiseThreshold
                cornerCut=C(j*len:(j+1)*len);
                xCut=X(j*len:(j+1)*len);
                yCut=Y(j*len:(j+1)*len);
                %Cross corr calculation with windowed data to remove a discontinuity at zero
                [crossY,~] = xcorr(tukeywin(length(cornerCut),.1).*cornerCut,tukeywin(length(yCut),.1).*yCut);
                [crossX,lags]=xcorr(tukeywin(length(cornerCut),.1).*cornerCut,tukeywin(length(xCut),.1).*xCut);
                
                crossX=abs(crossX);
                crossY=abs(crossY);
                %Finds the points within +/- 3/2 wavelength around the peak value
                peak=crossX(floor(find(crossX==max(crossX))-sampf/freq/8*peakLength):floor(find(crossX==max(crossX))+sampf/freq/8*peakLength));
                peakLags=lags(floor(find(crossX==max(crossX))-sampf/freq/8*peakLength):floor(find(crossX==max(crossX))+sampf/freq/8*peakLength))';
                %Fits peak with gaussian
                f=fit(peakLags,10^12*peak,fittype('gauss1'));
                delta_t_X=[delta_t_X f.b1/sampf];

                peak=crossY(floor(find(crossY==max(crossY))-sampf/freq/8/peakLength):floor(find(crossY==max(crossY))+sampf/freq/8/peakLength));
                peakLags=lags(floor(find(crossY==max(crossY))-sampf/freq/8/peakLength):floor(find(crossY==max(crossY))+sampf/freq/8/peakLength))';
                f=fit(peakLags,10^12*peak,fittype('gauss1'));
                delta_t_Y=[delta_t_Y f.b1/sampf];   
            else
                delta_t_X=[delta_t_X NaN];
                delta_t_Y=[delta_t_Y NaN];
            end
        end
                
        %% Calculations
        % Velocity and angle of incidence calculations using array
        % angle of incidence=arctan((time delay along Y)/(time delay along X))
        % velocity=4 km/ sqrt((time delay along X)^2 + (time delay along Y)^2)
        tX=mean(delta_t_X(find(not(isnan(delta_t_X)))));
        tY=mean(delta_t_Y(find(not(isnan(delta_t_Y)))));        
        errTX=std(delta_t_X(find(not(isnan(delta_t_X)))));
        errTY=std(delta_t_Y(find(not(isnan(delta_t_Y)))));
        
        ang=[ang atan2(tY,tX)*180/pi];        
        errAng=[errAng sqrt((tY^2*errTX^2+tX^2*errTY^2)/(tX^2+tY^2))*180/pi];
        
        velA=[velA 4e3./sqrt(tX.^2+tY.^2)];
 
        % Average phase velocity calculations.
        % velocity=-(amplitude of vertical velocity)/(amplitude of rotation)* sin(angle of incidence)
        seriesVS=[];
        for p=1:min([length(tempZ) length(tempRX)])
            if(abs(tempRX(p))>=noiseThreshold && abs(tempRX(p))>=ampThreshold) % Amplitude noiseThreshold
                if(not(isnan(abs(abs(tempZ(p)./tempRX(p).*sin(ang(a+1)*pi/180)))))) % NaN check
                    seriesVS=[seriesVS abs(tempZ(p)./tempRX(p).*sin(ang(a+1)*pi/180))];
                else
                    seriesVS=[seriesVS NaN];
                end                
            end
        end
        velS=[velS; mean(seriesVS(find(not(isnan(seriesVS)))))];
        
    end
    cutAng=ang(intersect(find(not(isnan(errAng))),find(not(0==errAng)))); %Removes zero error points and NaNs
    cutErrAng=errAng(intersect(find(not(isnan(errAng))),find(not(0==errAng))));
    avgAng=[avgAng; sum(cutAng./(cutErrAng.^2))/sum(1./(cutErrAng.^2))]
    avgErrAng=[avgErrAng; sqrt(1/sum(1./(cutErrAng.^2)))]
 %% Single Event Plotting

    %Velocity dispersion plot
    fig1=figure(6);
    hold on
    line4=plot(((0:length(velS)-1)*freqStep+startFreq)*1000,velS'/1000,'Color',[0 0.4470 0.7410]);
    line5=plot(((0:length(velS)-1)*freqStep+startFreq)*1000,velA'/1000,'Color',[0.8500 0.3250 0.0980],'LineStyle','--');
    ylabel('Phase Velocity (km/s)')
    xlabel('Frequency (mHz)')
    set(gca,'XTick',5*(0:100))
    set(gca,'YTick',.5*(0:100))    
    legend('Single Station','Array')
    xlim([20,75]);
    ylim([0,8]);
    set(line4,'LineWidth',1.5)
    set(line5,'LineWidth',1.5)
    set(gca,'FontSize',16)
    grid on
    box on
    hold off           
            
end

RayWaveSpec