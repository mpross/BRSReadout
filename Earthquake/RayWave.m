% close all
% Containers for final averaging
avgVel=containers.Map('ValueType','any','KeyType','double');
avgV=containers.Map('ValueType','any','KeyType','double');
avgVelErr=containers.Map('ValueType','any','KeyType','double');
avgVErr=containers.Map('ValueType','any','KeyType','double');
avgVelP=containers.Map('ValueType','any','KeyType','double');
avgVP=containers.Map('ValueType','any','KeyType','double');
avgVelErrP=containers.Map('ValueType','any','KeyType','double');
avgVErrP=containers.Map('ValueType','any','KeyType','double');
avgVelErrS=containers.Map('ValueType','any','KeyType','double');
avgVErrS=containers.Map('ValueType','any','KeyType','double');
%6.7 Vanuatu 4/6/16, 7.1 Atlantic 8/29/16, 7.8 New Zealand 11/13/16
%7.9 Papa New Guinea 12/17/16, 7.2 New Caledonia 8/12/16, 
%7.9 Papa New Guinea 1/22/17, 6.5 Botswana 4/3/17
fileName={'GPS1143962787_6_9Earthquake.mat','GPS1156480217_Atlantic.mat',...
    'GPS1163070017_7k_EQ5.mat','GPS1166005817_10k_EQ7_PapNG.mat','GPS1155000017_7k_EQ8_NewC.mat','PNG2EQData.mat'...
    'GPS1177682913_6_3_Alaska.mat','GPS1177676532_6_2_Alaska.mat','GPS1171992542_6_9_Fiji.mat'};
%We changed data file formats so this accommodates different formats
newArray=[1, 2, 2, 2, 2, 3, 3, 3, 3];
sampf =8;
noiseThreshold=0.5; % approx 5*BRS noise
phaseThreshold=20; % threshold for phase between seismometer signal and BRS signal
plIndex=0; %plotting index
avgAng=[];
avgErrAng=[];

 for m=0:5
    Astop1 = 20;
    Apass  = .5;
    Astop2 = 20;
    velA=[];
    velS=[];
    ang=[];     
    delta_t_X=[];
    delta_t_Y=[];
    errVelA=[];
    errAng=[];  
    errVelS=[];
    seriesX=[];
    seriesY=[];
    seriesZ=[];
    seriesRX=[];    
    r2zavg=[];
    r2rxavg=[];
    plIndex=plIndex+1;
   
    startFreq=0.025;
%     startFreq=0.06;
    freqStep=.005;
%     mockFreq=0.06;
    iter=floor((.1-startFreq)/freqStep);
%     iter=0;
    % Mock Data to produce 8*sqrt(2)/2 km/s with both methods at 45 degrees
    if m==-1
        mockFreq=0.06;
        BRSY= awgn(sqrt(2)*4*10^-8*sin(2*pi*mockFreq*(1:8*4000)/sampf).*exp(-((1:8*4000)/sampf-8*2000/sampf).^2/(50)^2),180);        
        ETMYZ = awgn(3.2*10^-4*sin(2*pi*mockFreq*(1:8*4000)/sampf).*exp(-((1:8*4000)/sampf-8*2000/sampf).^2/(50)^2),180);
        ETMXZ = awgn(3.2*10^-4*sin(2*pi*mockFreq*(1:8*4000)/sampf).*exp(-((1:8*4000)/sampf-8*2000/sampf).^2/(50)^2),180);
        ITMYZ = awgn(3.2*10^-4*sin(2*pi*mockFreq*(1-6:8*4000-6)/sampf).*exp(-((1-6:8*4000-6)/sampf-8*2000/sampf).^2/(50)^2),180);   
    else
        % Data Reading
        if (exist(fileName{m+1},'file')&& newArray(m+1)==1)
            myfile = load(fileName{m+1});
            mydata = myfile.mydata;
            rawETMXZ = mydata(:,3);
            rawETMYZ = mydata(:,6);
            rawITMYZ = mydata(:,9);
            rawBRSY= mydata(:,10);
        end   
        if (exist(fileName{m+1},'file')&& newArray(m+1)==2)
            myfile = load(fileName{m+1});
            mydata = myfile.rawdata8Hz1;
            rawBRSY= mydata(:,4);        
            rawETMXZ = mydata(:,2);
            rawETMYZ = mydata(:,1);
            rawITMYZ = mydata(:,3);       
        end
        if (exist(fileName{m+1},'file')&& newArray(m+1)==3)
            myfile = load(fileName{m+1});
            mydata = myfile.rawdata8Hz1;
            rawBRSY= mydata(:,4);        
            rawETMXZ = mydata(:,3);
            rawETMYZ = mydata(:,2);
            rawITMYZ = mydata(:,1);       
        end
        Sttime =9;
        Endtime=length(rawBRSY)-9;

        ETMXZ=1e-9 *rawETMXZ(Sttime:Endtime);
        ETMYZ=1e-9 *rawETMYZ(Sttime:Endtime);
        ITMYZ=1e-9 *rawITMYZ(Sttime:Endtime);
        BRSY=1e-9*rawBRSY(Sttime:Endtime);
%         ETMXZ=1e-9 *rawETMYZ(Sttime:Endtime);
%         ETMYZ=1e-9 *rawETMYZ(Sttime+8:Endtime+8);
%         ITMYZ=1e-9 *rawETMYZ(Sttime:Endtime);
%         BRSY=1/4e3*1e-9*rawETMYZ(Sttime+8:Endtime+8);
    end    
    
    localg = 9.8;
    BRSscale=1;  
    Navg =9;
    %% BRS parameters
    Mtotal = 4.5; % Total mass in kg
    Ibar = 0.61; % Moment of Inertia in kg m^2
    localg = 9.8; % local gravitation acceleration
    doffsetY = .5e-6; % CoM offset from pivot in m. d is plus if CoM below pivot
    doffsetX = 30e-6; 
    resfreq = 2*pi*7.7e-3;
    Qbar = 3000;
    Rot_time = transpose(1/sampf * (0:1:length(BRSY)-1));

    %% Response inversion
    [bb,aa] = butter(3,[2*0.01/sampf 2*1.00/sampf],'bandpass');

    % %STS response inversion filter
    STSInvertFilt = zpk(-2*pi*[pairQ(8.3e-3,0.7)],-2*pi*[0 0],1);
    STSInvertFilt = 1*STSInvertFilt/abs(freqresp(STSInvertFilt,2*pi*100));
    % 
    % %BRS response inversion filter
    BRSYInvertFilt = zpk(-2*pi*[pairQ(7.74e-3,3000)],-2*pi*[0 0],1);
    BRSYInvertFilt = 1*BRSYInvertFilt/abs(freqresp(BRSYInvertFilt,2*pi*100));
    % 
    % %Filters to differentiate and integrate
    DiffFilt = zpk(-2*pi*[0], -2*pi*2,1);
    DiffFilt = 1*DiffFilt/abs(freqresp(DiffFilt,2*pi*0.1592));
    % 
    IntFilt = zpk(-2*pi*[], -2*pi*5e-4,1);
    IntFilt = 1*IntFilt/abs(freqresp(IntFilt,2*pi*0.1592));
    % 
    % %Apply filters
    %
    T240cal_vel = lsim(STSInvertFilt,ETMYZ,Rot_time);
    ETMYZ_out=filter(bb,aa,ETMYZ);
    
    T240cal_vel = lsim(STSInvertFilt,ETMXZ,Rot_time);
    ETMXZ_out=filter(bb,aa,ETMXZ);

    T240cal_vel = lsim(STSInvertFilt,ITMYZ,Rot_time);
    ITMYZ_out=filter(bb,aa,ITMYZ);

    BRSYcal_out = lsim(BRSYInvertFilt,BRSY, Rot_time);
    BRSY_out=filter(bb,aa,BRSYcal_out);

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
    title('ETMYZ')

    figure(28)
    spectrogram(BRSY,1024*2,512*2,1024*2,8,'yaxis')
    ax = gca;
    ax.YScale = 'log';
    ylim([5e-3 0.5])
    title('BRSY')

    %Gives the filters time to ring down
    startTime=500*sampf;
    endTime=length(ETMXZ_out);
%     endTime=4000*sampf;
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
               a1=coeffvalues(myfit);
               %Fitting the BRS
               cut=1e12*seriesRX(j*fitLength:(j+1)*fitLength);
               g = fittype( @(a,b,cen_fr,x) a*sin(2*pi*cen_fr*x)+b*cos(2*pi*cen_fr*x), 'problem', 'cen_fr' );
               [myfit,st] = fit(tim,cut, g,'problem',freq,'StartPoint', [1, 1]);
               r2rx=[r2rx;st.rsquare];
               a2=coeffvalues(myfit);
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
        aTime=[];
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
%                 delta_t_X=[delta_t_X lags(find(crossX==max(crossX)))/sampf];

                peak=crossY(floor(find(crossY==max(crossY))-sampf/freq/8/peakLength):floor(find(crossY==max(crossY))+sampf/freq/8/peakLength));
                peakLags=lags(floor(find(crossY==max(crossY))-sampf/freq/8/peakLength):floor(find(crossY==max(crossY))+sampf/freq/8/peakLength))';
                f=fit(peakLags,10^12*peak,fittype('gauss1'));
                delta_t_Y=[delta_t_Y f.b1/sampf];   
%                 delta_t_Y=[delta_t_Y lags(find(crossY==max(crossY)))/sampf];
            else
                delta_t_X=[delta_t_X NaN];
                delta_t_Y=[delta_t_Y NaN];
            end
            aTime=[aTime j*len/sampf];
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
        velA=[velA 4e3./sqrt(tX^2+tY^2)];
        errVelA=[errVelA; sqrt((4e3)^2*(errTX^2+errTY^2)/(tX^2+tY^2)^3)];
 
        % Average phase velocity calculations.
        % velocity=-(amplitude of vertical velocity)/(amplitude of rotation)* sin(angle of incidence)
        N=0;
        avgPhi=0;
        avgK=0;
        avgEl=0;
        pltV=[];
        seriesVS=[];
        sTime=[];
        for p=1:min([length(tempZ) length(tempRX)])
            if(abs(tempRX(p))>=noiseThreshold && abs(tempRX(p))>=ampThreshold) % Amplitude noiseThreshold
                if(not(isnan(abs(abs(tempZ(p)./tempRX(p).*sin(ang(a+1)*pi/180)))))) % NaN check
                    seriesVS=[seriesVS abs(tempZ(p)./tempRX(p).*sin(ang(a+1)*pi/180))];
                else
                    seriesVS=[seriesVS NaN];
                end                
                sTime=[sTime p*fitLength/sampf];
            end
        end
        errVelS=[errVelS; std(seriesVS(find(not(isnan(seriesVS)))))];
        velS=[velS; mean(seriesVS(find(not(isnan(seriesVS)))))];
        r2zavg=[mean(r2z);r2zavg];
        r2rxavg=[mean(r2rx);r2rxavg];
        % Diagnostic velocity plot   
        figure(4)
        hold on
        line4=plot(sTime-2700+100,seriesVS,aTime-2700+100,4e3./sqrt((delta_t_X).^2+(delta_t_Y).^2)-500,'--');
        set(line4,'LineWidth',1.2)
        set(gca,'FontSize',12)
        title('Velocity Over Time (One freq bin)')
        legend('Single Station','Array')
        xlabel('Time (s)')
        ylabel('Phase Velocity (m/s)')
        ylim([0 1e4])
        grid on    
        
    end
    cutAng=ang(intersect(find(not(isnan(errAng))),find(not(0==errAng)))); %Removes zero error points and NaNs
    cutErrAng=errAng(intersect(find(not(isnan(errAng))),find(not(0==errAng))));
    avgAng=[avgAng; sum(cutAng./(cutErrAng.^2))/sum(1./(cutErrAng.^2))]
    avgErrAng=[avgErrAng; sqrt(1/sum(1./(cutErrAng.^2)))]
 %% Single Event Plotting

    %Another NaN check
    cInd1=find(not(isnan(velA)));
    cInd2=find(not(isnan(velS)));
    cInd=intersect(cInd1,cInd2);
    %Velocity dispersion plot for each event
    fig1=figure(6+plIndex);
    hold on
    line4=errorbar(((0:length(velS)-1)*freqStep+startFreq),velS'/1000,-errVelS'/1000,errVelS'/1000);
    line5=errorbar(((0:length(velA)-1)*freqStep+startFreq),velA/1000,-errVelA'/1000,errVelA'/1000,'--');
    ylabel('Velocity (km/s)')
    xlabel('Frequency (Hz)')
    if plIndex==1
        legend('Single Station-6.7 Vanuatu', 'Array-6.7 Vanuatu')
    end
    if plIndex==2
        legend('Single Station-7.1 Atlantic', 'Array-7.1 Atlantic')
    end
    if plIndex==3
        legend('Single Station-7.8 New Zealand', 'Array-7.8 New Zealand')
    end
    if plIndex==4
        legend('Single Station-7.9 Papa New Guinea','Array-7.9 Papa New Guinea')
    end
    if plIndex==5
        legend('Single Station-7.2 New Caledonia', 'Array-7.2 New Caledonia')
    end
    if plIndex==6
        legend('Single Station-7.9 Papa New Guinea 2nd','Array-7.9 Papa New Guinea 2nd')
    end
    if plIndex==7
        legend('Single Station-6.5 Botswana', 'Array-6.5 Botswana')
    end
    xlim([0.01,0.11]);
    ylim([0,8]);
    set(line4,'LineWidth',1.5)
    set(line5,'LineWidth',1.5)
    set(gca,'FontSize',12)
    grid on
    box on
    hold off    
        
    %Creating a vector of different length vectors corresponding to the
    %velocities of all the earthquakes that share a bin
    comV=intersect(cInd1,cell2mat(avgV.keys));
    comVel=intersect(cInd2,cell2mat(avgVel.keys));
    
    for k=1:length(comV)
        if(not(isnan(velS(comV(k))))&&not(isinf(velS(comV(k))))&&not(isnan(errVelS(comV(k))))&&not(0==errVelS(comV(k))))
            avgV(comV(k))=[avgV(comV(k)) velS(comV(k))];
            avgVErr(comV(k))=[avgVErr(comV(k)) errVelS(comV(k))'];
        end
    end
    for k=1:length(comVel)
        if(not(isnan(velA(comVel(k))))&& not(isinf(velA(comVel(k))))&&not(isnan(errVelA(comVel(k))))&&not(0==errVelA(comVel(k))))
            avgVel(comVel(k))=[avgVel(comVel(k)) velA(comVel(k))];
            avgVelErr(comVel(k))=[avgVelErr(comVel(k)) errVelA(comVel(k))'];
        end
    end
    newV=setxor(cInd1,intersect(cInd1,cell2mat(avgV.keys)));
    newVel=setxor(cInd2,intersect(cInd2,cell2mat(avgVel.keys)));
    for k=1:length(newV)
        if(not(isnan(velS(newV(k))))&&not(isinf(velS(newV(k))))&&not(isnan(errVelS(newV(k))))&&not(0==errVelS(newV(k))))
            avgV(newV(k))=velS(newV(k));
            avgVErr(newV(k))=errVelS(newV(k));
        end
    end
    for k=1:length(newVel)
        if(not(isnan(velA(newVel(k))))&&not(isinf(velA(newVel(k))))&&not(isnan(errVelA(newVel(k))))&&not(0==errVelA(newVel(k))))
            avgVel(newVel(k))=velA(newVel(k));
            avgVelErr(newVel(k))=errVelA(newVel(k));
        end
    end
        
end
%%
% Averaging down to one velocity per frequency bin
cInd1=cell2mat(avgV.keys);
for k=1:length(avgV.keys)
    N=length(avgV(cInd1(k)));
    avgVErrP(cInd1(k))=sqrt(1/sum(1./avgVErr(cInd1(k)).^2));
    avgVP(cInd1(k))=sum(abs(avgV(cInd1(k))./(avgVErr(cInd1(k)).^2)))/sum(1/avgVErrP(cInd1(k)).^2);
%     avgVErrP(cInd1(k))=std(avgV(cInd1(k)));
%     avgVErrP(cInd1(k))=sqrt(sum(avgVErr(cInd1(k)).^2))/length(avgVErr(cInd1(k)));
%     avgVP(cInd1(k))=mean(avgV(cInd1(k)));
end
cInd2=cell2mat(avgVel.keys);
for k=1:length(avgVel.keys)
    N=length(avgVel(cInd2(k)));
    avgVelErrP(cInd2(k))=sqrt(1/sum(1./avgVelErr(cInd2(k)).^2));
    avgVelP(cInd2(k))=sum(abs(avgVel(cInd2(k))./(avgVelErr(cInd2(k)).^2)))/sum(1/avgVelErrP(cInd2(k)).^2);
%     avgVelErrP(cInd2(k))=std(avgVel(cInd2(k)));
%     avgVelErrP(cInd2(k))=sqrt(sum(avgVelErr(cInd2(k)).^2))/length(avgVelErr(cInd2(k)));
%     avgVelP(cInd2(k))=mean(avgVel(cInd2(k)));
end

% Final dispersion plot for both methods
fig2=figure(21);
hold on
l=errorbar(((cInd1-1)*freqStep+startFreq)+10^-4,cell2mat(avgVP.values)/1000,-cell2mat(avgVErrP.values)/1000,cell2mat(avgVErrP.values)/1000);
ll=errorbar(((cInd2-1)*freqStep+startFreq),cell2mat(avgVelP.values)/1000,-cell2mat(avgVelErrP.values)/1000,cell2mat(avgVelErrP.values)/1000,'--');
ylabel('Average Phase Velocity (km/s)')
xlabel('Frequency (Hz)')
legend('Single Station','Array')
set(l,'LineWidth',1.2)
set(ll,'LineWidth',1.2)
set(gca,'FontSize',12)
set(gca,'XTick',0.005*(0:100))
set(gca,'YTick',.5*(0:100))
xlim([0.02,0.075]);
ylim([0,8]);
grid on
box on

RayWaveSpec