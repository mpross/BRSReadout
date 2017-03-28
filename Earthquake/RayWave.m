close all

avgVel=containers.Map('ValueType','any','KeyType','double');
avgV=containers.Map('ValueType','any','KeyType','double');
avgVelErr=containers.Map('ValueType','any','KeyType','double');
avgVErr=containers.Map('ValueType','any','KeyType','double');
avgVelP=containers.Map('ValueType','any','KeyType','double');
avgVP=containers.Map('ValueType','any','KeyType','double');
avgVelErrP=containers.Map('ValueType','any','KeyType','double');
avgVErrP=containers.Map('ValueType','any','KeyType','double');
fileName={'GPS1143962787_6_9Earthquake.mat','GPS1144888165_7_8Earthquake.mat','GPS1149581095_5_2Earthquake.mat',...
    'GPS1149331885_6_2Earthquake.mat','GPS1156069817_Burma.mat','GPS1156480217_Atlantic.mat','GPS1156782617_NewZealand.mat',...
    'GPS1157149817_Russia.mat','GPS1163070017_7k_EQ5.mat','GPS1163797217_7k_EQ6.mat','GPS1166005817_10k_EQ7_PapNG.mat','GPS1155000017_7k_EQ8_NewC.mat'};
% fileName={'GPS1166005817_10k_EQ7_PapNG_fast.mat','GPS1166005817_10k_EQ5_NZ_fast.mat','GPS1166005817_10k_EQ8_NewCal_fast.mat'};
newArray=[false, false, false, false, true, true, true, true, true, true, true, true,true];
sampf =8;
noiseThreshold=0.5; % approx 5*BRS noise
phaseThreshold=20; %20
plIndex=0;
ampIter=0;

for m=[0 5 8 10 11]
% for m=-1
% for m=8
% for m=[0 1 2]
    ampIter=ampIter+1;
    Astop1 = 20;
    Apass  = .5;
    Astop2 =20;
    velA=[];
    ang=[];
    angTim=[];
    velTim=[];
    lagsX=[];
    lagsY=[];      
    delta_t_X=[];
    delta_t_Y=[];
    errVelA=[];
    bootAng=[];    
    errVelS=[];
    seriesX=[];
    seriesY=[];
    seriesZ=[];
    seriesRX=[];    
    r2zavg=[];
    r2rxavg=[];
    ampX=[];
    ampY=[];
    ampZ=[];
    ampRX=[];
    ampTim=[];
    velS=[];
    plIndex=plIndex+1;
   
    startFreq=0.025;
%     startFreq=0.07;
    freqStep=.005;
    iter=floor((.1-startFreq)/freqStep);
%     iter=floor((.07-startFreq)/freqStep);
%     iter=1;
    % Mock Data to produce 8*sqrt(2)/2 km/s with both methods at 45 degrees
    if m==-1
        mockFreq=0.06;
        BRSY= awgn(sqrt(2)*4*10^-8*sin(2*pi*mockFreq*(1:8*4000)/sampf).*exp(-((1:8*4000)/sampf-8*2000/sampf).^2/(200)^2),180);        
        ETMYZ = awgn(3.2*10^-4*sin(2*pi*mockFreq*(1:8*4000)/sampf).*exp(-((1:8*4000)/sampf-8*2000/sampf).^2/(200)^2),180);
        ETMXZ = awgn(3.2*10^-4*sin(2*pi*mockFreq*(1:8*4000)/sampf).*exp(-((1:8*4000)/sampf-8*2000/sampf).^2/(200)^2),180);
        ITMYZ = awgn(3.2*10^-4*sin(2*pi*mockFreq*(1-6:8*4000-6)/sampf).*exp(-((1-6:8*4000-6)/sampf-8*2000/sampf).^2/(200)^2),180);   
    else
        % Data Reading
        if (exist(fileName{m+1},'file')&& newArray(m+1)==false)
            myfile = load(fileName{m+1});
            mydata = myfile.mydata;
            rawETMXZ = mydata(:,3);
            rawETMYZ = mydata(:,6);
            rawITMYZ = mydata(:,9);
            rawBRSY= mydata(:,10);
        end   
        if (exist(fileName{m+1},'file')&& newArray(m+1)==true)
            myfile = load(fileName{m+1});
            mydata = myfile.rawdata8Hz1;
            rawBRSY= mydata(:,4);        
            rawETMXZ = mydata(:,2);
            rawETMYZ = mydata(:,1);
            rawITMYZ = mydata(:,3);       
        end
%         if (exist(fileName{m+1},'file'))
%                 myfile = load(fileName{m+1});
%                 rawBRSY= myfile.fastBRSY.objs.data.x;        
%                 rawETMXZ = decimate(myfile.fastETMXZ.objs.data.x,2);
%                 rawETMYZ = decimate(myfile.fastETMYZ.objs.data.x,2);
%                 rawITMYZ = decimate(myfile.fastITMYZ.objs.data.x,2); 
%         end
        Sttime =01;
        Endtime=length(rawBRSY);

        ETMXZ=1e-9 *rawETMXZ(Sttime:Endtime);
        ETMYZ=1e-9 *rawETMYZ(Sttime:Endtime);
        ITMYZ=1e-9 *rawITMYZ(Sttime:Endtime);
        BRSY=1e-9*rawBRSY(Sttime:Endtime);
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
    [bb,aa] = butter(4,[2*0.02/sampf 2*0.100/sampf],'bandpass');
%     [bb,aa] = butter(2,[2*0.01/sampf, 2*1/sampf]);

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
    T240cal_disp = lsim(IntFilt,T240cal_vel,Rot_time);
%     ETMYZ_out=filter(bb,aa,T240cal_disp);
    ETMYZ_out=filter(bb,aa,filter(bb,aa,T240cal_vel));
    
    T240cal_vel = lsim(STSInvertFilt,ETMXZ,Rot_time);
    T240cal_disp = lsim(IntFilt,T240cal_vel,Rot_time);
%     ETMXZ_out=filter(bb,aa,T240cal_disp);
    ETMXZ_out=filter(bb,aa,filter(bb,aa,T240cal_vel));

    T240cal_vel = lsim(STSInvertFilt,ITMYZ,Rot_time);
    T240cal_disp = lsim(IntFilt,T240cal_vel,Rot_time);
%     ITMYZ_out=filter(bb,aa,T240cal_disp);
    ITMYZ_out=filter(bb,aa,filter(bb,aa,T240cal_vel));

    BRSYcal_out = lsim(BRSYInvertFilt,BRSY, Rot_time);
    BRSY_out=filter(bb,aa,filter(bb,aa,BRSYcal_out));

    tim=(1:length(ETMYZ_out))/sampf;

    ETMXZ_out=ETMXZ_out(300*sampf:length(ETMXZ_out));
    ETMYZ_out=ETMYZ_out(300*sampf:length(ETMYZ_out));
    ITMYZ_out=ITMYZ_out(300*sampf:length(ITMYZ_out));
    BRSY_out=BRSY_out(300*sampf:length(BRSY_out));
    errVelA=[];
    bootAng=[];
    startTime=1000*sampf;
    endTime=length(ETMXZ_out);
    ampThreshold=max(BRSY_out*1e9)/5;
      
    for a=0:iter
%     for a=0
%     for a=7
        %% Data Crunching  
        % Bandpass filtering to get data into frequency bins
        freq=(startFreq+a*freqStep);
        d = designfilt('bandpassiir', ...
          'StopbandFrequency1',(a-2)*freqStep+startFreq,'PassbandFrequency1', (a-1)*freqStep+startFreq, ...
          'PassbandFrequency2',(a+1)*freqStep+startFreq,'StopbandFrequency2', (a+2)*freqStep+startFreq, ...
          'StopbandAttenuation1',Astop1,'PassbandRipple', Apass, ...
          'StopbandAttenuation2',Astop2, ...
          'DesignMethod','butter','SampleRate',sampf);
        filtData=filter(d,ETMXZ_out-mean(ETMXZ_out));
        X=filtData(startTime:endTime);        
        filtData=filter(d,ETMYZ_out-mean(ETMYZ_out));
        Y=filtData(startTime:endTime);
        filtData=filter(d,ITMYZ_out-mean(ITMYZ_out));
        C=filtData(startTime:endTime);
        filtData=filter(d,ETMYZ_out-mean(ETMYZ_out));
        seriesZ=Y;
        filtData=filter(d,BRSY_out-mean(BRSY_out)); 
        seriesRX=filtData(startTime:endTime);
        
        % Fitting of the BRS and EndY seiesmometer signals to a*sin(w*t)+b*cos(w*t)
        FitData=zeros(length(seriesZ),1);
        tempX=[];
        tempY=[];
        tempZ=[];
        tempRX=[];
        r2z=[];
        r2rx=[];        
        fitLength=floor(1/(freq/sampf)/4);
        for j=1:floor(length(seriesRX)/fitLength)-2
           tim=(j*fitLength:(j+1)*fitLength)'./sampf;
           cut=1e9*seriesZ(j*fitLength:(j+1)*fitLength);
           g = fittype( @(a,b,cen_fr,x) a*sin(2*pi*cen_fr*x)+b*cos(2*pi*cen_fr*x), 'problem', 'cen_fr' );
           [myfit,st] = fit(tim,cut, g,'problem',freq,'StartPoint', [1, 1]);
           r2z=[r2z;st.rsquare];
           a1=coeffvalues(myfit);
           cut=1e12*seriesRX(j*fitLength:(j+1)*fitLength);
           g = fittype( @(a,b,cen_fr,x) a*sin(2*pi*cen_fr*x)+b*cos(2*pi*cen_fr*x), 'problem', 'cen_fr' );
           [myfit,st] = fit(tim,cut, g,'problem',freq,'StartPoint', [1, 1]);
           r2rx=[r2rx;st.rsquare];
           a2=coeffvalues(myfit);           
           FitData(j*fitLength:(j+1)*fitLength) = a2(1)*sin(2*pi*freq.*tim)+a2(2)*cos(2*pi*freq.*tim);
           % Extracting amplitude as a complex number
           if (abs(r2rx)>0) 
               if (abs(r2z)>0)
                   if(abs(angle(a1(2)+i*a1(1))-angle((a2(2)+i*a2(1))))<(phaseThreshold*pi/180) && ...
                           abs((a2(2)+a2(1)*i)/1e3)>=noiseThreshold && abs((a2(2)+a2(1)*i)/1e3)>=ampThreshold) % Phase and Amplitude Threshold
                       tempZ=[tempZ a1(2)+a1(1)*i];
                       tempRX=[tempRX (a2(2)+a2(1)*i)/1e3];
                   else
                       tempZ=[tempZ NaN];
                       tempRX=[tempRX NaN];
                   end
               end
           end
        end
        figure(6)
        histogram(((angle(tempRX)-angle(tempZ))*180/pi),20)
        legend('phase')
        % Diagnostic plots
        figure(1)
        plot((1:length(tempZ)),abs(tempZ)/1000,(1:length(tempRX)),abs(tempRX))
        
        figure(2)
        line1=plot((1:length(FitData))/sampf,seriesRX*1e9,(1:length(FitData))/sampf,FitData/1e3,'--',(1:length(FitData))/sampf,...
            noiseThreshold*(1:length(FitData))./(1:length(FitData)),(1:length(FitData))/sampf,ampThreshold*(1:length(FitData))./(1:length(FitData)));
        ylabel('Angle (nrad)')
        xlabel('Time (s)')
        legend('Data','Fit','Noise Threshold','Amp Threshold')
        set(line1,'LineWidth',1.2)
        set(gca,'FontSize',12)
        grid on

        figure(3)            
        line2=plot((startTime:endTime)/sampf,X,(startTime:endTime)/sampf,Y+1e-6,(startTime:endTime)/sampf,C+2e-6);
        legend('X','Y','C')
        set(line2,'LineWidth',1.2)
        set(gca,'FontSize',12)
        grid on
        
        % Time of arrival delay calculations using the cross correlation
        % between vertical components of each seismometer        
        delta_t_X=[];
        delta_t_Y=[];
%         len=100*sampf;
        len=floor(2/freq*sampf);
        for j=1:floor(length(X)/len)-2
            if(max(abs(tempRX(floor(j*len/fitLength):floor((j+1)*len/fitLength))))>=noiseThreshold &&...
                    max(abs(tempRX(floor(j*len/fitLength):floor((j+1)*len/fitLength))))>=ampThreshold) % Amplitude noiseThreshold
                cornerCut=C(j*len:(j+1)*len);
                xCut=X(j*len:(j+1)*len);
                yCut=Y(j*len:(j+1)*len);
%                 [crossY,~] = xcorr(hamming(length(cornerCut)).*cornerCut,hamming(length(yCut)).*yCut);
%                 [crossX,lags]=xcorr(hamming(length(cornerCut)).*cornerCut,hamming(length(xCut)).*xCut);
                [crossY,~] = xcorr(tukeywin(length(cornerCut),0.1).*cornerCut,tukeywin(length(yCut),0.1).*yCut);
                [crossX,lags]=xcorr(tukeywin(length(cornerCut),0.1).*cornerCut,tukeywin(length(xCut),0.1).*xCut);

                crossX=abs(crossX);
                crossY=abs(crossY);

                peak=crossX(floor(find(crossX==max(crossX))-sampf/freq/8/1.5):floor(find(crossX==max(crossX))+sampf/freq/8/1.5));
                peakLags=lags(floor(find(crossX==max(crossX))-sampf/freq/8/1.5):floor(find(crossX==max(crossX))+sampf/freq/8/1.5))';
%                 [myfit,s]=polyfit(peakLags,peak,2);
%                 delta_t_X=[delta_t_X -myfit(2)/(2*myfit(1))/sampf];
                
                f=fit(peakLags,10^12*peak,fittype('gauss1'));
                delta_t_X=[delta_t_X f.b1/sampf];
%                 delta_t_X=[delta_t_X lags(find(crossX==max(crossX)))/sampf];

                figure(4)
%                 line3=plot(lags/sampf,crossX);
                line3=plot(lags/sampf,crossX,peakLags/sampf,f.a1*exp(-((peakLags-f.b1)/f.c1).^2)/10^12);
%                 line3=plot(lags/sampf,crossX,peakLags/sampf,(myfit(1)*peakLags.^2+myfit(2)*peakLags+myfit(3)));
                ylabel('Cross Cor')
                xlabel('Lag (s)')
                legend('Data','Fit')
                set(line3,'LineWidth',1.2)
                set(gca,'FontSize',12)
                grid on

                peak=crossY(floor(find(crossY==max(crossY))-sampf/freq/8/1.5):floor(find(crossY==max(crossY))+sampf/freq/8/1.5));
                peakLags=lags(floor(find(crossY==max(crossY))-sampf/freq/8/1.5):floor(find(crossY==max(crossY))+sampf/freq/8/1.5))';
%                 [myfit,s]=polyfit(peakLags,peak,2);
%                 delta_t_Y=[delta_t_Y -myfit(2)/(2*myfit(1))/sampf];
                f=fit(peakLags,10^12*peak,fittype('gauss1'));
                delta_t_Y=[delta_t_Y f.b1/sampf];
%                 delta_t_Y=[delta_t_Y lags(find(crossY==max(crossY)))/sampf];
            end
        end
        
%         %% Bootstrapping
%         % Runs velocity calculations on a large number of randomly sampled
%         % sets to estimate errors.
%         tempBAng=[];
%         tempBVelA=[];
%         tempBVelS=[];
%         for k=0:1e4
%             array=[delta_t_X(find(not(isnan(delta_t_X)))); delta_t_Y(find(not(isnan(delta_t_Y))))]; % NaN check
%             bootArray=bootstrapData(array');
%             if length(bootArray)>0  % Empty check
%                 bootTX=bootArray(:,1)';
%                 bootTY=bootArray(:,2)';   
%             else
%                 bootTX=nan;
%                 bootTY=nan;
%             end
%             tempBAng=[tempBAng; atan2(mean(bootTY),mean(bootTX))*180/pi];      
%             tempBVelA=[tempBVelA; 4e3./sqrt((mean(bootTX)).^2+(mean(bootTY)).^2)];
%             
%             N=0;
%             sumZ=0;
%             sumRX=0;
%             avgVS=0;
%             array=[tempZ(find(not(isnan(tempZ)))); tempRX(find(not(isnan(tempRX))))];  % NaN check
%             bootArray=bootstrapData(array');
%             btempZ=bootArray(:,1)';
%             btempRX=bootArray(:,2)';
%             for p=1:min([length(btempZ) length(btempRX)])
%                 if(abs(btempRX(p))>=noiseThreshold) % Amplitude noiseThreshold
%                    if not(isnan(abs(btempZ(p)./btempRX(p).*sin(atan2(mean(bootTY),mean(bootTX))*180/pi))))
%                        avgVS=avgVS+abs(btempZ(p)./btempRX(p).*sin(atan2(mean(bootTY),mean(bootTX))*180/pi));                    
%                        N=N+1;
%                    end
%                 end                
%             end   
%             avgVS=avgVS/N;
%             if (not(isnan(avgVS)))
%                 tempBVelS=[tempBVelS; avgVS];
%             end
%         end
% %         errVelS=[errVelS; std(abs(tempZ(p)./tempRX(p).*sin(ang(a+1)*pi/180)))]  %std(tempBVelS')
% %         errVelA=[errVelA; std(4e3./sqrt((delta_t_X(find(not(isnan(delta_t_X))))).^2+(delta_t_Y(find(not(isnan(delta_t_Y)))))).^2)]  %std(tempBVelA')
%         bootAng=[bootAng; tempBAng'];
        
        %% Calculations
        % Velocity and angle of incidence calculations using array
        % angle of incidence=arctan((time delay along Y)/(time delay along X))
        % velocity=4 km/ sqrt((time delay along X)^2 + (time delay along Y)^2)
        tempAng=atan2((delta_t_Y(find(not(isnan(delta_t_Y))))),(delta_t_X(find(not(isnan(delta_t_X))))))*180/pi;
        tempVel=4e3./sqrt((delta_t_X(find(not(isnan(delta_t_X))))).^2+(delta_t_Y(find(not(isnan(delta_t_Y))))).^2);
        ang=[ang mean(tempAng(find(not(isinf(tempVel)))))];
        velA=[velA mean(abs(tempVel(find(not(isinf(tempVel))))))];        
%         errVelS=[errVelS; std(abs(tempZ(find(not(isnan(tempZ))))./tempRX(find(not(isnan(tempRX)))).*sin(mean(tempAng)*pi/180)))];  %std(tempBVelS')
        errVelA=[errVelA; std(abs(tempVel(find(not(isinf(tempVel))))))];  %std(tempBVelA')
        
%         tX=delta_t_X(find(not(isnan(delta_t_X))));
%         tY=delta_t_Y(find(not(isnan(delta_t_Y))));
%         errX=std(tX);
%         errY=std(tY);
%         mtX=mean(tX);
%         mtY=mean(tY);
%         c=cov(tX,tY);        
%         dX=-4e3*mtX/(mtX^2+mtY^2)^(3/2);
%         dY=-4e3*mtY/(mtX^2+mtY^2)^(3/2);
%         errVelA=[errVelA; sqrt(dX^2*errX^2+dY^2+errY^2)];
        
        % Average phase velocity calculations.
        % velocity=-(amplitude of vertical velocity)/(amplitude of rotation)* sin(angle of incidence)
        N=0;
        avgPhi=0;
        avgK=0;
        avgEl=0;
        pltV=[];
        seriesVS=[];
        for p=1:min([length(tempZ) length(tempRX)])
            if(abs(tempRX(p))>=noiseThreshold && abs(tempRX(p))>=ampThreshold) % Amplitude noiseThreshold
                if(not(isnan(abs(abs(tempZ(p)./tempRX(p).*sin(ang(a+1)*pi/180)))))) % NaN check
                    seriesVS=[seriesVS abs(tempZ(p)./tempRX(p).*sin(ang(a+1)*pi/180))];
                end
            end
        end
        seriesVS
        errVelS=[errVelS; std(seriesVS)];
        velS=[velS; mean(seriesVS)];
        r2zavg=[mean(r2z);r2zavg];
        r2rxavg=[mean(r2rx);r2rxavg];
        % Diagnostic velocity plot       
        figure(5)
        line4=plot((1:length(seriesVS))*fitLength,seriesVS,(1:length(delta_t_X))*50*sampf/2,4e3./sqrt((delta_t_X).^2+(delta_t_Y).^2));
        set(line4,'LineWidth',1.2)
        set(gca,'FontSize',12)
        legend('Single Station','Array')
        ylim([0 1e4])
        grid on     
        
        
        figure(16)
        histogram(delta_t_X,20)
        legend('deltaT_X')
        
        figure(17)
        histogram(delta_t_Y,20)
        legend('deltaT_Y')
        
        figure(18)
        histogram(abs(tempZ(find(not(isnan(tempZ))))./tempRX(find(not(isnan(tempRX))))*sin(mean(tempAng)*pi/180)),20)
        legend('velS')
        
        figure(19)
        histogram(4e3./sqrt((delta_t_X(find(not(isnan(delta_t_X))))).^2+(delta_t_Y(find(not(isnan(delta_t_Y))))).^2),20)
        legend('velA')
        
    end   

 %% Single Event Plotting


    cInd1=find(not(isnan(velA)));
    cInd2=find(not(isnan(velS)));
    cInd=intersect(cInd1,cInd2);
%     velS=velS(cInd);
%     velA=velA(cInd);
%     errVelS=errVelS(cInd);
%     errVelA=errVelA(cInd);
     
    fig1=figure(6+plIndex)
    hold on
    line4=errorbar(((0:length(velS)-1)*freqStep+startFreq),velS'/1000,-errVelS'/1000,errVelS'/1000);
    line5=errorbar(((0:length(velA)-1)*freqStep+startFreq),velA/1000,-errVelA'/1000,errVelA'/1000,'--');
    ylabel('Velocity (km/s)')
    xlabel('Frequency (Hz)')
    if plIndex==1
        legend('Single Station-6.7 Vanuatu', 'Array-6.7 Vanuatu')
    end
    if plIndex==2
        legend('Single Station-7.1 Pacific', 'Array-7.1 Pacific')
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
    xlim([0.01,0.11]);
    ylim([0,8]);
    set(line4,'LineWidth',1.5)
    set(line5,'LineWidth',1.5)
    set(gca,'FontSize',12)
    grid on
    box on
    hold off
    
    print(fig1,'-dpng',strcat('figure',num2str(plIndex+1),'.png'))

    comV=intersect(cInd1,cell2mat(avgV.keys));
    comVel=intersect(cInd2,cell2mat(avgVel.keys));
    for k=1:length(comV)
        if(not(isnan(velS(comV(k))))&&not(isinf(velS(comV(k))))&&not(isnan(errVelS(comV(k)))))
            avgV(comV(k))=[avgV(comV(k)) velS(comV(k))];
            avgVErr(comV(k))=[avgVErr(comV(k)) errVelS(comV(k))'];
        end
    end
    for k=1:length(comVel)
        if(not(isnan(velA(comVel(k))))&& not(isinf(velA(comVel(k))))&&not(isnan(errVelA(comVel(k)))))
            avgVel(comVel(k))=[avgVel(comVel(k)) velA(comVel(k))];
            avgVelErr(comVel(k))=[avgVelErr(comVel(k)) errVelA(comVel(k))'];
        end
    end
    newV=setxor(cInd1,intersect(cInd1,cell2mat(avgV.keys)));
    newVel=setxor(cInd2,intersect(cInd2,cell2mat(avgVel.keys)));
    for k=1:length(newV)
        if(not(isnan(velS(newV(k))))&&not(isinf(velS(newV(k))))&&not(isnan(errVelS(newV(k)))))
            avgV(newV(k))=velS(newV(k))
            avgVErr(newV(k))=errVelS(newV(k));
        end
    end
    for k=1:length(newVel)
        if(not(isnan(velA(newVel(k))))&&not(isinf(velA(newVel(k))))&&not(isnan(errVelA(newVel(k)))))
            avgVel(newVel(k))=velA(newVel(k))
            avgVelErr(newVel(k))=errVelA(newVel(k));
        end
    end
        
end
%%
cInd1=cell2mat(avgV.keys);
for k=1:length(avgV.keys)
    N=length(avgV(cInd1(k)));
    avgVP(cInd1(k))=sum(abs(avgV(cInd1(k))))/N;
    avgVErrP(cInd1(k))=sqrt(sum(avgVErr(cInd1(k)).^2))/N;
end
cInd2=cell2mat(avgVel.keys);
for k=1:length(avgVel.keys)
    N=length(avgVel(cInd2(k)));
    avgVelP(cInd2(k))=sum(abs(avgVel(cInd2(k))))/N;
    avgVelErrP(cInd2(k))=sqrt(sum(avgVelErr(cInd2(k)).^2))/N;
end

fig2=figure(20)
hold on
l=errorbar(((cInd1-1)*freqStep+startFreq)+10^-4,cell2mat(avgVP.values)/1000,-cell2mat(avgVErrP.values)/1000,cell2mat(avgVErrP.values)/1000);
ll=errorbar(((cInd2-1)*freqStep+startFreq),cell2mat(avgVelP.values)/1000,-cell2mat(avgVelErrP.values)/1000,cell2mat(avgVelErrP.values)/1000,'--');
ylabel('Average Phase Velocity (km/s)')
xlabel('Frequency (Hz)')
legend('Single Station','Array')
set(l,'LineWidth',1.2)
set(ll,'LineWidth',1.2)
set(gca,'FontSize',12)
set(gca,'XTick',0.01*(0:100))
set(gca,'YTick',.5*(0:100))
xlim([0.01,0.11]);
ylim([0,8]);
grid on
box on
email=1;
if email==1
    print(fig2,'-dpng','figure1.png')
%     print(fig2,'-dpng','figure2.png')

    setpref('Internet','E_mail',strcat(user,'@gmail.com'));
    setpref('Internet','SMTP_Server','smtp.gmail.com');
    setpref('Internet','SMTP_Username',user);
    setpref('Internet','SMTP_Password',password);
    props = java.lang.System.getProperties;
    props.setProperty('mail.smtp.auth','true');
    props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
    props.setProperty('mail.smtp.socketFactory.port','465');
    sendmail('mpross2@uw.edu','Rayleigh Wave Analysis Complete!','',{'figure1.png','figure2.png','figure3.png','figure4.png','figure5.png','figure6.png'});
end
