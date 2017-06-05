
%6.7 Vanuatu 4/6/16, 7.2 New Caledonia 8/12/16, 7.1 Atlantic 8/29/16, 
%7.8 New Zealand 11/13/16, 7.9 Papua New Guinea 12/17/16, 
%7.9 Papua New Guinea 1/22/17
fileName={'GPS1143961145_6_7_Vanuatu.mat','GPS1155000413_7_2_NewCaledonia.mat',...
    'GPS1156480214_7_1_Atlantic.mat','GPS1163070193_7_8_NewZealand.mat','GPS1166007087_7_9_PapuaNewGuinea.mat'...
    ,'GPS1169094640_7_9_PapuaNewGuinea2.mat'};

sampf =8;

for m=0:5
    velAS=[];
    velSS=[];
    
%     mockFreq=0.06;
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
        Sttime =9;
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
    
    %Gives filter time to ring down and removes non earthquake times
    startTime=2300*sampf;
    endTime=length(BRSY_out);    
 
    ETMXZ_out=ETMXZ_out(startTime:endTime);
    ETMYZ_out=ETMYZ_out(startTime:endTime);
    ITMYZ_out=ITMYZ_out(startTime:endTime);
    BRSY_out=BRSY_out(startTime:endTime);
    
    %Transfer function
    Navg2 = 3;
    [TfIYEY,~]=tfe2(ITMYZ_out,ETMYZ_out,1/sampf,Navg2,1,@hann);
    [TfIYEX,~]=tfe2(ITMYZ_out,ETMXZ_out,1/sampf,Navg2,1,@hann);
    [TfrXEY,freqs2]=tfe2(BRSY_out,ETMYZ_out,1/sampf,Navg2,1,@hann);
    
    avgFreqs=(0:30)*0.005+0.0225;
    
    avgXCS=[];
    avgYCS=[];
    avgYRS=[];    
    
    %Averaging 5 mHz wide bins
    for fI=1:floor(length(avgFreqs))-1
        %Finds points within bin
         aInd=intersect(find(freqs2>avgFreqs(fI)),find(freqs2<avgFreqs(fI+1)));
        %Removes points with standard deviation of delay time > 0.3 s
         if (abs(std(sqrt((angle(TfIYEX(aInd))/2/pi./freqs2(aInd)').^2+(angle(TfIYEY(aInd))/2/pi./freqs2(aInd)').^2)))<0.2)
            avgXCS=[avgXCS mean(angle(TfIYEX(aInd)))];
            avgYCS=[avgYCS mean(angle(TfIYEY(aInd)))];
            avgYRS=[avgYRS mean(abs(TfrXEY(aInd)))];

         else
            avgXCS=[avgXCS nan];
            avgYCS=[avgYCS nan];
            avgYRS=[avgYRS nan];
         end
            
    end
    avgFreqs=avgFreqs(2:end)-.0025;
    tXS=avgXCS/2/pi./avgFreqs;
    tYS=avgYCS/2/pi./avgFreqs;
    %Angle calculation
    angX = atan2(tYS,tXS);
    %Velocity calculations
    velSS = abs(avgYRS.*sin(angX));    
    velAS=4e3./sqrt(tXS.^2+tYS.^2);
    
    figure(51);
    hold on
    line4=plot(avgFreqs*1000,velSS/1000,'Color',[0 0.4470 0.7410]);
    line5=plot(avgFreqs*1000,velAS/1000,'Color',[0.8500 0.3250 0.0980],'LineStyle','--');
    ylabel(' Phase Velocity (km/s)')
    xlabel('Frequency (mHz)')
    legend('Single Station','Array')
    xlim([20,75]);
    ylim([0,8]);
    set(line4,'LineWidth',1.5)
    set(line5,'LineWidth',1.5)
    set(gca,'FontSize',16)
    set(gca,'XTick',5*(0:100))
    set(gca,'YTick',.5*(0:100))
    grid on
    box on
 end

