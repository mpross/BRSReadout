
%6.7 Vanuatu 4/6/16, 7.1 Atlantic 8/29/16, 7.8 New Zealand 11/13/16
%7.9 Papa New Guinea 12/17/16, 7.2 New Caledonia 8/12/16, 
%7.9 Papa New Guinea 1/22/17, 6.5 Botswana 4/3/17
fileName={'GPS1143962787_6_9Earthquake.mat','GPS1156480217_Atlantic.mat',...
    'GPS1163070017_7k_EQ5.mat','GPS1166005817_10k_EQ7_PapNG.mat','GPS1155000017_7k_EQ8_NewC.mat','PNG2EQData.mat'...
    'GPS1177682913_6_3_Alaska.mat','GPS1177676532_6_2_Alaska.mat','GPS1171992542_6_9_Fiji.mat'};
%We changed data file formats so this accommodates different formats
newArray=[1, 2, 2, 2, 2, 3, 3, 3, 3];
sampf =8;

 for m=0:5
% for m=[0 2 3 4 5]
% for m=4
    Astop1 = 20;
    Apass  = .5;
    Astop2 = 20;
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
        %% Uncomment to mock data using real seismometer signal. 
        % Must also change BRSYcal_out to BRSY in line 143 to avoid
        % filtering out BRS oscillation on a STS signal.
        % Expected velocity: 4km/s @ 90 deg
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
    [bb,aa] = butter(3,[2*0.01/sampf 2*.100/sampf],'bandpass');

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
    
    startTime=1;
    endTime=length(BRSY_out);
    
    if(endTime>length(BRSY_out))
        endTime=length(BRSY_out)-1;
    end
    if(startTime<0)
        startTime=1;
    end
 
    ETMXZ_out=ETMXZ_out(startTime:endTime);
    ETMYZ_out=ETMYZ_out(startTime:endTime);
    ITMYZ_out=ITMYZ_out(startTime:endTime);
    BRSY_out=BRSY_out(startTime:endTime);
    
    Navg2 = 1;
    [TfIYEY,~]=tfe2(ITMYZ_out,ETMYZ_out,1/sampf,Navg2,3,@hann);
    [TfIYEX,~]=tfe2(ITMYZ_out,ETMXZ_out,1/sampf,Navg2,3,@hann);
    [TfrXEY,freqs2]=tfe2(BRSY_out,ETMYZ_out,1/sampf,Navg2,3,@hann);
    
    avgFreqs=(0:30)*0.005+0.0225;
    
    avgXCS=[];
    avgYCS=[];
    avgYRS=[];    
    
    for fI=1:floor(length(avgFreqs))-1
            aInd=intersect(find(freqs2>avgFreqs(fI)),find(freqs2<avgFreqs(fI+1)));
         if (abs(std(sqrt((angle(TfIYEX(aInd))/2/pi./freqs2(aInd)').^2+(angle(TfIYEY(aInd))/2/pi./freqs2(aInd)').^2)))<0.3)
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
    
    angX = atan2(tYS,tXS);
    
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

