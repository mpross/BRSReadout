% Containers for final averaging
avgVelS=containers.Map('ValueType','any','KeyType','double');
avgVS=containers.Map('ValueType','any','KeyType','double');
avgVelErrS=containers.Map('ValueType','any','KeyType','double');
avgVErrS=containers.Map('ValueType','any','KeyType','double');
avgVelPS=containers.Map('ValueType','any','KeyType','double');
avgVPS=containers.Map('ValueType','any','KeyType','double');
avgVelErrPS=containers.Map('ValueType','any','KeyType','double');
avgVErrPS=containers.Map('ValueType','any','KeyType','double');
%6.7 Vanuatu 4/6/16, 7.1 Atlantic 8/29/16, 7.8 New Zealand 11/13/16
%7.9 Papa New Guinea 12/17/16, 7.2 New Caledonia 8/12/16, 
%7.9 Papa New Guinea 1/22/17, 6.5 Botswana 4/3/17
fileName={'GPS1143962787_6_9Earthquake.mat','GPS1156480217_Atlantic.mat','GPS1163070017_7k_EQ5.mat','GPS1166005817_10k_EQ7_PapNG.mat','GPS1155000017_7k_EQ8_NewC.mat','PNG2EQData.mat'};
%We changed data file formats so this accommodates different formats
newArray=[1, 2, 2, 2, 2, 3];
sampf =8;

 for m=0:5
    Astop1 = 20;
    Apass  = .5;
    Astop2 = 20;
    velAS=[];
    velSS=[];
    errAS=[]; 
    errSS=[];
    
    startFreq=0.025;
    freqStep=.005;
%     mockFreq=0.06;
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
% 
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
    [bb,aa] = butter(3,[2*0.01/sampf 2*0.100/sampf],'bandpass');

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
%     T240cal_disp = lsim(IntFilt,T240cal_vel,Rot_time);
    ETMYZ_out=filter(bb,aa,ETMYZ);
    
    T240cal_vel = lsim(STSInvertFilt,ETMXZ,Rot_time);
%     T240cal_disp = lsim(IntFilt,T240cal_vel,Rot_time);
    ETMXZ_out=filter(bb,aa,ETMXZ);

    T240cal_vel = lsim(STSInvertFilt,ITMYZ,Rot_time);
%     T240cal_disp = lsim(IntFilt,T240cal_vel,Rot_time);
    ITMYZ_out=filter(bb,aa,ITMYZ);

    BRSYcal_out = lsim(BRSYInvertFilt,BRSY, Rot_time);
    BRSY_out=filter(bb,aa,BRSYcal_out);

    time=(1:length(ETMYZ_out))/sampf;
 
    ETMXZ_out=ETMXZ_out(100*sampf:length(ETMXZ_out));
    ETMYZ_out=ETMYZ_out(100*sampf:length(ETMYZ_out));
    ITMYZ_out=ITMYZ_out(100*sampf:length(ITMYZ_out));
    BRSY_out=BRSY_out(100*sampf:length(BRSY_out));
    
    Navg2 = 1;
    [TfIYEY,~]=tfe2(ITMYZ_out,ETMYZ_out,1/sampf,Navg2,1,@hann);
    [TfIYEX,~]=tfe2(ITMYZ_out,ETMXZ_out,1/sampf,Navg2,1,@hann);
    [TfrXEY,freqs2]=tfe2(BRSY_out,ETMYZ_out,1/sampf,Navg2,1,@hann);
    
    avgFreqs=[];
    
    avgXCS=[];
    avgYCS=[];
    avgYRS=[];
    
    errXCS=[];
    errYCS=[];
    errYRS=[];
    
    for fI=1:floor(length(freqs2)/20)-1
        if(freqs2(fI*20)<0.1)
            avgFreqs=[avgFreqs mean(freqs2(fI*20:(fI+1)*20))];

            avgXCS=[avgXCS mean(angle(TfIYEX(fI*20:(fI+1)*20)))];
            avgYCS=[avgYCS mean(angle(TfIYEY(fI*20:(fI+1)*20)))];
            avgYRS=[avgYRS mean(abs(TfrXEY(fI*20:(fI+1)*20)))];

            errXCS=[errXCS std(angle(TfIYEX(fI*20:(fI+1)*20)))];
            errYCS=[errYCS std(angle(TfIYEY(fI*20:(fI+1)*20)))];
            errYRS=[errYRS std(abs(TfrXEY(fI*20:(fI+1)*20)))];
        end
    end
    
    tXS=avgXCS/2/pi./avgFreqs;
    tYS=avgYCS/2/pi./avgFreqs;
    
    errTXS=errXCS/2/pi./avgFreqs;
    errTYS=errYCS/2/pi./avgFreqs;
    
    angX = atan2(tYS,tXS);
    errAngS=sqrt(tYS.^2.*errTXS.^2+tXS.^2.*errTYS.^2)./(tXS.^2+tYS.^2);        
    
    errAS=sqrt((4e3)^2*(errTXS.^2+errTYS.^2)./((tXS.^2+tYS.^2).^3));
    errSS=sqrt((errYRS.*sin(angX)).^2+(errAngS.*avgYRS.*sin(angX)).^2);
    
    velSS = abs(avgYRS.*sin(angX));    
    velAS=4e3./sqrt(tXS.^2+tYS.^2);
    for t=1:length(velSS)
         if (errAS(t)>2e3 || errSS(t)>2e3)
            velSS(t)=NaN;
            velAS(t)=NaN; 
            errAS(t)=NaN;
            errSS(t)=NaN;
        end
    end
    
    figure(51);
    hold on
    errorbar(avgFreqs+10^-4,velSS/1000,-errSS/1000,errSS/1000);
    errorbar(avgFreqs,velAS/1000,-errAS/1000,errAS/1000,'--');
    ylabel('Phase Velocity (km/s)')
    xlabel('Frequency (Hz)')
    legend('Single Station','Array')
    xlim([0.01,0.11]);
    ylim([0,8]);
    set(gca,'FontSize',12)
    set(gca,'XTick',0.01*(0:100))
    set(gca,'YTick',.5*(0:100))
    grid on
    box on
    
    cInd1=find(not(isnan(velAS)));
    cInd2=find(not(isnan(velSS)));
    cInd=intersect(cInd1,cInd2);
    
    %Creating a vector of different length vectors corresponding to the
    %velocities of all the earthquakes that share a bin
    comV=intersect(cInd1,cell2mat(avgVS.keys));
    comVel=intersect(cInd2,cell2mat(avgVelS.keys));
    
    for k=1:length(comV)
        if(not(isnan(velSS(comV(k))))&&not(isinf(velSS(comV(k))))&&not(isnan(errSS(comV(k))))&&not(0==errSS(comV(k))))
            avgVS(comV(k))=[avgVS(comV(k)) velSS(comV(k))];
            avgVErrS(comV(k))=[avgVErrS(comV(k)) errSS(comV(k))'];
        end
    end
    for k=1:length(comVel)
        if(not(isnan(velAS(comVel(k))))&& not(isinf(velAS(comVel(k))))&&not(isnan(errAS(comVel(k))))&&not(0==errAS(comVel(k))))
            avgVelS(comVel(k))=[avgVelS(comVel(k)) velAS(comVel(k))];
            avgVelErrS(comVel(k))=[avgVelErrS(comVel(k)) errAS(comVel(k))'];
        end
    end
    newV=setxor(cInd1,intersect(cInd1,cell2mat(avgVS.keys)));
    newVel=setxor(cInd2,intersect(cInd2,cell2mat(avgVelS.keys)));
    for k=1:length(newV)
        if(not(isnan(velSS(newV(k))))&&not(isinf(velSS(newV(k))))&&not(isnan(errSS(newV(k))))&&not(0==errSS(newV(k))))
            avgVS(newV(k))=velSS(newV(k));
            avgVErrS(newV(k))=errSS(newV(k));
        end
    end
    for k=1:length(newVel)
        if(not(isnan(velAS(newVel(k))))&&not(isinf(velAS(newVel(k))))&&not(isnan(errAS(newVel(k))))&&not(0==errAS(newVel(k))))
            avgVelS(newVel(k))=velAS(newVel(k));
            avgVelErrS(newVel(k))=errAS(newVel(k));
        end
    end
 end
 cInd1=cell2mat(avgVS.keys);
for k=1:length(avgVS.keys)
    N=length(avgVS(cInd1(k)));
%     avgVErrPS(cInd1(k))=sqrt(1/sum(1./avgVErr(cInd1(k)).^2));
%     avgVP(cInd1(k))=sum(abs(avgVS(cInd1(k))./(avgVErr(cInd1(k)).^2)))/sum(1/avgVErrPS(cInd1(k)).^2);
    avgVErrPS(cInd1(k))=std(avgVS(cInd1(k)));
%     avgVErrPS(cInd1(k))=sqrt(sum(avgVErr(cInd1(k)).^2))/length(avgVErr(cInd1(k)));
    avgVPS(cInd1(k))=mean(avgVS(cInd1(k)));
end
cInd2=cell2mat(avgVelS.keys);
for k=1:length(avgVelS.keys)
    N=length(avgVelS(cInd2(k)));
%     avgVelErrPS(cInd2(k))=sqrt(1/sum(1./avgVelErr(cInd2(k)).^2));
%     avgVelP(cInd2(k))=sum(abs(avgVelS(cInd2(k))./(avgVelErr(cInd2(k)).^2)))/sum(1/avgVelErrPS(cInd2(k)).^2);
    avgVelErrPS(cInd2(k))=std(avgVelS(cInd2(k)));
%     avgVelErrPS(cInd2(k))=sqrt(sum(avgVelErr(cInd2(k)).^2))/length(avgVelErr(cInd2(k)));
    avgVelPS(cInd2(k))=mean(avgVelS(cInd2(k)));
end
% Final dispersion plot for both methods
fig2=figure(121);
hold on
l=errorbar(((cInd1-1)*freqStep+startFreq)+10^-4,cell2mat(avgVPS.values)/1000,-cell2mat(avgVErrPS.values)/1000,cell2mat(avgVErrPS.values)/1000);
ll=errorbar(((cInd2-1)*freqStep+startFreq),cell2mat(avgVelPS.values)/1000,-cell2mat(avgVelErrPS.values)/1000,cell2mat(avgVelErrPS.values)/1000,'--');
ylabel('Average Phase Velocity (km/s)')
xlabel('Frequency (Hz)')
legend('Single Station','Array')
set(l,'LineWidth',1.2)
set(ll,'LineWidth',1.2)
set(gca,'FontSize',12)
set(gca,'XTick',0.005*(0:100))
set(gca,'YTick',.5*(0:100))
xlim([0.01,0.11]);
ylim([0,8]);
grid on
box on