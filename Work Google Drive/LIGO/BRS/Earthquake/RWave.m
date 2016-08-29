close all
readData=true;
sampf = 8; % sampling frequency in Hz
quad='E';
if quad=='E'
    sgnX=1;
    sgnY=1;
    sgnRX=1;
    sgnZ=1;
end
if quad=='S'
    sgnX=-1;
    sgnY=1;
    sgnRX=1;
    sgnZ=1;
end
if quad=='W'
    sgnX=-1;
    sgnY=-1;
    sgnRX=1;
    sgnZ=1;
end
if quad=='N'
    sgnX=1;
    sgnY=-1;
    sgnRX=1;
    sgnZ=1;
end
% % 
% if exist('GPS1149581095_5_2Earthquake.mat','file') && readData==true
%     myfile = load('GPS1149581095_5_2Earthquake.mat');
if exist('GPS1143962787_6_9Earthquake.mat','file') && readData==true
    myfile = load('GPS1143962787_6_9Earthquake.mat');
% if exist('GPS1149331885_6_2Earthquake.mat','file') && readData==true
%     myfile = load('GPS1149331885_6_2Earthquake.mat');
% if exist('GPS1144888165_7_8Earthquake.mat','file') && readData==true
%     myfile = load('GPS1144888165_7_8Earthquake.mat');
    mydata = myfile.mydata;
    rawETMXX= mydata(:,1);
    rawETMXY = mydata(:,2);
    rawETMXZ = mydata(:,3);
    rawETMYX = mydata(:,4);
    rawETMYY = mydata(:,5);
    rawETMYZ = mydata(:,6);
    rawITMYX = mydata(:,7);
    rawITMYY = mydata(:,8);
    rawITMYZ = mydata(:,9);
    rawBRSY= mydata(:,10);
    rawBRSX=mydata(:,11);
%     rawWINDX=mydata(:,12);
%     rawWINDY=mydata(:,13);
end   
Sttime =01*sampf;
Endtime=length(rawBRSY)-0*sampf;
localg = 9.8;
% BRSscale=.85;
BRSscale=1;

ETMXX=1e-9 *rawETMXX(Sttime:Endtime);
ETMXY=1e-9 *rawETMXY(Sttime:Endtime);
ETMXZ=1e-9 *rawETMXZ(Sttime:Endtime);
ETMYX=1e-9 *rawETMYX(Sttime:Endtime);
% ETMYX=seed*cos(-130/180*pi);
ETMYY=1e-9 *rawETMYY(Sttime:Endtime);
% ETMYY=seed*sin(-130/180*pi);
ETMYZ=1e-9 *rawETMYZ(Sttime:Endtime);
ITMYX=1e-9 *rawITMYX(Sttime:Endtime);
ITMYY=1e-9 *rawITMYY(Sttime:Endtime);
ITMYZ=1e-9 *rawITMYZ(Sttime:Endtime);
BRSY=1e-9* BRSscale*rawBRSY(Sttime:Endtime);
BRSX=2.13e-9* BRSscale*rawBRSX(Sttime:Endtime);
Navg =9;
seed=randn(1,length(ETMYZ));
%% BRS parameters
Mtotal = 4.5; % Total mass in kg
Ibar = 0.61; % Moment of Inertia in kg m^2
localg = 9.8; % local gravitation acceleration
doffsetY = .5e-6; % CoM offset from pivot in m. d is plus if CoM below pivot
doffsetX = 30e-6; 
resfreq = 2*pi*7.7e-3;
Qbar = 3000;
Rot_time = transpose(1/sampf * (0:1:length(BRSY)-1));

%% torque computation
[bb,aa] = butter(4,[2*0.02/sampf 2*0.100/sampf],'bandpass');
% [bb,aa] = butter(4,[2*0.005/sampf, 2*.5/sampf]);
% %T240 response inversion filter
T240InvertFilt = zpk(-2*pi*[pairQ(8.2e-3,0.7)],-2*pi*[0 0],1);
T240InvertFilt = 1*T240InvertFilt/abs(freqresp(T240InvertFilt,2*pi*100));
% 
% %BRS response inversion filter

BRSXInvertFilt = zpk(-2*pi*[pairQ(8.9e-3,4000)],-2*pi*[0 0],1);
BRSXInvertFilt = 1*BRSXInvertFilt/abs(freqresp(BRSXInvertFilt,2*pi*100));

BRSYInvertFilt = zpk(-2*pi*[pairQ(7.7e-3,3000)],-2*pi*[0 0],1);
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
T240cal_vel = lsim(T240InvertFilt,ETMXZ,Rot_time);
T240cal_disp = lsim(IntFilt,T240cal_vel,Rot_time);
ETMXZ_out=filter(bb,aa,T240cal_disp);

% 
T240cal_accY = lsim(DiffFilt,ETMYY,Rot_time);
%
T240cal_vel = lsim(T240InvertFilt,ETMYZ,Rot_time);
T240cal_disp = lsim(IntFilt,T240cal_vel,Rot_time);
ETMYZ_out=filter(bb,aa,T240cal_disp);

% 
T240cal_vel = lsim(T240InvertFilt,ETMYX,Rot_time);
% T240cal_vel = lsim(T240InvertFilt,seed*cos(130/180*pi),Rot_time);
T240cal_disp = lsim(IntFilt,T240cal_vel,Rot_time);
ETMYX_out=filter(bb,aa,T240cal_disp);

T240cal_vel = lsim(T240InvertFilt,ETMYY,Rot_time);
% T240cal_vel = lsim(T240InvertFilt,seed*sin(130/180*pi),Rot_time);
T240cal_disp = lsim(IntFilt,T240cal_vel,Rot_time);
ETMYY_out=filter(bb,aa,T240cal_disp);

T240cal_vel = lsim(T240InvertFilt,ETMXZ,Rot_time);
T240cal_disp = lsim(IntFilt,T240cal_vel,Rot_time);
ETMXZ_out=filter(bb,aa,T240cal_disp);

T240cal_vel = lsim(T240InvertFilt,ITMYZ,Rot_time);
T240cal_disp = lsim(IntFilt,T240cal_vel,Rot_time);
ITMYZ_out=filter(bb,aa,T240cal_disp);

BRSYcal_out = lsim(BRSYInvertFilt,BRSY, Rot_time);
BRSY_out=filter(bb,aa,BRSYcal_out);


%% Single Station
seriesX=[];
seriesY=[];
seriesZ=[];
seriesRX=[];
ampX=[];
ampY=[];
ampZ=[];
ampRX=[];
ampTim=[];
phi=[];
k=[];
el=[];
v=[];
lPhi=[];
sigmaPhi=[];
sigmaV=[];
sigmaX=3e-8;
sigmaY=3e-8;
sigmaZ=3e-8;
% sigmaRX=8e-10;
sigmaRX=2e-10;
Astop1 = 10;
Apass  = 1;
Astop2 = 10;
freqStep1=0.001;
startFreq1=0.03;
startTime=01*sampf;
% endTime=800*sampf;
endTime=length(ETMYY);
% startTime=2000*sampf;
% endTime=2600*sampf;
% sgnZ=sign(ETMYZ_out(find(abs(ETMYZ_out-mean(ETMYZ_out))==max(abs(ETMYZ_out-mean(ETMYZ_out)))))-mean(ETMYZ_out));
% sgnX=sgnZ*sign(ETMYX_out(find(abs(ETMYZ_out-mean(ETMYZ_out))==max(abs(ETMYZ_out-mean(ETMYZ_out)))))-mean(ETMYX_out));
% sgnY=sgnZ*sign(ETMYY_out(find(abs(ETMYZ_out-mean(ETMYZ_out))==max(abs(ETMYZ_out-mean(ETMYZ_out)))))-mean(ETMYY_out));
% sgnRX=sgnZ*sign(BRSY_out(find(abs(ETMYZ_out-mean(ETMYZ_out))==max(abs(ETMYZ_out-mean(ETMYZ_out)))))-mean(BRSY_out));
% sgnZ=1;
for i=0:100
    freq1=(startFreq1+i*freqStep1);
%     [bb,aa] = butter(2,[(startFreq+i*freqStep)/sampf, (startFreq+(i+1)*freqStep)/sampf]);
    d = designfilt('bandpassiir', ...
      'StopbandFrequency1',freq1*.8,'PassbandFrequency1', freq1*.9, ...
      'PassbandFrequency2',freq1*1.1,'StopbandFrequency2', freq1*1.2, ...
      'StopbandAttenuation1',Astop1,'PassbandRipple', Apass, ...
      'StopbandAttenuation2',Astop2, ...
      'DesignMethod','butter','SampleRate',sampf);

    filtData=filter(d,ETMYX_out-mean(ETMYX_out));  
    seriesX=[seriesX filtData(startTime:endTime)];%+10^-5*i];
    filtData=filter(d,ETMYY_out-mean(ETMYY_out)); 
    seriesY=[seriesY filtData(startTime:endTime)];%+10^-5*i+10^-6];
    filtData=filter(d,ETMYZ_out-mean(ETMYZ_out)); 
    seriesZ=[seriesZ filtData(startTime:endTime)];%+10^-5*i+2*10^-6];
    filtData=filter(d,(BRSY_out-mean(BRSY_out))); 
    seriesRX=[seriesRX filtData(startTime:endTime)];%+10^-5*i+3*10^-6];  
       
    tempX=[];
    tempY=[];
    tempZ=[];
    tempRX=[];
%     for j=1:floor(length(seriesX(:,i+1))*(freq1/sampf))*2-2
%        cut=seriesX(floor(j/(freq1/sampf)/2):floor((j+1)/(freq1/sampf)/2),i+1);
%        tempX=[tempX max(abs(cut))];%-min(cut)];
%        cut=seriesY(floor(j/(freq1/sampf)/2):floor((j+1)/(freq1/sampf)/2),i+1);
%        tempY=[tempY max(abs(cut))];%-min(cut)];
%        cut=seriesZ(floor(j/(freq1/sampf)/2):floor((j+1)/(freq1/sampf)/2),i+1);
%        tempZ=[tempZ max(abs(cut))];%-min(cut)];
%        cut=seriesRX(floor(j/(freq1/sampf)/2):floor((j+1)/(freq1/sampf)/2),i+1);
%        tempRX=[tempRX max(abs(cut))];%-min(cut)];
%     end        
    for j=1:floor(length(seriesX(:,i+1))*(freq1/sampf))-2
       cut=seriesX(floor(j/(freq1/sampf)):floor((j+1)/(freq1/sampf)),i+1);
       tempX=[tempX sgnX*(max(cut)-min(cut))];
       cut=seriesY(floor(j/(freq1/sampf)):floor((j+1)/(freq1/sampf)),i+1);
       tempY=[tempY sgnY*(max(cut)-min(cut))];
       cut=seriesZ(floor(j/(freq1/sampf)):floor((j+1)/(freq1/sampf)),i+1);
       tempZ=[tempZ sgnZ*(max(cut)-min(cut))];
       cut=seriesRX(floor(j/(freq1/sampf)):floor((j+1)/(freq1/sampf)),i+1);
       tempRX=[tempRX sgnRX*(max(cut)-min(cut))];
    end    
%     for j=2:length((startTime:endTime))-1
%         if sign(seriesX(j,i+1)-seriesX(j-1,i+1))*sign(seriesX(j+1,i+1)-seriesX(j,i+1))<0
%             tempX=[tempX abs(seriesX(j,i+1))];
%         end
%         if sign(seriesY(j,i+1)-seriesY(j-1,i+1))*sign(seriesY(j+1,i+1)-seriesY(j,i+1))<0
%             tempY=[tempY abs(seriesY(j,i+1))];
%         end
%         if sign(seriesZ(j,i+1)-seriesZ(j-1,i+1))*sign(seriesZ(j+1,i+1)-seriesZ(j,i+1))<0
%             tempZ=[tempZ abs(seriesZ(j,i+1))];
%         end
%          if sign(seriesRX(j,i+1)-seriesRX(j-1,i+1))*sign(seriesRX(j+1,i+1)-seriesRX(j,i+1))<0
%             tempRX=[tempRX abs(seriesRX(j,i+1))];
%         end
%     end
%     ampTim=[0:55]/freq1;
%     ampX=[ampX; tempX(1:56)];%+10^-5*i];
%     ampY=[ampY; tempY(1:56)];%+10^-5*i+10^-6];
%     ampZ=[ampZ; tempZ(1:56)];%+10^-5*i+2*10^-6];
%     ampRX=[ampRX; tempRX(1:56)];%+10^-5*i+3*10^-6];
%     
%     tempX=smooth(tempX);
%     tempY=smooth(tempY);
%     tempZ=smooth(tempZ);
%     tempRX=smooth(tempRX);
    
    N=0;
    sumX=0;
    sumY=0;
    sumZ=0;
    sumRX=0;
    avgPhi=0;
    avgK=0;
    avgEl=0;
    avgV=0;
    sumSigmaPhi=0;
    sumSigmaV=0;
    threshold=2e-6;
    for l=1:min([length(tempX) length(tempY) length(tempZ) length(tempRX)])
        if(abs(tempZ(l))>=threshold)
%         if l>=20
            avgPhi=avgPhi+atan2(tempY(l),tempX(l))*180/pi;
            avgK=avgK+tempRX(l)./tempZ(l)./sin(atan2(tempY(l),tempX(l)));
            avgV=avgV+2*pi*freq1.*tempZ(l)./(tempRX(l)).*sin(atan2(tempY(l),tempX(l)));%(ang(i+1)*pi/180+180);%
%             avgV=avgV+2*pi*freq1.*tempY(l)./(tempRX(l)).*tempZ(l)./tempY(l).*sin(atan2(tempY(l),tempX(l)));
            avgEl=avgEl+acot(tempZ(l)/tempY(l).*sin(atan2(tempY(l),tempX(l))));
%             sumX=sumX+tempX(l);
%             sumY=sumY+tempY(l);
%             sumZ=sumZ+tempZ(l);
%             sumRX=sumRX+tempRX(l);
            sumSigmaPhi=sumSigmaPhi+(tempY(l)^2*sigmaX^2+tempX(l)^2*sigmaY^2)/(tempX(l)^2+tempY(l)^2)^2;
            dVdZ=2*pi*freq1./(tempRX(l)).*sin(atan2(tempY(l),tempX(l)));
            dVdRX=-2*pi*freq1.*tempZ(l)./(tempRX(l)^2).*sin(atan2(tempY(l),tempX(l)));
            dVdY=2*pi*freq1.*tempZ(l)./(tempRX(l)).*cos(atan2(tempY(l),tempX(l))).*(tempX(l)/(tempX(l)^2+tempY(l)^2));
            dVdX=-2*pi*freq1.*tempZ(l)./(tempRX(l)).*cos(atan2(tempY(l),tempX(l))).*(tempY(l)/(tempX(l)^2+tempY(l)^2));
            sumSigmaV=sumSigmaV+dVdX^2*sigmaX^2+dVdY^2*sigmaY^2+dVdZ^2*sigmaZ^2+dVdRX^2*sigmaRX^2;
            N=N+1;
        end
    end   
    avgX=sumX/N;
    avgY=sumY/N;
    avgZ=sumZ/N;
    avgRX=sumRX/N;
    
    avgPhi=avgPhi/N;
    avgK=avgK/N;
    avgV=avgV/N;
    avgEl=avgEl/N;
% 
%     phi=[phi; atan2(avgY,avgX)*180/pi];
%     k=[k; avgRX/avgZ/sin(atan2(avgY,avgX))];
%     v=[v; 2*pi*freq1*avgZ/(avgRX)*sin(atan2(avgY,avgX))];
%     el=[el; acot(avgZ/avgY*sin(atan2(avgY,avgX)))];
    phi=[phi; avgPhi];
    k=[k; avgK];
    v=[v; avgV];
    el=[el; avgEl];
    sigmaPhi=[sigmaPhi; sqrt(sumSigmaPhi)/N*180/pi];
    sigmaV=[sigmaV; sqrt(sumSigmaV)/N];
end

% 
% v=v(find(abs(v-mean(v))<=3*std(abs(v))));
% 
% phi=phi(find(abs(phi-mean(phi))<=3*std(abs(phi))));
% 
% el=el(find(abs(el-mean(el))<=3*std(abs(el))));



%% Array

vel=[];
ang=[];
angTim=[];
velTim=[];
lagsX=[];
lagsY=[];
startFreq2=startFreq1;
freqStep2=freqStep1;
sigmaTX=0.0231;
sigmaTY=0.0325;
sigmaAng=[];
sigmaVel=[];
for i=0:length(v)
    freq2=(startFreq2+i*freqStep2);
%     [bb,aa] = butter(2,[(0.032+i*0.005)/sampf, (0.032+(i+1)*0.005)/sampf]);
    d = designfilt('bandpassiir', ...
      'StopbandFrequency1',freq2*.8,'PassbandFrequency1', freq2*.9, ...
      'PassbandFrequency2',freq2*1.1,'StopbandFrequency2', freq2*1.2, ...
      'StopbandAttenuation1',Astop1,'PassbandRipple', Apass, ...
      'StopbandAttenuation2',Astop2, ...
      'DesignMethod','butter','SampleRate',sampf);

    filtData=filter(d,ETMXZ_out);
    X=filtData(startTime:endTime);
    filtData=filter(d,ETMYZ_out);
    Y=filtData(startTime:endTime);
    filtData=filter(d,ITMYZ_out);
    C=filtData(startTime:endTime);
%     injAng=150*pi/180;
%     X=filter(d,seed(ceil(100*cos(injAng))+150:length(ETMXZ_out)+ceil(100*cos(injAng))+150))';
%     Y=filter(d,seed(ceil(100*sin(injAng))+150:length(ETMXZ_out)+ceil(100*sin(injAng))+150))';
%     C=filter(d,seed(1+150:length(ETMXZ_out)+1+150))';
%     injVel=7000;
%     X=filter(d,seed(4e3/(sqrt(2)*injVel)*10*sampf+150:length(ETMXZ_out)+4e3/(sqrt(2)*injVel)*10*sampf+150))';
%     Y=filter(d,seed(4e3/(sqrt(2)*injVel)*10*sampf+150:length(ETMXZ_out)+4e3/(sqrt(2)*injVel)*10*sampf+150))';
%     C=filter(d,seed(1+150:length(ETMXZ_out)+1+150))';

%     X=ETMXZ_out(startTime*sampf:length(ETMXZ_out));
%     Y=ETMYZ_out(startTime*sampf:length(ETMYZ_out));
%     C=ITMYZ_out(startTime*sampf:length(ITMYZ_out));
    if max(C>=threshold)
        [crossY,~] = xcorr(C,Y);
        [crossX,lags]=xcorr(C,X);
        crossX=abs(crossX);
        crossY=abs(crossY);
        peak=crossX(floor(.9995*find(crossX==max(crossX))):floor(1.0005*find(crossX==max(crossX))));
        peakLags=lags(floor(.9995*find(crossX==max(crossX))):floor(1.0005*find(crossX==max(crossX))))';
        [fit,s]=polyfit(peakLags,peak,2);  
        delta_t_X=-fit(2)/(2*fit(1))/sampf;        
%         sigmaTX=std(peak-(fit(1)*peakLags.^2+fit(2)*peakLags+fit(3)));
%         errFit= sqrt(diag(inv(s.R)*inv(s.R'))./s.normr.^2./s.df);
%         sigmaTX=sqrt((-1/(2*fit(1))/sampf)^2*errFit(2)^2+(fit(2)/(2*fit(1)^2)/sampf)^2*errFit(1)^2);

        peak=crossY(floor(.9995*find(crossY==max(crossY))):floor(1.0005*find(crossY==max(crossY))));
        peakLags=lags(floor(.9995*find(crossY==max(crossY))):floor(1.0005*find(crossY==max(crossY))))';
        [fit,s]=polyfit(peakLags,peak,2);  
        delta_t_Y=-fit(2)/(2*fit(1))/sampf;       
%         sigmaTY=std(peak-(fit(1)*peakLags.^2+fit(2)*peakLags+fit(3)));
%         errFit= sqrt(diag(inv(s.R)*inv(s.R'))./s.normr.^2./s.df);
%         sigmaTY=sqrt((-1/(2*fit(1))/sampf)^2*errFit(2)^2+(fit(2)/(2*fit(1)^2)/sampf)^2*errFit(1)^2);
    else
        delta_t_X=nan;
        delta_t_Y=nan;
    end
%     delta_t_X=(lags(find(crossX==max(crossX))))/sampf;
%     delta_t_Y=(lags(find(crossY==max(crossY))))/sampf;
%     delta_t_X=(sum(crossX.*lags)/sum(crossX))/sampf;
%     delta_t_Y=(sum(crossY.*lags)/sum(crossY))/sampf;
      
    lagsX=[lagsX delta_t_X];
    lagsY=[lagsY delta_t_Y];
    
    ang=[ang atan2(delta_t_Y,delta_t_X)*180/pi];
    vel=[vel 4e3/sqrt((delta_t_X)^2+(delta_t_Y)^2)];
    sigmaAng=[sigmaAng sqrt(delta_t_X^2*sigmaTY^2+delta_t_Y^2*sigmaTX^2)/(delta_t_X^2+delta_t_Y^2)*180/pi];
    sigmaVel=[sigmaVel 2*4e3*sqrt(delta_t_X^2*sigmaTX^2+delta_t_Y^2*sigmaTY^2)/(delta_t_X^2+delta_t_Y^2)^2];
end

% vel=vel(find(abs(vel-mean(vel))<=abs(3*std(vel))));
% 
% ang=ang(find(abs(ang-mean(ang))<=abs(3*std(ang))));

% Navg3 = 1;
% [Ty,~] = tfe2(ETMYZ_out,ITMYZ_out ,1/sampf, Navg3);
% [Tx,F] = tfe2(ETMXZ_out,ITMYZ_out,1/sampf, Navg3);
% 
% F=F(find(F>=.01 & F<=.1));
% Tx=Tx(find(F>=.01 & F<=.1));
% Ty=Ty(find(F>=.01 & F<=.1));
% 
% delta_t_X=angle(Tx)/pi./F';
% delta_t_Y=angle(Ty)/pi./F';
% 
% ang=atan(delta_t_Y./delta_t_X)*180/pi;
% vel=4e3./sqrt((delta_t_X).^2+(delta_t_Y).^2);
% 
% vel=vel(find(isinf(vel)==false));
% 
% F=F(find(abs(vel-mean(vel))<=abs(3*std(vel))));
% ang=ang(find(abs(vel-mean(vel))<=abs(3*std(vel))));
% vel=vel(find(abs(vel-mean(vel))<=abs(3*std(vel))));
% figure(1)
% plot(Rot_time(startTime:length(BRSY_out)),seriesX,Rot_time(startTime:length(BRSY_out))...
%     ,seriesY,Rot_time(startTime:length(BRSY_out)),seriesZ,Rot_time(startTime:length(BRSY_out)),seriesRX)
% ylim([-1e-5 100e-5])
% xlim([400 length(ETMYZ_out)/sampf]);
% grid on

figure(2)
hold on
errorbar(((0:length(v)-1))*freqStep1+startFreq1,abs(v),-sigmaV,sigmaV)
errorbar(((0:length(vel)-1))*freqStep2+startFreq2,vel,-sigmaVel,sigmaVel)
ylabel('Velocity (m/s)')
xlabel('Frequency (Hz)')
legend('Single Station','Array')
grid on
% xlim([.01 .1])
ylim([0 1e4])
hold off

figure(3)
hold on
errorbar(((0:length(phi)-1))*freqStep1+startFreq1,phi,-sigmaPhi,sigmaPhi)
errorbar(((0:length(ang)-1))*freqStep2+startFreq2,ang,-sigmaAng,sigmaAng)
ylabel('Angle of Incidence (degrees)')
xlabel('Frequency (Hz)')
legend('Single Station','Array')
grid on
% xlim([.01 .1])
ylim([-200 200])
hold off
% figure(4)
% plot(((0:100)+1/2)*0.005+0.02,mean(k'))
% ylabel('Rayleigh Wave Wavenumber')
% xlabel('Frequency (Hz)')
% grid on
% % xlim([.02 .0675])
% % 
% figure(5)
% plot(((0:length(phi)-1))*freqStep1+startFreq1,el/pi)
% ylabel('Rayleigh Ellipticity (pi rad)')
% xlabel('Frequency (Hz)')
% grid on
% xlim([.02 .0675])
