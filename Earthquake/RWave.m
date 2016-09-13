close all
singVel=[0];
singAng=[0];
% for j=0:2
j=0;
    sampf =8;
    startFreq=0.03;
    freqStep=.005;
    iter=floor((.11-startFreq)/freqStep);
    [errFreq,transXErr,transYErr,transZErr,tiltErr]=RWaveMeasErr;
    if j==0
        [ETMXZ_out, ITMYZ_out, ETMYX_out, ETMYY_out, ETMYZ_out, BRSY_out]= RWaveDataIn('GPS1143962787_6_9Earthquake.mat');
    end
    %     RWaveDataIn('GPS1149323500_Quite.mat');
    if j==1
        [ETMXZ_out, ITMYZ_out, ETMYX_out, ETMYY_out, ETMYZ_out, BRSY_out]= RWaveDataIn('GPS1144888165_7_8Earthquake.mat');
    end
    if j==2
        [ETMXZ_out, ITMYZ_out, ETMYX_out, ETMYY_out, ETMYZ_out, BRSY_out]= RWaveDataIn('GPS1149581095_5_2Earthquake.mat');
    end
    if j==3
        [ETMXZ_out, ITMYZ_out, ETMYX_out, ETMYY_out, ETMYZ_out, BRSY_out]=RWaveDataIn('GPS1149331885_6_2Earthquake.mat');
    end 
    bootVel=[];
    bootAng=[];
    startTime=1;
    endTime=length(ETMXZ_out);
    [bb,aa] = butter(4,[2*0.01/sampf, 2*.1/sampf]);

    % seed=randn(1,length(ETMYZ_out));
    threshSeries=filter(bb,aa,ETMYZ_out);
    threshSeries=threshSeries(300*sampf:length(threshSeries));    
    threshold=rms(threshSeries)/2;
    
    [vel, ang,bootVel,bootAng]=RWaveArray(ETMXZ_out,ETMYZ_out,ITMYZ_out,sampf,threshold,startFreq,freqStep,iter,startTime,endTime);
    
    [v,phi,el,k,bootV,bootPhi,bootEl,bootK]=...
    RWaveSingle(ETMYX_out,ETMYY_out,ETMYZ_out,BRSY_out,...
        'S',errFreq,transXErr,transYErr,transZErr,tiltErr,sampf,ang,threshold,startFreq,freqStep,iter,startTime,endTime);

%     cInd1=find((abs(v)+4*std(bootV')'>=abs(vel)'-4*std(bootVel')'));
%     cInd2=find((abs(v)-4*std(bootV')'<=abs(vel)'+4*std(bootVel')'));
%     cInd3=find(std(bootVel')<=1000);
%     cInd4=find(std(bootV')<=1000);
%     cInd5=find(std(bootVel')>=1);
%     cInd6=find(std(bootV')>=1);
%     cInd12=intersect(cInd1,cInd2);
%     cInd123=intersect(cInd12,cInd3);
%     cInd1234=intersect(cInd123,cInd4);
%     cInd12345=intersect(cInd1234,cInd5);
%     cInd=intersect(cInd12345,cInd6);
    cInd=find(abs(v)>=0);
    v=v(cInd);
    vel=vel(cInd);
    ang=ang(cInd);
    phi=phi(cInd);
    bootV=bootV(cInd,:);
    bootVel=bootVel(cInd,:);
    
    singVel=[singVel; mean(v)];
    singAng=[singAng; mean(ang)];

    % figure(1)
    % plot(Rot_time(startTime:length(BRSY_out)),seriesX,Rot_time(startTime:length(BRSY_out))...
    %     ,seriesY,Rot_time(startTime:length(BRSY_out)),seriesZ,Rot_time(startTime:length(BRSY_out)),seriesRX)
    % ylim([-1e-5 100e-5])
    % xlim([400 length(ETMYZ_out)/sampf]);
    % grid on

    figure(2)
    hold on
    % errorbar(((0:length(v)-1))*freqStep+startFreq,abs(v),-sigmaV,sigmaV)
    % errorbar(((0:length(vel)-1))*freqStep+startFreq,vel,-sigmaVel,sigmaVel)
    l=errorbar((cInd-1)*freqStep+startFreq,abs(v),-std(bootV'),std(bootV'));
    ll=errorbar((cInd-1)*freqStep+startFreq,vel,-std(bootVel'),std(bootVel'),'--');
    ylabel('Velocity (m/s)')
    xlabel('Frequency (Hz)')
    legend('Single Station Vanuatu', 'Array Vanuatu','Single Station Ecuador','Array Ecuador','Single Station California','Array California')
    grid on
    % xlim([.01 .1])
    ylim([0 1e4])
    set(gca,'FontSize',12)
    set(l,'LineWidth',1.2)
    set(ll,'LineWidth',1.2)
    hold off

    % figure(3)
    % hold on
    % % errorbar(((0:length(phi)-1))*freqStep+startFreq,phi,-sigmaPhi,sigmaPhi)
    % % errorbar(((0:length(ang)-1))*freqStep+startFreq,ang,-sigmaAng,sigmaAng)
    % errorbar(((0:length(phi)-1))*freqStep+startFreq,phi,-std(bootPhi'),std(bootPhi'))
    % errorbar(((0:length(ang)-1))*freqStep+startFreq,ang,-std(bootAng'),std(bootAng'))
    % ylabel('Angle of Incidence (degrees)')
    % xlabel('Frequency (Hz)')
    % legend('Single Station','Array','Single Station Bootstrapping','Array Bootstrapping')
    % grid on
    % % xlim([.01 .1])
    % ylim([-200 200])
    % hold off
    % figure(4)
    % plot(((0:100)+1/2)*0.005+0.02,mean(k'))
    % ylabel('Rayleigh Wave Wavenumber')
    % xlabel('Frequency (Hz)')
    % grid on
    % % xlim([.02 .0675])
    % % 
    % figure(5)
    % plot(((0:length(phi)-1))*freqStep+startFreq,el/pi)
    % ylabel('Rayleigh Ellipticity (pi rad)')
    % xlabel('Frequency (Hz)')
    % grid on
    % xlim([.02 .0675])
% end
%%
% figure(6)
% n=5;
% r=(0:n)'/n;
% singAng=sort(singAng);
% singVel=sort(singVel);
% singAng=[singAng; 360];
% singVel=[singVel; 0];
% X = r*cos(singAng'*pi/180-pi);
% Y = r*sin(singAng'*pi/180-pi);
% C=r./r*(abs(singVel.*cos(singAng*pi/180-pi)))';
% pcolor(X,Y,C)
% axis equal tight
% colorbar