close all
singVel=[0];
singAng=[0];
% for j=0:2
j=5;
    sampf =8;
    startFreq=0.03;
    if j==2
        freqStep=.01;
    else
        freqStep=.005;
    end
    iter=floor((.12-startFreq)/freqStep);
    [errFreq,transXErr,transYErr,transZErr,tiltErr]=RWaveMeasErr;
    if j==0
        [ETMXZ_out, ITMYZ_out, ETMYX_out, ETMYY_out, ETMYZ_out, BRSY_out]= RWaveDataIn('GPS1143962787_6_9Earthquake.mat',false);
    end
%         [ETMXZ_out, ITMYZ_out, ETMYX_out, ETMYY_out, ETMYZ_out, BRSY_out]= RWaveDataIn('GPS1149323500_Quite.mat',new);
    if j==1
        [ETMXZ_out, ITMYZ_out, ETMYX_out, ETMYY_out, ETMYZ_out, BRSY_out]= RWaveDataIn('GPS1144888165_7_8Earthquake.mat',false);
    end
    if j==2
        [ETMXZ_out, ITMYZ_out, ETMYX_out, ETMYY_out, ETMYZ_out, BRSY_out]= RWaveDataIn('GPS1149581095_5_2Earthquake.mat',false);
    end
    if j==3
        [ETMXZ_out, ITMYZ_out, ETMYX_out, ETMYY_out, ETMYZ_out, BRSY_out]=RWaveDataIn('GPS1149331885_6_2Earthquake.mat',false);
    end 
    if j==4
        [ETMXZ_out, ITMYZ_out, ETMYX_out, ETMYY_out, ETMYZ_out, BRSY_out]=RWaveDataIn('GPS1156069817_Burma.mat',true);
    end 
    if j==5
        [ETMXZ_out, ITMYZ_out, ETMYX_out, ETMYY_out, ETMYZ_out, BRSY_out]=RWaveDataIn('GPS1156480217_Atlantic.mat',true);
    end
    
    ETMXZ_out=ETMXZ_out(300*sampf:length(ETMXZ_out));
    ETMYZ_out=ETMYZ_out(300*sampf:length(ETMYZ_out));
    ITMYZ_out=ITMYZ_out(300*sampf:length(ITMYZ_out));
    
    ETMYX_out=ETMYX_out(300*sampf:length(ETMYX_out));
    ETMYY_out=ETMYY_out(300*sampf:length(ETMYY_out));
    BRSY_out=BRSY_out(300*sampf:length(BRSY_out));
    bootVel=[];
    bootAng=[];
    if j==2
        startTime=150*sampf;
    else
        startTime=300*sampf;
    end
%     startTime=1;
%     
    if j==2
        endTime=750*sampf;
    else
        endTime=length(ETMXZ_out);
    end
    % seed=randn(1,length(ETMYZ_out));
       
    threshold=rms(ETMYZ_out);
    
    [vel, ang,bootVel,bootAng]=RWaveArray(ETMXZ_out,ETMYZ_out,ITMYZ_out,sampf,threshold,startFreq,freqStep,iter,startTime,endTime);
    
    [v,phi,el,k,bootV,bootPhi,bootEl,bootK]=...
    RWaveSingle(ETMYX_out,ETMYY_out,ETMYZ_out,BRSY_out,...
        'S',errFreq,transXErr,transYErr,transZErr,tiltErr,sampf,ang,bootAng,...
        threshold,startFreq,freqStep,iter,startTime,endTime);

%     cInd1=find((abs(v)+2*std(bootV')'>=abs(vel)'-2*std(bootVel')'));
%     cInd2=find((abs(v)-2*std(bootV')'<=abs(vel)'+2*std(bootVel')'));
%     cInd3=find(std(bootVel')<=10000);
%     cInd4=find(std(bootV')<=10000);
%     cInd5=find(std(bootVel')>=1);
%     cInd6=find(std(bootV')>=1);
%     cInd12=intersect(cInd1,cInd2);
%     cInd123=intersect(cInd12,cInd3);
%     cInd1234=intersect(cInd123,cInd4);
%     cInd12345=intersect(cInd1234,cInd5);
%     cInd=intersect(cInd12345,cInd6);
%     cInd=intersect(cInd5,cInd6);
    cInd=find(abs(v)>=0);
    v=v(cInd);
    vel=vel(cInd);
    ang=ang(cInd);
    phi=phi(cInd);
    bootV=bootV(cInd,:);
    bootVel=bootVel(cInd,:);
    bootAng=bootAng(cInd,:);
     
    singVel=[singVel; mean(v)];
    singAng=[singAng; mean(ang)];

    % figure(1)
    % plot(Rot_time(startTime:length(BRSY_out)),seriesX,Rot_time(startTime:length(BRSY_out))...
    %     ,seriesY,Rot_time(startTime:length(BRSY_out)),seriesZ,Rot_time(startTime:length(BRSY_out)),seriesRX)
    % ylim([-1e-5 100e-5])
    % xlim([400 length(ETMYZ_out)/sampf]);
    % grid on
%     [C,F]=coh2(ETMYZ_out,BRSY_out,1/8,9);
%     figure(8)
%     plot(F,C)
%     xlim([.01 .1])
%     
%     C=C(find(F<.1));
%     F=F(find(F<.1));
%     C=C(find(F>.01));
%     F=F(find(F>.01));
%     
%     coInd=[];
%     for k=1:length(v)
%         [ind ind]=min(abs(k*freqStep+startFreq-F));
%         coInd=[coInd; ind];
%     end
%     F=F(coInd);
%     C=C(coInd);
               
    figure(2)
    hold on
    % errorbar(((0:length(v)-1))*freqStep+startFreq,abs(v),-sigmaV,sigmaV)
    % errorbar(((0:length(vel)-1))*freqStep+startFreq,vel,-sigmaVel,sigmaVel)
    l=errorbar(((cInd-1)*freqStep+startFreq),abs(v),-std(bootV'),std(bootV'));
    ll=errorbar(((cInd-1)*freqStep+startFreq),vel,-std(bootVel'),std(bootVel'),'--');
    ylabel('Velocity (m/s)')
    xlabel('Frequency (Hz)')
%     legend('Single Station', 'Array','Single Station Ecuador','Array Ecuador','Single Station California','Array California')
    grid on
    box on
    % xlim([.01 .1])
    ylim([0 8e3])
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