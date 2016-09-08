for j=0:2
    sampf =8;
    startFreq=0.015;
    freqStep=.005;
    iter=(.1-startFreq)/freqStep;
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
%     RWaveDataIn('GPS1149331885_6_2Earthquake.mat');
    bootVel=[];
    bootAng=[];
    threshold=rms(ETMYZ_out)/2;
    % seed=randn(1,length(ETMYZ_out));

    [vel, ang, sigmaVel,sigmaAng,bootVel,bootAng]=RWaveArray(ETMXZ_out,ETMYZ_out,ITMYZ_out,sampf,threshold,startFreq,freqStep,iter);
    freqStep=0.1/length(ang);
    [v,phi,el,k,sigmaV,sigmaPhi,bootV,bootPhi,bootEl,bootK]=...
    RWaveSingle(ETMYX_out,ETMYY_out,ETMYZ_out,BRSY_out,...
        'S',errFreq,transXErr,transYErr,transZErr,tiltErr,sampf,ang,threshold,startFreq,freqStep,iter);

    cInd1=find((abs(v)+std(bootV')'>=abs(vel)'-std(bootVel')'));
    cInd2=find((abs(v)-std(bootV')'<=abs(vel)'+std(bootVel')'));
    cInd3=find(std(bootVel')<=1000);
    cInd4=find(std(bootV')<=1000);
    cInd12=intersect(cInd1,cInd2);
    cInd123=intersect(cInd12,cInd3);
    cInd=intersect(cInd123,cInd4);
    v=v(cInd);
    vel=vel(cInd);
    bootV=bootV(cInd,:);
    bootVel=bootVel(cInd,:);

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
    l=errorbar(cInd*freqStep+startFreq,abs(v),-std(bootV'),std(bootV'))
    ll=errorbar(cInd*freqStep+startFreq,vel,-std(bootVel'),std(bootVel'))
    ylabel('Velocity (m/s)')
    xlabel('Frequency (Hz)')
    legend('Single Station','Array','Single Station Bootstrapping','Array Bootstrapping')
    grid on
    % xlim([.01 .1])
    ylim([0 1e4])
    set(gca,'FontSize',12)
    set(l,'LineWidth',2)
    set(ll,'LineWidth',2)
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
end
