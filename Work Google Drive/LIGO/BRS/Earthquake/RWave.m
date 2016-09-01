sampf =8;
freqStep=0.001;
startFreq=0.03;
[errFreq,transXErr,transYErr,transZErr,tiltErr]=RWaveMeasErr;
[ETMXZ_out, ITMYZ_out, ETMYX_out, ETMYY_out, ETMYZ_out, BRSY_out]=...    
    RWaveDataIn('GPS1143962787_6_9Earthquake.mat');
%       RWaveDataIn('GPS1144888165_7_8Earthquake.mat');
%     RWaveDataIn('GPS1149323500_Quite.mat');
%     RWaveDataIn('GPS1149331885_6_2Earthquake.mat');       
%     RWaveDataIn('GPS1149581095_5_2Earthquake.mat');      
bootVel=[];
bootAng=[];
% seed=randn(1,length(ETMYZ_out));

[v,phi,el,k,sigmaV,sigmaPhi,bootV,bootPhi,bootEl,bootK]=...
RWaveSingle(ETMYX_out,ETMYY_out,ETMYZ_out,BRSY_out,...
    'E',errFreq,transXErr,transYErr,transZErr,tiltErr,sampf);

[vel, ang, sigmaVel,sigmaAng,bootVel,bootAng]=RWaveArray(ETMXZ_out,ETMYZ_out,ITMYZ_out,sampf);


% figure(1)
% plot(Rot_time(startTime:length(BRSY_out)),seriesX,Rot_time(startTime:length(BRSY_out))...
%     ,seriesY,Rot_time(startTime:length(BRSY_out)),seriesZ,Rot_time(startTime:length(BRSY_out)),seriesRX)
% ylim([-1e-5 100e-5])
% xlim([400 length(ETMYZ_out)/sampf]);
% grid on

figure(2)
hold on
errorbar(((0:length(v)-1))*freqStep+startFreq,abs(v),-sigmaV,sigmaV)
errorbar(((0:length(vel)-1))*freqStep+startFreq,vel,-sigmaVel,sigmaVel)
errorbar(((0:length(phi)-1))*freqStep+startFreq,mean(bootV'),-std(bootV'),std(bootV'))
errorbar(((0:length(phi)-1))*freqStep+startFreq,mean(bootVel'),-std(bootVel'),std(bootVel'))
ylabel('Velocity (m/s)')
xlabel('Frequency (Hz)')
legend('Single Station','Array','Single Station Bootstrapping','Array Bootstrapping')
grid on
% xlim([.01 .1])
ylim([0 1e4])
hold off

figure(3)
hold on
errorbar(((0:length(phi)-1))*freqStep+startFreq,phi,-sigmaPhi,sigmaPhi)
errorbar(((0:length(ang)-1))*freqStep+startFreq,ang,-sigmaAng,sigmaAng)
errorbar(((0:length(phi)-1))*freqStep+startFreq,mean(bootPhi'),-std(bootPhi'),std(bootPhi'))
errorbar(((0:length(phi)-1))*freqStep+startFreq,mean(bootAng'),-std(bootAng'),std(bootAng'))
ylabel('Angle of Incidence (degrees)')
xlabel('Frequency (Hz)')
legend('Single Station','Array','Single Station Bootstrapping','Array Bootstrapping')
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
% plot(((0:length(phi)-1))*freqStep+startFreq,el/pi)
% ylabel('Rayleigh Ellipticity (pi rad)')
% xlabel('Frequency (Hz)')
% grid on
% xlim([.02 .0675])
