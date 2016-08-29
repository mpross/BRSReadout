sampf =8;
bootV=[];
bootPhi=[];
bootAng=[];
bootVel=[];
[errFreq,transXErr,transYErr,transZErr,tiltErr]=RWaveMeasErr;
[ETMXZ_out, ITMYZ_out, ETMYX_out, ETMYY_out, ETMYZ_out, BRSY_out]=...
      RWaveDataIn('GPS1144888165_7_8Earthquake.mat');
%     RWaveDataIn('GPS1143962787_6_9Earthquake.mat');
%     RWaveDataIn('GPS1149323500_Quite.mat');
%     RWaveDataIn('GPS1149331885_6_2Earthquake.mat');       
%     RWaveDataIn('GPS1149581095_5_2Earthquake.mat');      
   
seed=randn(1,length(ETMYZ));
for i=1:100
    [v,phi,el,k,sigmaV,sigmaPhi]=...
    RWaveSingle(bootstrapData(ETMYX_out),bootstrapData(ETMYY_out),...
        bootstrapData(ETMYZ_out),bootstrapData(BRSY_out),...
        'S',errFreq,transXErr,transYErr,transZErr,tiltErr,sampf);

    [vel, ang, sigmaVel,sigmaAng]=RWaveArray(...
        bootstrapData(ETMXZ_out),bootstrapData(ETMYZ_out),bootstrapData(ITMYZ_out),sampf);

    bootV=[bootV; v'];
    bootPhi=[bootPhi; phi'];
    bootVel=[bootVel; vel'];
    bootAng=[bootAng; ang'];
end

[v,phi,el,k,sigmaV,sigmaPhi]=...
RWaveSingle(ETMYX_out,ETMYY_out,ETMYZ_out,BRSY_out,...
    'S',errFreq,transXErr,transYErr,transZErr,tiltErr,sampf);

[vel, ang, sigmaVel,sigmaAng]=RWaveArray(ETMXZ_out,ETMYZ_out,ITMYZ_out,sampf);


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
