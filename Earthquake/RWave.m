close all
singVel=[0];
singAng=[0];
avgVel=containers.Map('ValueType','any','KeyType','double');
avgV=containers.Map('ValueType','any','KeyType','double');
avgVelErr=containers.Map('ValueType','any','KeyType','double');
avgVErr=containers.Map('ValueType','any','KeyType','double');
fileName={'GPS1143962787_6_9Earthquake.mat','GPS1144888165_7_8Earthquake.mat','GPS1149581095_5_2Earthquake.mat',...
    'GPS1149331885_6_2Earthquake.mat','GPS1156069817_Burma.mat','GPS1156480217_Atlantic.mat','GPS1156782617_NewZealand.mat',...
    'GPS1157149817_Russia.mat','GPS1163070017_7k_EQ5.mat','GPS1163797217_7k_EQ6.mat','GPS1166005817_10k_EQ7_PapNG.mat','GPS1155000017_7k_EQ8_NewC.mat'};
newArray=[false, false, false, false, true, true, true, true, true, true, true, true, true];
sampf =8;
% startArray=[1000, 800, 1, 750, 3000, 2000, 3000, 1, 2500, 1, 4500].*sampf;
% endArray=[2000, 1700, 1000, 1000, 4200, 4500, 4500, 3200, 5500, 1, 5500].*sampf;
% startArray=[500, 800, 1, 750, 3000, 2500, 2500, 1].*sampf;
% endArray=[1250, 1700, 1000, 1000, 4200, 3250, 3250, 3200].*sampf;
% for j=0:7
% for j=[0 1 4 5 6]
for j=[0 5 8 10]
% for j=[0]
% for j=10
    startFreq=0.025;
    freqStep=.005;
    iter=floor((.1-startFreq)/freqStep);
    [errFreq,transXErr,transYErr,transZErr,tiltErr]=RWaveMeasErr;
    [ETMXZ_out, ITMYZ_out, ETMYX_out, ETMYY_out, ETMYZ_out, BRSY_out]=RWaveDataIn(fileName{j+1},newArray(j+1));
    
    ETMXZ_out=ETMXZ_out(300*sampf:length(ETMXZ_out));
    ETMYZ_out=ETMYZ_out(300*sampf:length(ETMYZ_out));
    ITMYZ_out=ITMYZ_out(300*sampf:length(ITMYZ_out));
    
    ETMYX_out=ETMYX_out(300*sampf:length(ETMYX_out));
    ETMYY_out=ETMYY_out(300*sampf:length(ETMYY_out));
    BRSY_out=BRSY_out(300*sampf:length(BRSY_out));
    bootVel=[];
    bootAng=[];
    startTime=1;
%     startTime=startArray(j+1);
%     endTime=endArray(j+1);
    endTime=length(ETMXZ_out);
%     % seed=randn(1,length(ETMYZ_out));
    threshold=rms(ETMYZ_out(startTime:endTime))*0;
    
    [vel, ang,bootVel,bootAng]=RWaveArray(ETMXZ_out,ETMYZ_out,ITMYZ_out,BRSY_out,sampf,threshold,startFreq,freqStep,iter,startTime,endTime);
    
    [v,phi,el,k,bootV,bootPhi,bootEl,bootK]=...
    RWaveSingle(ETMYX_out,ETMYY_out,ETMYZ_out,BRSY_out,...
        'S',errFreq,transXErr,transYErr,transZErr,tiltErr,sampf,ang,bootAng,...
        threshold,startFreq,freqStep,iter,startTime,endTime);

%     cInd1=find((abs(v)+2*std(bootV')'>=abs(vel)'-2*std(bootVel')'));
%     cInd2=find((abs(v)-2*std(bootV')'<=abs(vel)'+2*std(bootVel')'));
%     cInd3=find(std(bootVel')<=10000);
%     cInd4=find(std(bootV')<=10000);
    cInd5=find(std(bootVel')>=1);
    cInd6=find(std(bootV')>=1);
%     cInd12=intersect(cInd1,cInd2);
%     cInd123=intersect(cInd12,cInd3);
%     cInd1234=intersect(cInd123,cInd4);
%     cInd12345=intersect(cInd1234,cInd5);
%     cInd=intersect(cInd12345,cInd6);
    cInd=intersect(cInd5,cInd6);
%     cInd=find(abs(v)>=0);
    v=v(cInd);
    vel=vel(cInd);
    ang=ang(cInd);
%     phi=phi(cInd);
    bootV=bootV(cInd,:);
    bootVel=bootVel(cInd,:);
    bootAng=bootAng(cInd,:);
     
    singVel=[singVel; mean(v)];
    singAng=[singAng; mean(ang)]

   
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
        figure(5)
        hold on
        % errorbar(((0:length(v)-1))*freqStep+startFreq,abs(v),-sigmaV,sigmaV)
        % errorbar(((0:length(vel)-1))*freqStep+startFreq,vel,-sigmaVel,sigmaVel)
        l=errorbar(((cInd-1)*freqStep+startFreq),abs(v),-std(bootV'),std(bootV'));
        ll=errorbar(((cInd-1)*freqStep+startFreq),vel,-std(bootVel'),std(bootVel'),'--');
        ylabel('Velocity (m/s)')
        xlabel('Frequency (Hz)')
        legend('Single Station', 'Array')
        grid on
        box on
        % xlim([.01 .1])
%         ylim([0 8e3])
%         set(gca,'FontSize',12)
%         set(l,'LineWidth',1.2)
%         set(ll,'LineWidth',1.2)
        hold off
        
        com=intersect(cInd,cell2mat(avgV.keys));
        for k=1:length(com)
            avgV(com(k))=[avgV(com(k)) v(k)];
            avgVel(com(k))=[avgVel(com(k)) vel(k)];
            avgVErr(com(k))=[avgVErr(com(k)) std(bootV(k,:)')^2];
            avgVelErr(com(k))=[avgVelErr(com(k)) std(bootVel(k,:)')^2];
        end
        new=setxor(cInd,intersect(cInd,cell2mat(avgV.keys)));
        for k=1:length(new)
            avgV(new(k))=v(k);
            avgVel(new(k))=vel(k);
            avgVErr(new(k))=std(bootV(k,:)')^2;
            avgVelErr(new(k))=std(bootVel(k,:)')^2;
        end
        
    end
%     figure(3)
%     hold on
%     % errorbar(((0:length(phi)-1))*freqStep+startFreq,phi,-sigmaPhi,sigmaPhi)
%     % errorbar(((0:length(ang)-1))*freqStep+startFreq,ang,-sigmaAng,sigmaAng)
%     errorbar(((0:length(phi)-1))*freqStep+startFreq,phi,-std(bootPhi'),std(bootPhi'))
%     errorbar(((0:length(ang)-1))*freqStep+startFreq,ang,-std(bootAng'),std(bootAng'))
%     ylabel('Angle of Incidence (degrees)')
%     xlabel('Frequency (Hz)')
%     legend('Single Station','Array','Single Station Bootstrapping','Array Bootstrapping')
%     grid on
%     % xlim([.01 .1])
%     ylim([-200 200])
%     hold off
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
cInd=cell2mat(avgV.keys);
for k=1:length(avgV.keys)
    N=length(avgV(cInd(k)));
    avgV(cInd(k))=sum(abs(avgV(cInd(k))))/N;
    avgVel(cInd(k))=sum(abs(avgVel(cInd(k))))/N;
    avgVErr(cInd(k))=sqrt(sum(avgVErr(cInd(k))))/N;
    avgVelErr(cInd(k))=sqrt(sum(avgVelErr(cInd(k))))/N;
end


figure(6)
hold on
l=errorbar(((cInd-1)*freqStep+startFreq)+10^-4,cell2mat(avgV.values)/1000,-cell2mat(avgVErr.values)/1000,cell2mat(avgVErr.values)/1000);
ll=errorbar(((cInd-1)*freqStep+startFreq),cell2mat(avgVel.values)/1000,-cell2mat(avgVelErr.values)/1000,cell2mat(avgVelErr.values)/1000,'--');
ylabel('Average Phase Velocity (km/s)')
xlabel('Frequency (Hz)')
legend('Single Station','Array')
set(l,'LineWidth',1.2)
set(ll,'LineWidth',1.2)
set(gca,'FontSize',12)
set(gca,'XTick',0.005*(0:100))
set(gca,'YTick',.5*(0:100))
grid on
box on
