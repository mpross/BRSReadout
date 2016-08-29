function [v,phi,el,k,sigmaV,sigmaPhi,bootV,bootPhi,bootEl,bootK]=...
    RWaveSingle(ETMYX_out,ETMYY_out,ETMYZ_out,BRSY_out,...
    quad,errFreq,transXErr,transYErr,transZErr,tiltErr,sampf)
    %% Single Station
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

    bootV=[];
    bootPhi=[];
    bootEl=[];
    bootK=[];
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
%     sigmaX=3e-8;
%     sigmaY=3e-8;
%     sigmaZ=3e-8;
%     % sigmaRX=8e-10;
%     sigmaRX=2e-10;
    Astop1 = 10;
    Apass  = 1;
    Astop2 = 10;
    freqStep1=0.001;
    startFreq1=0.03;
    startTime=01*sampf;
    % endTime=800*sampf;
    endTime=length(ETMYZ_out);
    % startTime=2000*sampf;
    % endTime=2600*sampf;
    % sgnZ=sign(ETMYZ_out(find(abs(ETMYZ_out-mean(ETMYZ_out))==max(abs(ETMYZ_out-mean(ETMYZ_out)))))-mean(ETMYZ_out));
    % sgnX=sgnZ*sign(ETMYX_out(find(abs(ETMYZ_out-mean(ETMYZ_out))==max(abs(ETMYZ_out-mean(ETMYZ_out)))))-mean(ETMYX_out));
    % sgnY=sgnZ*sign(ETMYY_out(find(abs(ETMYZ_out-mean(ETMYZ_out))==max(abs(ETMYZ_out-mean(ETMYZ_out)))))-mean(ETMYY_out));
    % sgnRX=sgnZ*sign(BRSY_out(find(abs(ETMYZ_out-mean(ETMYZ_out))==max(abs(ETMYZ_out-mean(ETMYZ_out)))))-mean(BRSY_out));
    % sgnZ=1;
    for i=0:100
        tempBV=[];
        tempBPhi=[];
        tempBEl=[];
        tempBK=[];
        
        freq1=(startFreq1+i*freqStep1);
        [indx indx]=min(abs(errFreq-freq1));
        sigmaX=transXErr(indx);
        sigmaY=transYErr(indx);
        sigmaZ=transZErr(indx);
        sigmaRX=tiltErr(indx);
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
        for i=1:100
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
                btempX=bootstrapData(tempX);
                btempY=bootstrapData(tempY);
                btempZ=bootstrapData(tempZ);
                btempRX=bootstrapData(tempRX);
                for l=1:min([length(btempX) length(btempY) length(btempZ) length(btempRX)])
                    if(abs(btempZ(l))>=threshold)
                        avgPhi=avgPhi+atan2(btempY(l),btempX(l))*180/pi;
                        avgK=avgK+btempRX(l)./btempZ(l)./sin(atan2(btempY(l),btempX(l)));
                        avgV=avgV+2*pi*freq1.*btempZ(l)./(btempRX(l)).*sin(atan2(btempY(l),btempX(l)));%(ang(i+1)*pi/180+180);%
            %             avgV=avgV+2*pi*freq1.*btempY(l)./(btempRX(l)).*btempZ(l)./btempY(l).*sin(atan2(btempY(l),btempX(l)));
                        avgEl=avgEl+acot(btempZ(l)/btempY(l).*sin(atan2(btempY(l),btempX(l))));
                        N=N+1;
                    end
                end   

                avgPhi=avgPhi/N;
                avgK=avgK/N;
                avgV=avgV/N;
                avgEl=avgEl/N;

            tempBV=[tempBV; avgV];
            tempBPhi=[tempBPhi; avgPhi];
            tempBEl=[tempBEl; avgEl];
            tempBK=[tempBK; avgK];
        end
        bootV=[bootV; tempBV'];
        bootPhi=[bootPhi; tempBPhi'];
        bootEl=[bootEl; tempBEl'];
        bootK=[bootK; tempBK'];
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
                sumSigmaPhi=sumSigmaPhi+(tempY(l).^2.*sigmaX.^2+tempX(l).^2.*sigmaY.^2)/(tempX(l).^2+tempY(l).^2)^2;
                dVdZ=2*pi*freq1./(tempRX(l)).*sin(atan2(tempY(l),tempX(l)));
                dVdRX=-2*pi*freq1.*tempZ(l)./(tempRX(l)^2).*sin(atan2(tempY(l),tempX(l)));
                dVdY=2*pi*freq1.*tempZ(l)./(tempRX(l)).*cos(atan2(tempY(l),tempX(l))).*(tempX(l)/(tempX(l)^2+tempY(l)^2));
                dVdX=-2*pi*freq1.*tempZ(l)./(tempRX(l)).*cos(atan2(tempY(l),tempX(l))).*(tempY(l)/(tempX(l)^2+tempY(l)^2));
                sumSigmaV=sumSigmaV+dVdX.^2.*sigmaX.^2+dVdY.^2.*sigmaY.^2+dVdZ.^2.*sigmaZ.^2+dVdRX.^2.*sigmaRX.^2;
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
        sigmaV=[sigmaV, sqrt(sumSigmaV)/N];
    end

    % 
    % v=v(find(abs(v-mean(v))<=3*std(abs(v))));
    % 
    % phi=phi(find(abs(phi-mean(phi))<=3*std(abs(phi))));
    % 
    % el=el(find(abs(el-mean(el))<=3*std(abs(el))));

end