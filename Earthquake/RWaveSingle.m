function [v,phi,el,k,bootV,bootPhi,bootEl,bootK]=...
    RWaveSingle(ETMYX_out,ETMYY_out,ETMYZ_out,BRSY_out,...
    quad,errFreq,transXErr,transYErr,transZErr,tiltErr,sampf,ang,bootAng,...
    threshold,startFreq,freqStep,iter,startTime,endTime)

    %% Single Station
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
    Astop1 = 5;
    Apass  = .5;
    Astop2 = 5;
    for o=0:iter
        tempBV=[];
        tempBPhi=[];
        tempBEl=[];
        tempBK=[];
        
        freq1=(startFreq+o*freqStep);
        [indx indx]=min(abs(errFreq-freq1));
        sigmaX=transXErr(indx);
        sigmaY=transYErr(indx);
        sigmaZ=transZErr(indx);
        sigmaRX=tiltErr(indx);
        d = designfilt('bandpassiir', ...
          'StopbandFrequency1',(o-2)*freqStep+startFreq,'PassbandFrequency1', (o-1)*freqStep+startFreq, ...
          'PassbandFrequency2',(o+1)*freqStep+startFreq,'StopbandFrequency2', (o+2)*freqStep+startFreq, ...
          'StopbandAttenuation1',Astop1,'PassbandRipple', Apass, ...
          'StopbandAttenuation2',Astop2, ...
          'DesignMethod','butter','SampleRate',sampf);

        filtData=filter(d,ETMYX_out-mean(ETMYX_out));
        seriesX=[seriesX filtData(startTime:endTime)];
        filtData=filter(d,ETMYY_out-mean(ETMYY_out)); 
        seriesY=[seriesY filtData(startTime:endTime)];
        filtData=filter(d,ETMYZ_out-mean(ETMYZ_out));
        seriesZ=[seriesZ filtData(startTime:endTime)];
        localThreshold1=max(abs(filtData))*.7*0;
        filtData=filter(d,BRSY_out-mean(BRSY_out)); 
        seriesRX=[seriesRX filtData(startTime:endTime)];        
        localThreshold2=max(abs(filtData))*.7*0;
        tempX=[];
        tempY=[];
        tempZ=[];
        tempRX=[];      
%         for j=1:floor(length(seriesX(:,i+1))*(freq1/sampf))-4
%            cut=seriesX(floor(j/(freq1/sampf)):floor((j+1)/(freq1/sampf)),i+1);
%            tempX=[tempX (max(cut)-min(cut))];
%            cut=seriesY(floor(j/(freq1/sampf)):floor((j+1)/(freq1/sampf)),i+1);
%            tempY=[tempY (max(cut)-min(cut))];
%            cut=seriesZ(floor(j/(freq1/sampf)):floor((j+1)/(freq1/sampf)),i+1);
%            tempZ=[tempZ (max(cut)-min(cut))];
%            cut=seriesRX(floor(j/(freq1/sampf)):floor((j+1)/(freq1/sampf)),i+1);
%            tempRX=[tempRX (max(cut)-min(cut))];
%         end 
        for j=1:floor(length(seriesX(:,o+1))*(freq1/sampf))-4
           tim=(floor(j/(freq1/sampf)):floor((j+1)/(freq1/sampf)))'./8;
           cut=1e9*seriesX(floor(j/(freq1/sampf)):floor((j+1)/(freq1/sampf)),o+1);           
           tempX=[tempX (max(cut)-min(cut))];
           cut=1e9*seriesY(floor(j/(freq1/sampf)):floor((j+1)/(freq1/sampf)),o+1);
           tempY=[tempY (max(cut)-min(cut))];
           cut=1e9*seriesZ(floor(j/(freq1/sampf)):floor((j+1)/(freq1/sampf)),o+1);
           myfit = fit(tim,cut, 'a*sin(2*pi*freq1*x)+b*cos(2*pi*freq1*x)','StartPoint', [1, 1, freq1]);
           a=coeffvalues(myfit);
           tempZ=[tempZ a(2)+i*a(1)];
           cut=1e9*seriesRX(floor(j/(freq1/sampf)):floor((j+1)/(freq1/sampf)),o+1);
           myfit = fit(tim,cut, 'a*sin(2*pi*freq1*x)+b*cos(2*pi*freq1*x)','StartPoint', [1, 1, freq1]);
           a=coeffvalues(myfit);
           tempRX=[tempRX a(2)+i*a(1)];
        end
        for p=1:1e4
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
            array=[tempX; tempY; tempZ; tempRX]; 
            bootArray=bootstrapData(array');
            btempX=bootArray(:,1)';
            btempY=bootArray(:,2)';
            btempZ=bootArray(:,3)';
            btempRX=bootArray(:,4)';
            for l=1:min([length(btempX) length(btempY) length(btempZ) length(btempRX)])
                if(abs(btempZ(l))>=threshold && abs(btempZ(l))>=localThreshold1 && abs(btempRX(l))>=localThreshold2)
                    avgPhi=avgPhi+atan2(btempY(l),btempX(l))*180/pi;
                    avgK=avgK+btempRX(l)./btempZ(l)./sin(atan2(btempY(l),btempX(l)));
                    avgV=avgV+2*pi*freq1.*abs(btempZ(l))./abs(btempRX(l)).*cos(angle(btempZ(l))-angle(btempRX(l))).*sin(bootAng(3,floor(1 + (length(bootAng)-1).*rand(1)))*pi/180);%ang(i+1)*pi/180%(atan2(btempY(l),btempX(l)));%
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
        avgPhi=0;
        avgK=0;
        avgEl=0;
        avgV=0;
        pltV=[];
        for l=1:min([length(tempX) length(tempY) length(tempZ) length(tempRX)])
            if(abs(tempZ(l))>=threshold && abs(tempRX(l))>=localThreshold2 && abs(tempZ(l))>=localThreshold1)
    %         if l>=20
                avgPhi=avgPhi+atan2(tempY(l),tempX(l))*180/pi;
                avgK=avgK+tempRX(l)./tempZ(l)./sin(atan2(tempY(l),tempX(l)));
                avgV=avgV+2*pi*freq1.*abs(tempZ(l))./abs(tempRX(l)).*cos(angle(tempZ(l))-angle(tempRX(l))).*sin(ang(o+1)*pi/180);%(atan2(tempY(l),tempX(l)));%
                pltV=[pltV; 2*pi*freq1.*tempZ(l)./(tempRX(l)).*sin(ang(o+1)*pi/180)];
                avgEl=avgEl+acot(tempZ(l)/tempY(l).*sin(atan2(tempY(l),tempX(l))));                
                N=N+1;
            end
        end   
%         if i==7
%             figure(5)
%             plot(pltV)
%         end
        avgPhi=avgPhi/N;
        avgK=avgK/N;
        avgV=avgV/N;
        avgEl=avgEl/N;
        
        phi=[phi; avgPhi];
        k=[k; avgK];
        v=[v; avgV];
        el=[el; avgEl];
    end
end