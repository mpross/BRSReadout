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
    r2zavg=[];
    r2rxavg=[];
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
        localThreshold1=max(abs(filtData))*.0;
        filtData=filter(d,BRSY_out-mean(BRSY_out)); 
        seriesRX=[seriesRX filtData(startTime:endTime)];        
        localThreshold2=max(abs(filtData))*.0;
        FitData=seriesZ;
        tempX=[];
        tempY=[];
        tempZ=[];
        tempRX=[];
        r2z=[];
        r2rx=[];
%         for j=1:floor(length(seriesX(:,i+1))*(freq1/sampf))-4
%            cut=seriesX(j*fitLength:(j+1)*fitLength,i+1);
%            tempX=[tempX (max(cut)-min(cut))];
%            cut=seriesY(j*fitLength:(j+1)*fitLength,i+1);
%            tempY=[tempY (max(cut)-min(cut))];
%            cut=seriesZ(j*fitLength:(j+1)*fitLength,i+1);
%            tempZ=[tempZ (max(cut)-min(cut))];
%            cut=seriesRX(j*fitLength:(j+1)*fitLength,i+1);
%            tempRX=[tempRX (max(cut)-min(cut))];
%         end 
        fitLength=floor(1/(freq1/sampf)/2);
        for j=1:floor(length(seriesZ(:,o+1))/fitLength)-2
           tim=(j*fitLength:(j+1)*fitLength)'./8;
%            cut=1e9*seriesX(j*fitLength:(j+1)*fitLength,o+1);           
%            tempX=[tempX (max(cut)-min(cut))];
%            cut=1e9*seriesY(j*fitLength:(j+1)*fitLength,o+1);
%            tempY=[tempY (max(cut)-min(cut))];
           cut=1e9*seriesZ(j*fitLength:(j+1)*fitLength,o+1);
           g = fittype( @(a,b,cen_fr,x) a*sin(2*pi*cen_fr*x)+b*cos(2*pi*cen_fr*x), 'problem', 'cen_fr' );
           [myfit,st] = fit(tim,cut, g,'problem',freq1,'StartPoint', [1, 1]);
           r2z=[r2z;st.rsquare];
           a1=coeffvalues(myfit);
           FitData(j*fitLength:(j+1)*fitLength) = a1(1)*sin(2*pi*freq1.*tim)+a1(2)*cos(2*pi*freq1.*tim);          
           cut=1e12*seriesRX(j*fitLength+floor(sampf/(4*freq1)):(j+1)*fitLength+floor(sampf/(4*freq1)),o+1);
           g = fittype( @(a,b,cen_fr,x) a*sin(2*pi*cen_fr*x)+b*cos(2*pi*cen_fr*x), 'problem', 'cen_fr' );
           [myfit,st] = fit(tim,cut, g,'problem',freq1,'StartPoint', [1, 1]);
           r2rx=[r2rx;st.rsquare];
           a2=coeffvalues(myfit);
           if (abs(r2rx)>0) 
               if (abs(r2z)>0)
                   if(abs(cos(angle(a1(2)+i*a1(1))-angle((a2(2)+i*a2(1)))))>cos(10*pi/180))
                       tempZ=[tempZ a1(2)+i*a1(1)];
                       tempRX=[tempRX (a2(2)+i*a2(1))/1e3];
                       tim=(j*fitLength:(j+1)*fitLength)'./8;
                       cut=1e9*seriesX(j*fitLength:(j+1)*fitLength,o+1);           
                       tempX=[tempX (max(cut)-min(cut))];
                       cut=1e9*seriesY(j*fitLength:(j+1)*fitLength,o+1);
                       tempY=[tempY (max(cut)-min(cut))];
                       cut=1e9*seriesZ(j*fitLength:(j+1)*fitLength,o+1);
                   end
               end
           end
        end
%         figure(8)
%         lll=plot((1:length(seriesZ)),seriesZ*1e6,(1:length(FitData)),FitData/1e3,'--')
%         ylabel('Vertical Displacement (um)')
%         xlabel('Time (s)')
%         legend('Data','Fit')
%         set(lll,'LineWidth',1.2)
%         set(gca,'FontSize',12)
%         grid on
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
                    avgV=avgV+2*pi*freq1.*abs(btempZ(l))./abs(btempRX(l)).*sin(bootAng(3,floor(1 + (length(bootAng)-1).*rand(1)))*pi/180);%ang(i+1)*pi/180%(atan2(btempY(l),btempX(l)));%*cos(angle(btempZ(l))-angle(btempRX(l)))
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
                avgV=avgV+2*pi*freq1.*abs(tempZ(l))./abs(tempRX(l)).*sin(ang(o+1)*pi/180);%(atan2(tempY(l),tempX(l)));%.*cos(angle(tempZ(l))-angle(tempRX(l)))
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
        r2zavg=[mean(r2z);r2zavg];
        r2rxavg=[mean(r2rx);r2rxavg];
    end
    if o==0
        hold on
        figure(8)
        plot(r2zavg)
        plot(r2rxavg)
        ylabel('R^2 Value')
        xlabel('Time Cut')
        legend('Vertical','Rotational')
        set(l,'LineWidth',1.2)
        set(ll,'LineWidth',1.2)
        set(gca,'FontSize',12)
        grid on
    end
end