function [vel, ang, bootVel,bootAng]=RWaveArray(ETMXZ_out,ETMYZ_out,ITMYZ_out,BRSY_out,sampf,threshold,startFreq,freqStep,iter,startTime,endTime)
%% Array
    Astop1 = 5;
    Apass  = .5;
    Astop2 =5;
    vel=[];
    ang=[];
    angTim=[];
    velTim=[];
    lagsX=[];
    lagsY=[];      
    delta_t_X=[];
    delta_t_Y=[];
    bootVel=[];
    bootAng=[];
    for i=0:iter
        freq2=(startFreq+i*freqStep);
        d = designfilt('bandpassiir', ...
          'StopbandFrequency1',(i-2)*freqStep+startFreq,'PassbandFrequency1', (i-1)*freqStep+startFreq, ...
          'PassbandFrequency2',(i+1)*freqStep+startFreq,'StopbandFrequency2', (i+2)*freqStep+startFreq, ...
          'StopbandAttenuation1',Astop1,'PassbandRipple', Apass, ...
          'StopbandAttenuation2',Astop2, ...
          'DesignMethod','butter','SampleRate',sampf);
        filtData=filter(d,ETMXZ_out-mean(ETMXZ_out));
        X=filtData(startTime:endTime);        
        filtData=filter(d,ETMYZ_out-mean(ETMYZ_out));
        Y=filtData(startTime:endTime);
        filtData=filter(d,ITMYZ_out-mean(ITMYZ_out));
        C=filtData(startTime:endTime);
        localThreshold=.1;
        
        filtData=filter(d,BRSY_out-mean(BRSY_out)); 
        seriesRX=filtData(startTime:endTime);        
        
        r2rx=[];
        tempRX=[];
        
        fitLength=floor(1/(freq2/sampf)/2);
        for j=1:floor(length(seriesRX)/fitLength)-2
           tim=(j*fitLength:(j+1)*fitLength)'./8;
%            cut=1e9*seriesX(j*fitLength:(j+1)*fitLength,o+1);           
%            tempX=[tempX (max(cut)-min(cut))];
%            cut=1e9*seriesY(j*fitLength:(j+1)*fitLength,o+1);
%            tempY=[tempY (max(cut)-min(cut))];           
           cut=1e12*seriesRX(j*fitLength:(j+1)*fitLength);
           g = fittype( @(a,b,cen_fr,x) a*sin(2*pi*cen_fr*x)+b*cos(2*pi*cen_fr*x), 'problem', 'cen_fr' );
           [myfit,st] = fit(tim,cut, g,'problem',freq2,'StartPoint', [1, 1]);
           r2rx=[r2rx;st.rsquare];
           a2=coeffvalues(myfit);
           if (abs(r2rx)>0) 
                    tempRX=[tempRX (a2(2)+i*a2(1))/1e3];
                    tim=(j*fitLength:(j+1)*fitLength)'./8;%                        
           end
        end
%         if i==2
%             figure(11)            
%             plot((startTime:endTime)/sampf,X,(startTime:endTime)/sampf,Y+1e-6,(startTime:endTime)/sampf,C+2e-6)
%             legend('X','Y','C')
%         end
%         if i==6
%             figure(12)            
%             plot((startTime:endTime)/sampf,X,(startTime:endTime)/sampf,Y+1e-6,(startTime:endTime)/sampf,C+2e-6)
%             legend('X','Y','C')
%         end
    %     injAng=150*pi/180;
    %     X=filter(d,seed(ceil(100*cos(injAng))+150:length(ETMXZ_out)+ceil(100*cos(injAng))+150))';
    %     Y=filter(d,seed(ceil(100*sin(injAng))+150:length(ETMXZ_out)+ceil(100*sin(injAng))+150))';
    %     C=filter(d,seed(1+150:length(ETMXZ_out)+1+150))';
    %     injVel=7000;
    %     X=filter(d,seed(4e3/(sqrt(2)*injVel)*10*sampf+150:length(ETMXZ_out)+4e3/(sqrt(2)*injVel)*10*sampf+150))';
    %     Y=filter(d,seed(4e3/(sqrt(2)*injVel)*10*sampf+150:length(ETMXZ_out)+4e3/(sqrt(2)*injVel)*10*sampf+150))';
    %     C=filter(d,seed(1+150:length(ETMXZ_out)+1+150))';

    delta_t_X=[];
    delta_t_Y=[];
    if endTime<=6500
        len=250*sampf;
    else
        len=250*sampf;
    end
    for j=1:floor(length(X)/len)-1
%         if ((max(C(j*len:(j+1)*len))-min(C(j*len:(j+1)*len)))>=threshold && (max(C(j*len:(j+1)*len))-min(C(j*len:(j+1)*len)))>=localThreshold)
        if(max(abs(tempRX(floor(j*len/fitLength):floor((j+1)*len/fitLength))))>=localThreshold)

            [crossY,~] = xcorr(C(j*len:(j+1)*len),Y(j*len:(j+1)*len));
            [crossX,lags]=xcorr(C(j*len:(j+1)*len),X(j*len:(j+1)*len));

            crossX=abs(crossX);
            crossY=abs(crossY);
            
            peak=crossX(floor(find(crossX==max(crossX))-1/freq2):floor(find(crossX==max(crossX))+1/freq2));
            peakLags=lags(floor(find(crossX==max(crossX))-1/freq2):floor(find(crossX==max(crossX))+1/freq2))';
            [myfit,s]=polyfit(peakLags,peak,2);  
            delta_t_X=[delta_t_X -myfit(2)/(2*myfit(1))/sampf];
            
            peak=crossY(floor(find(crossY==max(crossY))-1/freq2):floor(find(crossY==max(crossY))+1/freq2));
            peakLags=lags(floor(find(crossY==max(crossY))-1/freq2):floor(find(crossY==max(crossY))+1/freq2))';
            [myfit,s]=polyfit(peakLags,peak,2);
            delta_t_Y=[delta_t_Y -myfit(2)/(2*myfit(1))/sampf];

        end
    end   
    tempBAng=[];
    tempBVel=[];
    for k=0:1e4
        array=[delta_t_X; delta_t_Y];
        bootArray=bootstrapData(array');
        if length(bootArray)>0
            bootTX=bootArray(:,1)';
            bootTY=bootArray(:,2)';   
            bootTX=mean(bootTX);
            bootTY=mean(bootTY);
        else
            bootTX=nan;
            bootTY=nan;
        end
        tempBAng=[tempBAng; atan2(bootTY,bootTX)*180/pi];      
        tempBVel=[tempBVel; 4e3./sqrt((bootTX).^2+(bootTY).^2)];
    end     
    bootVel=[bootVel; tempBVel'];
    bootAng=[bootAng; tempBAng'];
%     if i==7
%         figure(4)
%         plot(4e3./sqrt((delta_t_X).^2+(delta_t_Y).^2));
%     end
    delta_t_X=mean(delta_t_X);
    delta_t_Y=mean(delta_t_Y);
    {atan2(delta_t_Y,delta_t_X)*180/pi,4e3./sqrt((delta_t_X).^2+(delta_t_Y).^2)}
    ang=[ang atan2(delta_t_Y,delta_t_X)*180/pi];    
    vel=[vel 4e3./sqrt((delta_t_X).^2+(delta_t_Y).^2)];
    end
end