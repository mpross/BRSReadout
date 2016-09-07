function [vel, ang, sigmaVel,sigmaAng,bootVel,bootAng]=RWaveArray(ETMXZ_out,ETMYZ_out,ITMYZ_out,sampf,threshold,startFreq,freqStep,iter)
%% Array
    Astop1 = 10;
    Apass  = 1;
    Astop2 = 10;
    vel=[];
    ang=[];
    angTim=[];
    velTim=[];
    lagsX=[];
    lagsY=[];    
    sigmaTX=0.0231;
    sigmaTY=0.0325;
    sigmaAng=[];
    sigmaVel=[];    
   startTime=01*sampf;
    % endTime=800*sampf;
    endTime=length(ETMYZ_out);
    delta_t_X=[];
    delta_t_Y=[];
    bootVel=[];
    bootAng=[];
    for i=0:iter
        freq2=(startFreq+i*freqStep);
    %     [bb,aa] = butter(2,[(0.032+i*0.005)/sampf, (0.032+(i+1)*0.005)/sampf]);
        d = designfilt('bandpassiir', ...
          'StopbandFrequency1',(i-2)*freqStep+startFreq,'PassbandFrequency1', (i-1)*freqStep+startFreq, ...
          'PassbandFrequency2',(i+1)*freqStep+startFreq,'StopbandFrequency2', (i+2)*freqStep+startFreq, ...
          'StopbandAttenuation1',Astop1,'PassbandRipple', Apass, ...
          'StopbandAttenuation2',Astop2, ...
          'DesignMethod','butter','SampleRate',sampf);

        filtData=filter(d,ETMXZ_out);
        X=filtData(startTime:endTime);
        filtData=filter(d,ETMYZ_out);
        Y=filtData(startTime:endTime);
        filtData=filter(d,ITMYZ_out);
        C=filtData(startTime:endTime);
    %     injAng=150*pi/180;
    %     X=filter(d,seed(ceil(100*cos(injAng))+150:length(ETMXZ_out)+ceil(100*cos(injAng))+150))';
    %     Y=filter(d,seed(ceil(100*sin(injAng))+150:length(ETMXZ_out)+ceil(100*sin(injAng))+150))';
    %     C=filter(d,seed(1+150:length(ETMXZ_out)+1+150))';
    %     injVel=7000;
    %     X=filter(d,seed(4e3/(sqrt(2)*injVel)*10*sampf+150:length(ETMXZ_out)+4e3/(sqrt(2)*injVel)*10*sampf+150))';
    %     Y=filter(d,seed(4e3/(sqrt(2)*injVel)*10*sampf+150:length(ETMXZ_out)+4e3/(sqrt(2)*injVel)*10*sampf+150))';
    %     C=filter(d,seed(1+150:length(ETMXZ_out)+1+150))';

    %     X=ETMXZ_out(startTime*sampf:length(ETMXZ_out));
    %     Y=ETMYZ_out(startTime*sampf:length(ETMYZ_out));
    %     C=ITMYZ_out(startTime*sampf:length(ITMYZ_out));
%     for j=1:floor(length(X)*(freq2/sampf))-2
    delta_t_X=[];
    delta_t_Y=[];
    len=3000;
    for j=1:floor(length(X)/len)-1
        if (max(C(j*len:(j+1)*len))-min(C(j*len:(j+1)*len)))>=threshold
%             if (max(C(j*len:(j+1)*len))-min(C(j*len:(j+1)*len)))>=max(abs(C))/3
%         if (max(C)-min(C))>=threshold    
%             [crossY,~] = xcorr(C(floor(j/(freq2/sampf)):floor((j+1)/(freq2/sampf))...
%                 ,Y(floor(j/(freq2/sampf)):floor((j+1)/(freq2/sampf)))));
%             [crossX,lags]=xcorr(C(floor(j/(freq2/sampf)):floor((j+1)/(freq2/sampf))...
%                 ,X(floor(j/(freq2/sampf)):floor((j+1)/(freq2/sampf)))));
            [crossY,~] = xcorr(C(j*len:(j+1)*len),Y(j*len:(j+1)*len));
            [crossX,lags]=xcorr(C(j*len:(j+1)*len),X(j*len:(j+1)*len));
%             [crossY,~] = xcorr(C,Y);
%             [crossX,lags]=xcorr(C,X);
            crossX=abs(crossX);
            crossY=abs(crossY);
%             peak=crossX(floor(.95*find(crossX==max(crossX))):floor(1.05*find(crossX==max(crossX))));
%             peakLags=lags(floor(.95*find(crossX==max(crossX))):floor(1.05*find(crossX==max(crossX))))';
            peak=crossX(floor(find(crossX==max(crossX))-1.5/freq2):floor(find(crossX==max(crossX))+1.5/freq2));
            peakLags=lags(floor(find(crossX==max(crossX))-1.5/freq2):floor(find(crossX==max(crossX))+1.5/freq2))';
            [fit,s]=polyfit(peakLags,peak,2);  
            delta_t_X=[delta_t_X -fit(2)/(2*fit(1))/sampf];
%             [ind ind]=max(crossY);
%             delta_t_X=[delta_t_X lags(ind)/sampf];
%             sigmaTX=sum((peak-(fit(1)*peakLags.^2+fit(2)*peakLags+fit(3))).^2);
%             errFit= sqrt(diag(inv(s.R)*inv(s.R'))./s.normr.^2./s.df);
%             sigmaTX=sqrt((-1/(2*fit(1))/sampf)^2*errFit(2)^2+(fit(2)/(2*fit(1)^2)/sampf)^2*errFit(1)^2);

%             peak=crossY(floor(.95*find(crossY==max(crossY))):floor(1.05*find(crossY==max(crossY))));
%             peakLags=lags(floor(.95*find(crossY==max(crossY))):floor(1.05*find(crossY==max(crossY))))'; 
            peak=crossY(floor(find(crossY==max(crossY))-1.5/freq2):floor(find(crossY==max(crossY))+1.5/freq2));
            peakLags=lags(floor(find(crossY==max(crossY))-1.5/freq2):floor(find(crossY==max(crossY))+1.5/freq2))';
            [fit,s]=polyfit(peakLags,peak,2);
            delta_t_Y=[delta_t_Y -fit(2)/(2*fit(1))/sampf];
%             [ind ind]=max(crossY);
%             delta_t_Y=[delta_t_Y lags(ind)/sampf];
%             sigmaTY=sum((peak-(fit(1)*peakLags.^2+fit(2)*peakLags+fit(3))).^2);
%             errFit= sqrt(diag(inv(s.R)*inv(s.R'))./s.normr.^2./s.df);
%             sigmaTY=sqrt((-1/(2*fit(1))/sampf)^2*errFit(2)^2+(fit(2)/(2*fit(1)^2)/sampf)^2*errFit(1)^2);
%         else
%             delta_t_X=[delta_t_X nan];
%             delta_t_Y=[delta_t_Y nan;
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
%     
    bootVel=[bootVel; tempBVel'];
    bootAng=[bootAng; tempBAng'];
%     
%         bootVel=[bootVel; nan];
%         bootAng=[bootAng; nan];
    %     delta_t_X=(lags(find(crossX==max(crossX))))/sampf;
    %     delta_t_Y=(lags(find(crossY==max(crossY))))/sampf;
    %     delta_t_X=(sum(crossX.*lags)/sum(crossX))/sampf;
    %     delta_t_Y=(sum(crossY.*lags)/sum(crossY))/sampf;

%         lagsX=[lagsX delta_t_X];
%         lagsY=[lagsY delta_t_Y];
%         sigmaTX=std(delta_t_X);
%         sigmaTY=std(delta_t_Y);
        sigmaTX=0;
        sigmaTY=0;
        delta_t_X=mean(delta_t_X);
        delta_t_Y=mean(delta_t_Y);
        ang=[ang atan2(delta_t_Y,delta_t_X)*180/pi];    
        vel=[vel 4e3./sqrt((delta_t_X).^2+(delta_t_Y).^2)];
        sigmaAng=[sigmaAng sqrt(delta_t_X.^2.*sigmaTY.^2+delta_t_Y.^2.*sigmaTX.^2)/(delta_t_X.^2+delta_t_Y.^2)*180/pi];
        sigmaVel=[sigmaVel 2*4e3*sqrt(delta_t_X.^2.*sigmaTX.^2+delta_t_Y.^2.*sigmaTY.^2)/(delta_t_X.^2+delta_t_Y.^2).^2];
    end
    % vel=vel(find(abs(vel-mean(vel))<=abs(3*std(vel))));
% 
% ang=ang(find(abs(ang-mean(ang))<=abs(3*std(ang))));

% Navg3 = 1;
% [Ty,~] = tfe2(ETMYZ_out,ITMYZ_out ,1/sampf, Navg3);
% [Tx,F] = tfe2(ETMXZ_out,ITMYZ_out,1/sampf, Navg3);
% 
% F=F(find(F>=.01 & F<=.1));
% Tx=Tx(find(F>=.01 & F<=.1));
% Ty=Ty(find(F>=.01 & F<=.1));
% 
% delta_t_X=angle(Tx)/pi./F';
% delta_t_Y=angle(Ty)/pi./F';
% 
% ang=atan(delta_t_Y./delta_t_X)*180/pi;
% vel=4e3./sqrt((delta_t_X).^2+(delta_t_Y).^2);
% 
% vel=vel(find(isinf(vel)==false));
% 
% F=F(find(abs(vel-mean(vel))<=abs(3*std(vel))));
% ang=ang(find(abs(vel-mean(vel))<=abs(3*std(vel))));
% vel=vel(find(abs(vel-mean(vel))<=abs(3*std(vel))));
end