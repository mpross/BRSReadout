function [vel, ang, sigmaVel,sigmaAng,bootVel,bootAng]=RWaveArray(ETMXZ_out,ETMYZ_out,ITMYZ_out,sampf)
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
    startFreq2=0.03;
    freqStep2=.001;
    sigmaTX=0.0231;
    sigmaTY=0.0325;
    sigmaAng=[];
    sigmaVel=[];
    threshold=2e-6;
    startTime=01*sampf;
    % endTime=800*sampf;
    endTime=length(ETMYZ_out);
    delta_t_X=[];
    delta_t_Y=[];
    bootVel=[];
    bootAng=[];
    for i=0:100
        freq2=(startFreq2+i*freqStep2);
    %     [bb,aa] = butter(2,[(0.032+i*0.005)/sampf, (0.032+(i+1)*0.005)/sampf]);
        d = designfilt('bandpassiir', ...
          'StopbandFrequency1',freq2*.8,'PassbandFrequency1', freq2*.9, ...
          'PassbandFrequency2',freq2*1.1,'StopbandFrequency2', freq2*1.2, ...
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
    len=5000;
    for j=1:floor(length(X)/len)-1          
        if max(C>=threshold)          
%             [crossY,~] = xcorr(C(floor(j/(freq2/sampf)):floor((j+1)/(freq2/sampf))...
%                 ,Y(floor(j/(freq2/sampf)):floor((j+1)/(freq2/sampf)))));
%             [crossX,lags]=xcorr(C(floor(j/(freq2/sampf)):floor((j+1)/(freq2/sampf))...
%                 ,X(floor(j/(freq2/sampf)):floor((j+1)/(freq2/sampf)))));
            [crossY,~] = xcorr(C(j*len:(j+1)*len),Y(j*len:(j+1)*len));
            [crossX,lags]=xcorr(C(j*len:(j+1)*len),X(j*len:(j+1)*len));
            crossX=abs(crossX);
            crossY=abs(crossY);
%             peak=crossX(floor(.95*find(crossX==max(crossX))):floor(1.05*find(crossX==max(crossX))));
%             peakLags=lags(floor(.95*find(crossX==max(crossX))):floor(1.05*find(crossX==max(crossX))))';
            peak=crossX(floor(find(crossX==max(crossX))-1.5/freq2):floor(find(crossX==max(crossX))+1.5/freq2));
            peakLags=lags(floor(find(crossX==max(crossX))-1.5/freq2):floor(find(crossX==max(crossX))+1.5/freq2))';
            [fit,s]=polyfit(peakLags,peak,2);  
            delta_t_X=[delta_t_X -fit(2)/(2*fit(1))/sampf];        
    %         sigmaTX=std(peak-(fit(1)*peakLags.^2+fit(2)*peakLags+fit(3)));
    %         errFit= sqrt(diag(inv(s.R)*inv(s.R'))./s.normr.^2./s.df);
    %         sigmaTX=sqrt((-1/(2*fit(1))/sampf)^2*errFit(2)^2+(fit(2)/(2*fit(1)^2)/sampf)^2*errFit(1)^2);

%             peak=crossY(floor(.95*find(crossY==max(crossY))):floor(1.05*find(crossY==max(crossY))));
%             peakLags=lags(floor(.95*find(crossY==max(crossY))):floor(1.05*find(crossY==max(crossY))))'; 
            peak=crossY(floor(find(crossY==max(crossY))-1.5/freq2):floor(find(crossY==max(crossY))+1.5/freq2));
            peakLags=lags(floor(find(crossY==max(crossY))-1.5/freq2):floor(find(crossY==max(crossY))+1.5/freq2))';
            [fit,s]=polyfit(peakLags,peak,2);  
            delta_t_Y=[delta_t_Y -fit(2)/(2*fit(1))/sampf];       
    %         sigmaTY=std(peak-(fit(1)*peakLags.^2+fit(2)*peakLags+fit(3)));
    %         errFit= sqrt(diag(inv(s.R)*inv(s.R'))./s.normr.^2./s.df);
    %         sigmaTY=sqrt((-1/(2*fit(1))/sampf)^2*errFit(2)^2+(fit(2)/(2*fit(1)^2)/sampf)^2*errFit(1)^2);
%         else
%             delta_t_X=[delta_t_X nan];
%             delta_t_Y=[delta_t_Y nan];
        end
    end   
    tempBAng=[];
    tempBVel=[];
    for k=0:10000
        bootTX=bootstrapData(delta_t_X);
        bootTY=bootstrapData(delta_t_Y);
        tempBAng=[tempBAng; mean(atan2(bootTY,bootTX)*180/pi)];
        tempBVel=[tempBVel; mean(4e3./sqrt((bootTX).^2+(bootTY).^2))];
    end
    bootVel=[bootVel; tempBVel'];
    bootAng=[bootAng; tempBAng'];
    %     delta_t_X=(lags(find(crossX==max(crossX))))/sampf;
    %     delta_t_Y=(lags(find(crossY==max(crossY))))/sampf;
    %     delta_t_X=(sum(crossX.*lags)/sum(crossX))/sampf;
    %     delta_t_Y=(sum(crossY.*lags)/sum(crossY))/sampf;

%         lagsX=[lagsX delta_t_X];
%         lagsY=[lagsY delta_t_Y];
        delta_t_X=mean(delta_t_X);
        delta_t_Y=mean(delta_t_Y);
        sigmaTX=std(delta_t_X);
        sigmaTY=std(delta_t_Y);
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