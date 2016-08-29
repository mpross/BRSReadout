function [vel, ang, sigmaVel,sigmaAng]=RWaveArray(ETMXZ_out,ETMYZ_out,ITMYZ_out,sampf)
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
        if max(C>=threshold)
            [crossY,~] = xcorr(C,Y);
            [crossX,lags]=xcorr(C,X);
            crossX=abs(crossX);
            crossY=abs(crossY);
            peak=crossX(floor(.9995*find(crossX==max(crossX))):floor(1.0005*find(crossX==max(crossX))));
            peakLags=lags(floor(.9995*find(crossX==max(crossX))):floor(1.0005*find(crossX==max(crossX))))';
            [fit,s]=polyfit(peakLags,peak,2);  
            delta_t_X=-fit(2)/(2*fit(1))/sampf;        
    %         sigmaTX=std(peak-(fit(1)*peakLags.^2+fit(2)*peakLags+fit(3)));
    %         errFit= sqrt(diag(inv(s.R)*inv(s.R'))./s.normr.^2./s.df);
    %         sigmaTX=sqrt((-1/(2*fit(1))/sampf)^2*errFit(2)^2+(fit(2)/(2*fit(1)^2)/sampf)^2*errFit(1)^2);

            peak=crossY(floor(.9995*find(crossY==max(crossY))):floor(1.0005*find(crossY==max(crossY))));
            peakLags=lags(floor(.9995*find(crossY==max(crossY))):floor(1.0005*find(crossY==max(crossY))))';
            [fit,s]=polyfit(peakLags,peak,2);  
            delta_t_Y=-fit(2)/(2*fit(1))/sampf;       
    %         sigmaTY=std(peak-(fit(1)*peakLags.^2+fit(2)*peakLags+fit(3)));
    %         errFit= sqrt(diag(inv(s.R)*inv(s.R'))./s.normr.^2./s.df);
    %         sigmaTY=sqrt((-1/(2*fit(1))/sampf)^2*errFit(2)^2+(fit(2)/(2*fit(1)^2)/sampf)^2*errFit(1)^2);
        else
            delta_t_X=nan;
            delta_t_Y=nan;
        end
    %     delta_t_X=(lags(find(crossX==max(crossX))))/sampf;
    %     delta_t_Y=(lags(find(crossY==max(crossY))))/sampf;
    %     delta_t_X=(sum(crossX.*lags)/sum(crossX))/sampf;
    %     delta_t_Y=(sum(crossY.*lags)/sum(crossY))/sampf;

        lagsX=[lagsX delta_t_X];
        lagsY=[lagsY delta_t_Y];

        ang=[ang atan2(delta_t_Y,delta_t_X)*180/pi];
        vel=[vel 4e3/sqrt((delta_t_X)^2+(delta_t_Y)^2)];
        sigmaAng=[sigmaAng sqrt(delta_t_X^2*sigmaTY^2+delta_t_Y^2*sigmaTX^2)/(delta_t_X^2+delta_t_Y^2)*180/pi];
        sigmaVel=[sigmaVel 2*4e3*sqrt(delta_t_X^2*sigmaTX^2+delta_t_Y^2*sigmaTY^2)/(delta_t_X^2+delta_t_Y^2)^2];
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