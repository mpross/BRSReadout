close all

avgVel=containers.Map('ValueType','any','KeyType','double');
avgV=containers.Map('ValueType','any','KeyType','double');
avgVelErr=containers.Map('ValueType','any','KeyType','double');
avgVErr=containers.Map('ValueType','any','KeyType','double');
fileName={'GPS1143962787_6_9Earthquake.mat','GPS1144888165_7_8Earthquake.mat','GPS1149581095_5_2Earthquake.mat',...
    'GPS1149331885_6_2Earthquake.mat','GPS1156069817_Burma.mat','GPS1156480217_Atlantic.mat','GPS1156782617_NewZealand.mat',...
    'GPS1157149817_Russia.mat','GPS1163070017_7k_EQ5.mat','GPS1163797217_7k_EQ6.mat','GPS1166005817_10k_EQ7_PapNG.mat','GPS1155000017_7k_EQ8_NewC.mat'};
newArray=[false, false, false, false, true, true, true, true, true, true, true, true,true];
sampf =8;
threshold=2;

for j=[0 5 8 10 11]
    Astop1 = 5;
    Apass  = .5;
    Astop2 =5;
    velA=[];
    ang=[];
    angTim=[];
    velTim=[];
    lagsX=[];
    lagsY=[];      
    delta_t_X=[];
    delta_t_Y=[];
    bootVelA=[];
    bootAng=[];    
    bootVelS=[];
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
    velS=[];    
   
    startFreq=0.025;
    freqStep=.005;
    iter=floor((.1-startFreq)/freqStep);
    sampf = 8; 
    
    %% Data Reading
    if (exist(fileName{j+1},'file')&& newArray(j+1)==false)
        myfile = load(fileName{j+1});
        mydata = myfile.mydata;
        rawETMXZ = mydata(:,3);
        rawETMYZ = mydata(:,6);
        rawITMYZ = mydata(:,9);
        rawBRSY= mydata(:,10);
    end   
    if (exist(fileName{j+1},'file')&& newArray(j+1)==true)
        myfile = load(fileName{j+1});
        mydata = myfile.rawdata8Hz1;
        rawBRSY= mydata(:,4);        
        rawETMXZ = mydata(:,2);
        rawETMYZ = mydata(:,1);
        rawITMYZ = mydata(:,3);       
    end 
    Sttime =01;
    Endtime=length(rawBRSY);
    localg = 9.8;
    BRSscale=1;

    ETMXZ=1e-9 *rawETMXZ(Sttime:Endtime);
    ETMYZ=1e-9 *rawETMYZ(Sttime:Endtime);
    ITMYZ=1e-9 *rawITMYZ(Sttime:Endtime);
    BRSY=1e-9*rawBRSY(Sttime:Endtime);
    Navg =9;
    %% BRS parameters
    Mtotal = 4.5; % Total mass in kg
    Ibar = 0.61; % Moment of Inertia in kg m^2
    localg = 9.8; % local gravitation acceleration
    doffsetY = .5e-6; % CoM offset from pivot in m. d is plus if CoM below pivot
    doffsetX = 30e-6; 
    resfreq = 2*pi*7.7e-3;
    Qbar = 3000;
    Rot_time = transpose(1/sampf * (0:1:length(BRSY)-1));

    %% Response inversion
    % [bb,aa] = butter(4,[2*0.02/sampf 2*0.100/sampf],'bandpass');
    [bb,aa] = butter(4,[2*0.001/sampf, 2*1/sampf]);
    
    % %T240 response inversion filter
    T240InvertFilt = zpk(-2*pi*[pairQ(8.2e-3,0.7)],-2*pi*[0 0],1);
    T240InvertFilt = 1*T240InvertFilt/abs(freqresp(T240InvertFilt,2*pi*100));
    % 
    % %BRS response inversion filter
    BRSYInvertFilt = zpk(-2*pi*[pairQ(7.74e-3,3000)],-2*pi*[0 0],1);
    BRSYInvertFilt = 1*BRSYInvertFilt/abs(freqresp(BRSYInvertFilt,2*pi*100));
    % 
    % %Filters to differentiate and integrate
    DiffFilt = zpk(-2*pi*[0], -2*pi*2,1);
    DiffFilt = 1*DiffFilt/abs(freqresp(DiffFilt,2*pi*0.1592));
    % 
    IntFilt = zpk(-2*pi*[], -2*pi*5e-4,1);
    IntFilt = 1*IntFilt/abs(freqresp(IntFilt,2*pi*0.1592));
    % 
    % %Apply filters
    %
    T240cal_vel = lsim(T240InvertFilt,ETMYZ,Rot_time);
    T240cal_disp = lsim(IntFilt,T240cal_vel,Rot_time);
%     ETMYZ_out=filter(bb,aa,T240cal_disp);
    ETMYZ_out=filter(bb,aa,T240cal_vel);

    T240cal_vel = lsim(T240InvertFilt,ETMXZ,Rot_time);
    T240cal_disp = lsim(IntFilt,T240cal_vel,Rot_time);
%     ETMXZ_out=filter(bb,aa,T240cal_disp);
    ETMXZ_out=filter(bb,aa,T240cal_vel);

    T240cal_vel = lsim(T240InvertFilt,ITMYZ,Rot_time);
    T240cal_disp = lsim(IntFilt,T240cal_vel,Rot_time);
%     ITMYZ_out=filter(bb,aa,T240cal_disp);
    ITMYZ_out=filter(bb,aa,T240cal_vel);

    BRSYcal_out = lsim(BRSYInvertFilt,BRSY, Rot_time);
    BRSY_out=filter(bb,aa,BRSYcal_out);

    tim=(1:length(ETMYZ_out))/8;
    
    ETMXZ_out=ETMXZ_out(300*sampf:length(ETMXZ_out));
    ETMYZ_out=ETMYZ_out(300*sampf:length(ETMYZ_out));
    ITMYZ_out=ITMYZ_out(300*sampf:length(ITMYZ_out));
    BRSY_out=BRSY_out(300*sampf:length(BRSY_out));
    bootVelA=[];
    bootAng=[];
    startTime=1;
    endTime=length(ETMXZ_out);
    threshold=rms(ETMYZ_out(startTime:endTime))*0;
    
  
    for i=0:iter
        %% Data Crunching  
        freq=(startFreq+i*freqStep);
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
        filtData=filter(d,ETMYZ_out-mean(ETMYZ_out));
        seriesZ=Y;
        localThreshold1=max(abs(filtData))*.7*1e9*0;
        filtData=filter(d,BRSY_out-mean(BRSY_out)); 
        seriesRX=filtData(startTime:endTime);    

        FitData=zeros(length(seriesZ),1);
        tempX=[];
        tempY=[];
        tempZ=[];
        tempRX=[];
        r2z=[];
        r2rx=[];        
        fitLength=floor(1/(freq/sampf)/2);
        for j=1:floor(length(seriesRX)/fitLength)-2
           tim=(j*fitLength:(j+1)*fitLength)'./8;
           cut=1e9*seriesZ(j*fitLength:(j+1)*fitLength);
           g = fittype( @(a,b,cen_fr,x) a*sin(2*pi*cen_fr*x)+b*cos(2*pi*cen_fr*x), 'problem', 'cen_fr' );
           [myfit,st] = fit(tim,cut, g,'problem',freq,'StartPoint', [1, 1]);
           r2z=[r2z;st.rsquare];
           a1=coeffvalues(myfit);
           if i==7
                FitData(j*fitLength:(j+1)*fitLength) = a1(1)*sin(2*pi*freq.*tim)+a1(2)*cos(2*pi*freq.*tim);
           end
           cut=1e12*seriesRX(j*fitLength:(j+1)*fitLength);
           g = fittype( @(a,b,cen_fr,x) a*sin(2*pi*cen_fr*x)+b*cos(2*pi*cen_fr*x), 'problem', 'cen_fr' );
           [myfit,st] = fit(tim,cut, g,'problem',freq,'StartPoint', [1, 1]);
           r2rx=[r2rx;st.rsquare];
           a2=coeffvalues(myfit);
           j
           if (abs(r2rx)>0) 
               if (abs(r2z)>0)
                   if(abs(cos(angle(a1(2)+i*a1(1))-angle((a2(2)+i*a2(1)))))>cos(10*pi/180))
                       tempZ=[tempZ a1(2)+i*a1(1)];
                       tempRX=[tempRX (a2(2)+i*a2(1))/1e3];
                   end
               end
           end
        end

        if i==7
            figure(8)
            line1=plot((1:length(FitData)),seriesZ*1e9,(1:length(FitData)),FitData,'--');
            ylabel('Vertical Displacement (um)')
            xlabel('Time (s)')
            legend('Data','Fit')
            set(line1,'LineWidth',1.2)
            set(gca,'FontSize',12)
            grid on
            length(seriesZ)
            
            figure(12)            
            line2=plot((startTime:endTime)/sampf,X,(startTime:endTime)/sampf,Y+1e-6,(startTime:endTime)/sampf,C+2e-6);
            legend('X','Y','C')
            set(line2,'LineWidth',1.2)
            set(gca,'FontSize',12)
            grid on
        end

        delta_t_X=[];
        delta_t_Y=[];
        len=250*sampf;
        for j=1:floor(length(X)/len)-1
            j
            if(max(abs(tempRX(floor(j*len/fitLength):floor((j+1)*len/fitLength))))>=threshold)

                [crossY,~] = xcorr(C(j*len:(j+1)*len),Y(j*len:(j+1)*len));
                [crossX,lags]=xcorr(C(j*len:(j+1)*len),X(j*len:(j+1)*len));

                crossX=abs(crossX);
                crossY=abs(crossY);

                peak=crossX(floor(find(crossX==max(crossX))-1/freq):floor(find(crossX==max(crossX))+1/freq));
                peakLags=lags(floor(find(crossX==max(crossX))-1/freq):floor(find(crossX==max(crossX))+1/freq))';
                [myfit,s]=polyfit(peakLags,peak,2);  
                delta_t_X=[delta_t_X -myfit(2)/(2*myfit(1))/sampf];

                if i==7
                    figure(18)
                    line3=plot(peakLags,peak,peakLags,myfit(1)*peakLags.^2+myfit(2)*peakLags+myfit(3));
                    legend('Cross Cor','Fit')
                    set(line3,'LineWidth',1.2)
                    set(gca,'FontSize',12)
                    grid on
                end

                peak=crossY(floor(find(crossY==max(crossY))-1/freq):floor(find(crossY==max(crossY))+1/freq));
                peakLags=lags(floor(find(crossY==max(crossY))-1/freq):floor(find(crossY==max(crossY))+1/freq))';
                [myfit,s]=polyfit(peakLags,peak,2);
                delta_t_Y=[delta_t_Y -myfit(2)/(2*myfit(1))/sampf];           

            end
        end
        
        %% Bootstrapping
        tempBAng=[];
        tempBVelA=[];
        tempBVelS=[];
        for k=0:1e2
            k
            array=[delta_t_X; delta_t_Y];
            bootArray=bootstrapData(array');
            if length(bootArray)>0
                bootTX=bootArray(:,1)';
                bootTY=bootArray(:,2)';   
            else
                bootTX=nan;
                bootTY=nan;
            end
            tempBAng=[tempBAng; mean(atan2(bootTY,bootTX)*180/pi)];      
            tempBVelA=[tempBVelA; mean(4e3./sqrt((bootTX).^2+(bootTY).^2))];
            
            N=0;
            sumZ=0;
            sumRX=0;
            avgVS=0;
            array=[tempZ; tempRX]; 
            bootArray=bootstrapData(array');
            btempZ=bootArray(:,1)';
            btempRX=bootArray(:,2)';
            for p=1:min([length(btempZ) length(btempRX)])
                if(abs(btempZ(p))>=threshold)
                   avgVS=avgVS+abs(btempZ(p))./abs(btempRX(p)).*sin(tempBAng)*pi/180;                    
                   N=N+1;
                end
            end   
            avgVS=avgVS/N;
            tempBVelS=[tempBVelS; avgVS];
        end
        bootVelS=[bootVelS; tempBVelS'];
        bootVelA=[bootVelA; tempBVelA'];
        bootAng=[bootAng; tempBAng'];
        
        %% Calculations
        tempAng=atan2(delta_t_Y,delta_t_X)*180/pi;
        tempVel=4e3./sqrt((delta_t_X).^2+(delta_t_Y).^2);
        {mean(tempAng),mean(tempVel)}
        ang=[ang mean(tempAng)];
        velA=[velA mean(tempVel)];  
        
        N=0;
        avgPhi=0;
        avgK=0;
        avgEl=0;
        avgVS=0;
        pltV=[];
        for p=1:min([length(tempZ) length(tempRX)])
            if(abs(tempRX(p))>=threshold)
                if(not(isnan(abs(tempZ(p))./abs(tempRX(p)).*sin(ang(i+1)*pi/180))))
                    avgVS=avgVS+abs(tempZ(p))./abs(tempRX(p)).*sin(ang(i+1)*pi/180);
                    pltV=[pltV abs(tempZ(p))./abs(tempRX(p)).*sin(ang(i+1)*pi/180)];
                    N=N+1;
                end
            end
        end   

        avgVS=avgVS/N
        velS=[velS; avgVS];
        r2zavg=[mean(r2z);r2zavg];
        r2rxavg=[mean(r2rx);r2rxavg];
               
        figure(4)
        line4=plot((1:length(delta_t_X))*250*sampf,4e3./sqrt((delta_t_X).^2+(delta_t_Y).^2),(1:length(pltV))*fitLength,pltV);
        set(line4,'LineWidth',1.2)
        set(gca,'FontSize',12)
        grid on        
        
    end   

 %% Single Event Plotting


    cInd5=find(std(bootVelA')>0);
    cInd6=find(std(bootVelS')>0);
    cInd=intersect(cInd5,cInd6);
    velS=velS(cInd);
    velA=velA(cInd);
    bootVelS=bootVelS(cInd,:);
    bootVelA=bootVelA(cInd,:);
     
    singVel=[singVel; mean(velS)];
    singAng=[singAng; mean(ang)]

    figure(5)
    hold on
    errorbar(((cInd-1)*freqStep+startFreq),abs(velS),-std(bootVelS'),std(bootVelS'));
    errorbar(((cInd-1)*freqStep+startFreq),velA,-std(bootVelA'),std(bootVelA'),'--');
    ylabel('Velocity (m/s)')
    xlabel('Frequency (Hz)')
    legend('Single Station-6.7 Vanuatu', 'Array-6.7 Vanuatu','Single Station-7.1 Pacific', 'Array-7.1 Pacific',...
        'Single Station-7.8 New Zealand', 'Array-7.8 New Zealand','Single Station-7.9 Papa New Guinea', 'Array-7.9 Papa New Guinea',...
        'Single Station-7.2 New Caledonia', 'Array-7.2 New Caledonia')
    grid on
    box on
    hold off

    com=intersect(cInd,cell2mat(avgV.keys));
    for k=1:length(com)
        avgV(com(k))=[avgV(com(k)) velS(k)];
        avgVel(com(k))=[avgVel(com(k)) velA(k)];
        avgVErr(com(k))=[avgVErr(com(k)) std(bootVelS(k,:)')^2];
        avgVelErr(com(k))=[avgVelErr(com(k)) std(bootVelA(k,:)')^2];
    end
    new=setxor(cInd,intersect(cInd,cell2mat(avgV.keys)));
    for k=1:length(new)
        avgV(new(k))=velS(k);
        avgVel(new(k))=velA(k);
        avgVErr(new(k))=std(bootVelS(k,:)')^2;
        avgVelErr(new(k))=std(bootVelA(k,:)')^2;
    end
        
end

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

