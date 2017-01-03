function [ETMXZ_out, ITMYZ_out, ETMYX_out, ETMYY_out, ETMYZ_out, BRSY_out]=RWaveDataIn(fileName,new)   
%     close all
    sampf = 8; % sampling frequency in Hz

    % % 
    % if exist('GPS1149581095_5_2Earthquake.mat','file')
    %     myfile = load('GPS1149581095_5_2Earthquake.mat');
%     if exist('GPS1143962787_6_9Earthquake.mat','file')
%         myfile = load('GPS1143962787_6_9Earthquake.mat');
    % if exist('GPS1149331885_6_2Earthquake.mat','file')
    %     myfile = load('GPS1149331885_6_2Earthquake.mat');
    % if exist('GPS1144888165_7_8Earthquake.mat','file')
    %     myfile = load('GPS1144888165_7_8Earthquake.mat');
    if (exist(fileName,'file')&& new==false)
        myfile = load(fileName);
        mydata = myfile.mydata;
        rawETMXX= mydata(:,1);
        rawETMXY = mydata(:,2);
        rawETMXZ = mydata(:,3);
        rawETMYX = mydata(:,4);
        rawETMYY = mydata(:,5);
        rawETMYZ = mydata(:,6);
        rawITMYX = mydata(:,7);
        rawITMYY = mydata(:,8);
        rawITMYZ = mydata(:,9);
        rawBRSY= mydata(:,10);
        rawBRSX=mydata(:,11);
    %     rawWINDX=mydata(:,12);
    %     rawWINDY=mydata(:,13);
    end   
    if (exist(fileName,'file')&& new==true)
        myfile = load(fileName);
        mydata = myfile.rawdata8Hz1;
        rawBRSY= mydata(:,4);
        rawETMXX= zeros([1 length(rawBRSY)]);
        rawETMXY =  zeros([1 length(rawBRSY)]);
        rawETMXZ = mydata(:,2);
        rawETMYX = zeros([1 length(rawBRSY)]);
        rawETMYY = zeros([1 length(rawBRSY)]);
        rawETMYZ = mydata(:,1);
        rawITMYX =  zeros([1 length(rawBRSY)]);
        rawITMYY =  zeros([1 length(rawBRSY)]);
        rawITMYZ = mydata(:,3);        
        rawBRSX= zeros([1 length(rawBRSY)]);
    %     rawWINDX=mydata(:,12);
    %     rawWINDY=mydata(:,13);
    end 
    Sttime =01;
    Endtime=length(rawBRSY);
    localg = 9.8;
    % BRSscale=.85;
    BRSscale=1;

    ETMXX=1e-9 *rawETMXX(Sttime:Endtime);
    ETMXY=1e-9 *rawETMXY(Sttime:Endtime);
    ETMXZ=1e-9 *rawETMXZ(Sttime:Endtime);
    ETMYX=1e-9 *rawETMYX(Sttime:Endtime);
    % ETMYX=seed*cos(-130/180*pi);
    ETMYY=1e-9 *rawETMYY(Sttime:Endtime);
    % ETMYY=seed*sin(-130/180*pi);
    ETMYZ=1e-9 *rawETMYZ(Sttime:Endtime);
    ITMYX=1e-9 *rawITMYX(Sttime:Endtime);
    ITMYY=1e-9 *rawITMYY(Sttime:Endtime);
    ITMYZ=1e-9 *rawITMYZ(Sttime:Endtime);
    BRSY=1e-9*rawBRSY(Sttime:Endtime);
    BRSX=2.13e-9*rawBRSX(Sttime:Endtime);
    Navg =9;
    seed=randn(1,length(ETMYZ));
    %% BRS parameters
    Mtotal = 4.5; % Total mass in kg
    Ibar = 0.61; % Moment of Inertia in kg m^2
    localg = 9.8; % local gravitation acceleration
    doffsetY = .5e-6; % CoM offset from pivot in m. d is plus if CoM below pivot
    doffsetX = 30e-6; 
    resfreq = 2*pi*7.7e-3;
    Qbar = 3000;
    Rot_time = transpose(1/sampf * (0:1:length(BRSY)-1));

    %% torque computation
    % [bb,aa] = butter(4,[2*0.02/sampf 2*0.100/sampf],'bandpass');
    [bb,aa] = butter(4,[2*0.01/sampf, 2*.3/sampf]);
    
    % %T240 response inversion filter
    T240InvertFilt = zpk(-2*pi*[pairQ(8.2e-3,0.7)],-2*pi*[0 0],1);
    T240InvertFilt = 1*T240InvertFilt/abs(freqresp(T240InvertFilt,2*pi*100));
    % 
    % %BRS response inversion filter
    
    BRSXInvertFilt = zpk(-2*pi*[pairQ(8.9e-3,4000)],-2*pi*[0 0],1);
    BRSXInvertFilt = 1*BRSXInvertFilt/abs(freqresp(BRSXInvertFilt,2*pi*100));

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
    T240cal_vel = lsim(T240InvertFilt,ETMXZ,Rot_time);
    T240cal_disp = lsim(IntFilt,T240cal_vel,Rot_time);
    ETMXZ_out=filter(bb,aa,T240cal_disp);

    % 
    T240cal_accY = lsim(DiffFilt,ETMYY,Rot_time);
    %
    T240cal_vel = lsim(T240InvertFilt,ETMYZ,Rot_time);
    T240cal_disp = lsim(IntFilt,T240cal_vel,Rot_time);
    ETMYZ_out=filter(bb,aa,T240cal_disp);

    % 
    T240cal_vel = lsim(T240InvertFilt,ETMYX,Rot_time);
    % T240cal_vel = lsim(T240InvertFilt,seed*cos(130/180*pi),Rot_time);
    T240cal_disp = lsim(IntFilt,T240cal_vel,Rot_time);
    ETMYX_out=filter(bb,aa,T240cal_disp);

    T240cal_vel = lsim(T240InvertFilt,ETMYY,Rot_time);
    % T240cal_vel = lsim(T240InvertFilt,seed*sin(130/180*pi),Rot_time);
    T240cal_disp = lsim(IntFilt,T240cal_vel,Rot_time);
    ETMYY_out=filter(bb,aa,T240cal_disp);

    T240cal_vel = lsim(T240InvertFilt,ETMXZ,Rot_time);
    T240cal_disp = lsim(IntFilt,T240cal_vel,Rot_time);
    ETMXZ_out=filter(bb,aa,T240cal_disp);

    T240cal_vel = lsim(T240InvertFilt,ITMYZ,Rot_time);
    T240cal_disp = lsim(IntFilt,T240cal_vel,Rot_time);
    ITMYZ_out=filter(bb,aa,T240cal_disp);

    BRSYcal_out = lsim(BRSYInvertFilt,BRSY, Rot_time);
    BRSY_out=filter(bb,aa,BRSYcal_out);

%     
%     [AETMYZ, ~] = asd2(ETMYZ_out,1/sampf, Navg, 1, @hann);
% %     [AT240Xaccel_cal, ~] = asd2(T240cal_accX,1/sampf, Navg, 1, @hann);
%     [AT240Yaccel_cal, ~] = asd2(T240cal_accY,1/sampf, Navg, 1, @hann);
%     [ABRSY, Af] = asd2(BRSY_out,1/sampf, Navg, 1, @hann);
    tim=(1:length(ETMYZ_out))/8;
%     figure(1)
%     ll=loglog(Af,AETMYZ,Af,ABRSY,Af,AT240Yaccel_cal.*Mtotal*localg/Ibar*abs(doffsetY)./(2*pi*Af').^2);
%     set(ll,'LineWidth',2);
%     set(gca,'FontSize',18)
%     xlim([3e-3 4])
%     ylim([1e-11 1e-2])
%     set(gca,'YTick',10.^(-10:-2))
%     xlabel('frequency (Hz)');
%     ylabel('ASD (m/s/rt(Hz) or rad/rt(Hz))');
%     title('STS');
%     legend('STS_Z','BRS','BRS Possible Acceleration Coupling')
%     grid on

    figure(2)
    lll=plot(tim-1000,1e6*ETMYZ_out+1000,tim-1000,BRSY_out*1e10);
    hold off
    set(lll,'LineWidth',.5);
    set(gca,'FontSize',18)
    set(gca,'YTick',250.*(-10:20))
%     xlim([0 max(tim)-3000])
%     ylim([-3e-5 1e-4])
    legend('STS_Z ETMY','BRS ETMY')
    ylabel('Displacement (um or 0.1 nrad)')
    xlabel('Time (s)')
    grid on

end