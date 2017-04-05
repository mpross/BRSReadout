printFigs = false;

fileload = true;

sampf = 8;

if (exist('GPS1149323500_Quite','file')&& fileload==true)

        myfile = load('GPS1149323500_Quite.mat');

         mydata = myfile.rawdata8Hz1;

        rawBRSY= mydata(:,4);

        rawETMXZ = mydata(:,2);

        rawETMYZ = mydata(:,1);

        rawITMYZ = mydata(:,3);

        rawETMYY = mydata(:,5);

end



sttime = 000*sampf + 1;

edtime = sttime + 2000*sampf;

%edtime = length(rawBRSY);



BRSY = 1e-9*rawBRSY(sttime:edtime);

ETMYZ = 1e-9*rawETMYZ(sttime:edtime);

ETMYY = 1e-9*rawETMYY(sttime:edtime);

eqtime = [1:length(BRSY)]*1/8-1/8;



%% Filters

[bb,aa] = butter(4,[2*0.03/sampf, 2*.3/sampf]);

    % %T240 response inversion filter

STSInvertFilt = zpk(-2*pi*[pairQ(8.2e-3,0.7)],-2*pi*[0 0],1);

STSInvertFilt = 1*STSInvertFilt/abs(freqresp(STSInvertFilt,2*pi*100));



BRSYInvertFilt = zpk(-2*pi*[pairQ(7.74e-3,3000)],-2*pi*[0 0],1);

BRSYInvertFilt = 1*BRSYInvertFilt/abs(freqresp(BRSYInvertFilt,2*pi*100));

% 

% %Filters to differentiate and integrate

DiffFilt = zpk(-2*pi*[0], -2*pi*2,1);

DiffFilt = 1*DiffFilt/abs(freqresp(DiffFilt,2*pi*0.1592));

% 

IntFilt = zpk(-2*pi*[], -2*pi*5e-4,1);

IntFilt = 1*IntFilt/abs(freqresp(IntFilt,2*pi*0.1592));



% Apply filters

BRSYcal_out = lsim(BRSYInvertFilt,BRSY, eqtime);

BRSY_out=filter(bb,aa,BRSYcal_out);



EYZcal_vel = lsim(STSInvertFilt,ETMYZ,eqtime);

EYZcal_disp = lsim(IntFilt,EYZcal_vel,eqtime);

EYZ_disp_filt=filter(bb,aa,EYZcal_disp);

EYZ_vel_filt=filter(bb,aa,EYZcal_vel);



EYYcal_vel = lsim(STSInvertFilt,ETMYY,eqtime);

EYYcal_disp = lsim(IntFilt,EYYcal_vel,eqtime);

EYY_disp_filt=filter(bb,aa,EYYcal_disp);

EYY_vel_filt=filter(bb,aa,EYYcal_vel);



%% ASD

Navg=5;

Navg2=9;



[AEYY,freqs]=asd2(ETMYY,1/sampf,Navg,1,@hann);

[AEYZ,~]=asd2(ETMYZ,1/sampf,Navg,1,@hann);

[ABRSY,~]=asd2(BRSY,1/sampf,Navg,1,@hann);

[ABRSYcal,~]=asd2(BRSYcal_out,1/sampf,Navg,1,@hann);



% [FCohYY,freqs2]=coh2(BRSY,ETMYY,1/sampf,Navg2,1,@hann);

% [FCohYZ,~]=coh2(BRSY,ETMYZ,1/sampf,Navg2,1,@hann);



%% Plots



figure(1)

    subplot('Position',[0.1, 0.59, 0.8, 0.35])

    l = plot(eqtime(sampf*250:end)-250, EYZ_vel_filt(sampf*250:end));

    legend('Z velocity','Location','SouthEast');

    set(l,'LineWidth',1);

    set(gca,'FontSize',15)

    %xlabel('time (s)');

%     xlim([0 4000])

    ylim([-5e-5 5e-5])

    set(gca,'xtick',500*[0:10])

    set(gca,'ytick',2e-5*[-3:3])

    set(gca,'xticklabel',[])

    ylabel('Velocity (m/s)');

    grid on



   subplot('Position',[0.1, 0.15, 0.8, 0.35])

    ll = plot(eqtime(sampf*250:end)-250, BRSY_out(sampf*250:end));

    legend('BRS output','Location','SouthEast');

    set(ll,'LineWidth',1);

    set(gca,'FontSize',15)

%     xlim([0 4000])

    ylim([-2e-8 2e-8])

    %xlabel('time (s)');

    set(gca,'xtick',500*[0:10])

    set(gca,'ytick',1e-8*[-2:2])

    ylabel('angle (rad)');

    grid on

    

%     subplot('Position',[0.1, 0.1, 0.8, 0.25])
% 
%     lll = plot(eqtime(sampf*250:end)-250, EYY_disp_filt(sampf*250:end));
% 
%     legend('Y displacement','Location','SouthEast');
% 
%     set(lll,'LineWidth',1);
% 
%     set(gca,'FontSize',15)
% 
%     xlim([0 4800])
% 
%     set(gca,'xtick',500*[0:10])
% 
%     xlabel('time (s)');
% 
%     ylabel('displacement (m)');
% 
%     ylim([-7e-4 7e-4])
% 
%     set(gca,'ytick',2e-4*[-3:3])
% 
%     grid on



figure(2)

    ll = loglog(freqs,AEYZ,freqs,AEYY,'--',freqs, ABRSY,':', freqs, ABRSYcal,'-.');

    set(ll,'LineWidth',2);

    xlim([5e-3 5])

    ylim([1e-10 1e-3])

    set(gca,'YTick',10.^(-10:-3))

    %set(gca,'XTick',10.^(-2:2))

    set(gca,'FontSize',15)

    legend('STS Z velocity','STS Y velocity','BRS Y raw','BRS Y calibrated','Location','Northeast')

    xlabel('frequency (Hz)')

    ylabel('ASD (rad/\surd(Hz) or m/s/\surd(Hz))')

    grid on

    

% figure(3)
% 
%     lll = semilogx(freqs2,FCohYY,freqs2,FCohYZ);
% 
%     set(lll,'LineWidth',1);
% 
%     set(gca,'FontSize',15)
% 
%     xlim([5e-3 10])
% 
%     xlabel('frequency (Hz)')
% 
%     ylabel('Coherence')
% 
%     legend('Between BRS and EY Y','BRS Y and EY Z')
% 
%     grid on
% 
% 

if printFigs

%     figure(1)

%     FillPage('w')

%     saveas(gcf,['BRS_PapNG_EQ_Dec17_2016_TT.pdf'])



    figure(2)

    FillPage('w')

    saveas(gcf,['BRS_PapNG_EQ_Dec17_2016_ASD_Quiet.pdf'])

    

end