printFigs = false;

fileload = true;

sampf = 8;

% if (exist('GPS1166005817_10k_EQ7_PapNG','file')&& fileload==true)

        myfile = load('GPS1166005817_10k_EQ7_PapNG.mat');

        mydata = myfile.rawdata8Hz1;
%         mydata=myfile.mydata;
% 
        rawBRSY= mydata(:,4);

        rawETMXZ = mydata(:,2);

        rawETMYZ = mydata(:,1);

        rawITMYZ = mydata(:,3);

        rawETMYY = mydata(:,5);
% 
%         rawBRSY= mydata(:,10);
% 
%         rawETMXZ = mydata(:,3);
% 
%         rawETMYZ = mydata(:,6);
% 
%         rawITMYZ = mydata(:,9);
% 
%         rawETMYY = mydata(:,5);
        
%         rawBRSY= mydata(:,4);
% 
%         rawETMXZ = mydata(:,3);
% 
%         rawETMYZ = mydata(:,2);
% 
%         rawITMYZ = mydata(:,1);
% % 
%         rawETMYY = mydata(:,5);

% end



sttime =000*sampf + 1;

edtime = sttime + 2000*sampf;

% edtime = length(rawBRSY);



BRSY = 1e-9*rawBRSY(sttime:edtime);

ETMYZ = 1e-9*rawETMYZ(sttime:edtime);

ETMYY = 1e-9*rawETMYY(sttime:edtime);

eqtime = [1:length(BRSY)]*1/8-1/8;



%% Filters

[bb,aa] = butter(4,[2*0.001/sampf, 2*3/sampf]);

    % %T240 response inversion filter

STSInvertFilt = zpk(-2*pi*[pairQ(8.2e-3,0.7)],-2*pi*[0 0],1);

STSInvertFilt = 1*STSInvertFilt/abs(freqresp(STSInvertFilt,2*pi*100));



BRSYInvertFilt = zpk(-2*pi*[pairQ(7.73e-3,3000)],-2*pi*[0 0],1);

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

[ABRSY,~]=asd2(BRSY_out,1/sampf,Navg,1,@hann);



% [FCohYY,freqs2]=coh2(BRSY,ETMYY,1/sampf,Navg2,1,@hann);

% [FCohYZ,~]=coh2(BRSY,ETMYZ,1/sampf,Navg2,1,@hann);



%% Plots

figure(4)
plot(BRSY_out)

figure(1)
    subplot('Position',[0.1, 0.7, 0.8, 0.25]) 

    l = plot(eqtime(sampf*250:end)-250, EYZ_vel_filt(sampf*250:end)/1e-6);

    legend('Z velocity','Location','SouthEast');
    
    title('M7.9 Papua New Guinea')

    set(l,'LineWidth',1);

    set(gca,'FontSize',15)

    %xlabel('time (s)');

    xlim([0 4750])

    ylim([-175 175])

    set(gca,'xtick',500*[0:10])

    set(gca,'ytick',50*[-40:40])

    set(gca,'xticklabel',[])

    ylabel('Velocity (\mum/s)');

    grid on



   subplot('Position',[0.1, 0.4, 0.8, 0.25])

    ll = plot(eqtime(sampf*250:end)-250, BRSY_out(sampf*250:end)/1e-9);

    legend('BRS output','Location','SouthEast');

    set(ll,'LineWidth',1);

    set(gca,'FontSize',15)

    xlim([0 4750])

    ylim([-40 40])

    xlabel('Time (s)');

    set(gca,'xtick',500*[0:10])

    set(gca,'ytick',10*[-40:40])
    
    set(gca,'xticklabel',[])

    ylabel('Angle (nrad)');

    grid on

    

    subplot('Position',[0.1, 0.1, 0.8, 0.25])

    lll = plot(eqtime(sampf*250:end)-250, EYY_disp_filt(sampf*250:end)/1e-6);

    legend('Y displacement','Location','SouthEast');

    set(lll,'LineWidth',1);

    set(gca,'FontSize',15)

    xlim([0 4750])

    set(gca,'xtick',500*[0:10])

    xlabel('Time (s)');

    ylabel('Displacement (\mum)');

    ylim([-700 700])

    set(gca,'ytick',200*[-3:3])

    grid on



figure(2)

    ll = loglog(freqs,AEYZ,freqs,AEYY,'--',freqs, ABRSY);

    set(ll,'LineWidth',2);

    xlim([5e-3 5])

    ylim([1e-10 1e-5])

    set(gca,'YTick',10.^(-10:-5))

    %set(gca,'XTick',10.^(-2:2))

    set(gca,'FontSize',15)

    legend('Z velocity','Y velocity','BRS output','BRS Y calibrated','Location','Northeast')

    xlabel('Frequency (Hz)')

    ylabel('ASD (rad/$\sqrt{Hz}$ or m/s/$\sqrt{Hz}$)','Interpreter','Latex')

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