close all

x=1:4096;

refMu=zeros(1000,1);
mainMu=zeros(1000,1);

startIndexLeft=200;
startIndexRight=2000;
spacing=60;
width=10;
refAmp=2^12*0.85;
mainAmp=2^12*0.6;
gap=30;

patternLength=20*(width+spacing)+gap;

firstRefFrame=synthPattern(x, refAmp, startIndexLeft, width, spacing, gap);
firstMainFrame=synthPattern(x, mainAmp, startIndexRight, width, spacing, gap);


refCor=xcorr(firstRefFrame(startIndexLeft-4*width: startIndexLeft-width+patternLength),...
    firstRefFrame(startIndexLeft-4*width: startIndexLeft-width+patternLength));

mainCor=xcorr(firstMainFrame(startIndexRight-4*width: startIndexRight-width+patternLength),...
    firstMainFrame(startIndexRight-4*width: startIndexRight-width+patternLength));

figure(1)
subplot(3,2,[1,2]);
plot1 = plot(x, awgn(firstRefFrame+firstMainFrame,25,'measured'));
hold on
plot1line1= plot([startIndexLeft-width startIndexLeft-width], [0 2^12]);
plot1line2= plot([startIndexRight-width startIndexRight-width], [0 2^12]);
hold off
xlim([1 4096])
ylim([0 2^12])
xlabel('Pixel Number')
ylabel('Intensity')
legend('Data','Start Index Left','Start Index Right')

fitLength=20;
fitCor=refCor(floor(length(refCor)/2)-fitLength: floor(length(refCor)/2)+fitLength);

fitX1=(-fitLength: fitLength);
fitX2=fitX1.^2;
fit1=ones(length(fitX1),1);
fitX=[fitX2' fitX1' fit1];

w=inv(fitX'*fitX)*fitX'*log(fitCor)';

mu=w(2)/(2*w(1));
sigma=sqrt(1/abs(w(1)));
A=exp(w(3)+mu^2/sigma^2);

subplot(3,2,3); 
plot2=plot((1:length(refCor))-length(refCor)/2+startIndexLeft-4*width, refCor);
hold on
fit2=plot((-fitLength: fitLength)+startIndexLeft-4*width,A*exp(-((-fitLength: fitLength)+mu).^2/sigma^2));
center2= plot([mu+startIndexLeft-4*width mu+startIndexLeft-4*width], [0 5e9]);
hold off
xlim([startIndexLeft-200 startIndexLeft+200])
xlabel('Offsets')
ylabel('Reference Cross Corr')
legend('Data','Fit','Center')

fitLength=20;
fitCor=mainCor(floor(length(mainCor)/2)-fitLength: floor(length(mainCor)/2)+fitLength);

fitX1=(-fitLength: fitLength);
fitX2=fitX1.^2;
fit1=ones(length(fitX1),1);
fitX=[fitX2' fitX1' fit1];

w=inv(fitX'*fitX)*fitX'*log(fitCor)';

mu=w(2)/(2*w(1));
sigma=sqrt(1/abs(w(1)));
A=exp(w(3)+mu^2/sigma^2);

subplot(3,2,4);
plot3=plot((1:length(mainCor))-length(mainCor)/2+startIndexRight-4*width, mainCor);
hold on
fit3=plot((-fitLength: fitLength)+startIndexRight-4*width,A*exp(-((-fitLength: fitLength)+mu).^2/sigma^2));
center3= plot([mu+startIndexLeft-4*width mu+startIndexLeft-4*width], [0 2.5e9]);
hold off
xlim([startIndexRight-200 startIndexRight+200])
xlabel('Offsets')
ylabel('Main Cross Corr')
legend('Data','Fit','Center')

subplot(3,2,[5,6]);
angle1= plot(refMu);
hold on
angle2= plot(mainMu);
angle3= plot(mainMu-refMu);
hold off
xlim([1 1000])
ylim([-200 200])
xlabel('Time')
ylabel('Angle')
legend('Ref','Main','Diff')

refMu=[];
mainMu=[];

t=0;
while true
    
    refSignal=100*sin(2*pi*0.01*t);
    mainSignal=100*sin(2*pi*0.1*t)+refSignal;
    
    refFrame=synthPattern(x, refAmp, startIndexLeft+refSignal, width, spacing, gap);
    mainFrame=synthPattern(x, mainAmp, startIndexRight+mainSignal, width, spacing, gap);
    
    set(plot1,'XData',x)
    set(plot1,'YData',awgn(refFrame+mainFrame,25,'measured'))
    
    set(plot1line1,'XData',[startIndexLeft-4*width+refSignal startIndexLeft-4*width+refSignal])
    set(plot1line1,'YData',[0 2^12])
    
    set(plot1line2,'XData',[startIndexRight-4*width+mainSignal startIndexRight-4*width+mainSignal])
    set(plot1line2,'YData',[0 2^12])
    
    refCor=xcorr(refFrame(startIndexLeft-4*width+floor(refSignal): startIndexLeft-width+patternLength+floor(refSignal)),...
        firstRefFrame(startIndexLeft-4*width: startIndexLeft-width+patternLength));
    
    mainCor=xcorr(mainFrame(startIndexRight-4*width+floor(mainSignal): startIndexRight-width+patternLength+floor(mainSignal)),...
        firstMainFrame(startIndexRight-4*width: startIndexRight-width+patternLength));
    
    set(plot2,'XData',(1:length(refCor))-length(refCor)/2+startIndexLeft-4*width+floor(refSignal))
    set(plot2,'YData',refCor)
    
    fitLength=20;
    fitCor=refCor(floor(length(refCor)/2)-fitLength: floor(length(refCor)/2)+fitLength);

    fitX1=(-fitLength: fitLength);
    fitX2=fitX1.^2;
    fit1=ones(length(fitX1),1);
    fitX=[fitX2' fitX1' fit1];

    w=inv(fitX'*fitX)*fitX'*log(fitCor)';

    mu=w(2)/(2*w(1));
    sigma=sqrt(1/abs(w(1)));
    A=exp(w(3)+mu^2/sigma^2);
    
    
    refMu=[refMu; mu+refSignal];
    
    set(fit2,'XData',(-fitLength: fitLength)+startIndexLeft-4*width+refSignal)
    set(fit2,'YData',A*exp(-((-fitLength: fitLength)+mu).^2/sigma^2))
    
    set(center2,'XData',[mu+startIndexLeft-4*width+refSignal mu+startIndexLeft-4*width+refSignal]);
    set(center2,'YData', [0 5e9]);
    
    set(plot3,'XData',(1:length(mainCor))-length(mainCor)/2+startIndexRight-4*width+floor(mainSignal))
    set(plot3,'YData',mainCor)
    
    fitLength=20;
    fitCor=mainCor(floor(length(mainCor)/2)-fitLength: floor(length(mainCor)/2)+fitLength);

    fitX1=(-fitLength: fitLength);
    fitX2=fitX1.^2;
    fit1=ones(length(fitX1),1);
    fitX=[fitX2' fitX1' fit1];

    w=inv(fitX'*fitX)*fitX'*log(fitCor)';

    mu=w(2)/(2*w(1));
    sigma=sqrt(1/abs(w(1)));
    A=exp(w(3)+mu^2/sigma^2);
    
    mainMu=[mainMu; mu+mainSignal];
    
    set(fit3,'XData',(-fitLength: fitLength)+startIndexRight-4*width+mainSignal)
    set(fit3,'YData',A*exp(-((-fitLength: fitLength)+mu).^2/sigma^2))    
    
    set(center3,'XData',[mu+startIndexRight-4*width+mainSignal mu+startIndexRight-4*width+mainSignal]);
    set(center3,'YData', [0 2.5e9]);
    
    set(angle1,'XData',(1:length(refMu)))
    set(angle1,'YData',refMu)
    
    set(angle2,'XData',(1:length(mainMu)))
    set(angle2,'YData',mainMu)    
    
    set(angle3,'XData',(1:length(mainMu)))
    set(angle3,'YData',mainMu-refMu)
    
    refreshdata
    drawnow
    
    t=t+0.1;
    pause(0.1)
end
% 
% figure(1)
% frames(k)=getframe(gcf);

% v = VideoWriter('dispersionSwarm.avi');
% v.FrameRate=2;
% open(v);
% writeVideo(v,frames);
% close(v);