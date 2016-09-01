function [errFreq,transXErr,transYErr,transZErr,tiltErr]=RWaveMeasErr
sampf=8;
[ETMXZ_out, ITMYZ_out, ETMYX_out, ETMYY_out, ETMYZ_out, BRSY_out]=...
    RWaveDataIn('GPS1149323500_Quite.mat');

[ATransX, ~] = asd2(ETMYX_out,1/sampf, 3, 1, @hann);
[ATransY, ~] = asd2(ETMYY_out,1/sampf, 3, 1, @hann);
[ATransZ, ~] = asd2(ETMYZ_out,1/sampf, 3, 1, @hann);
[ATilt, Af] = asd2(BRSY_out,1/sampf, 3, 1, @hann);

transXErr=smooth(ATransX,10);
transYErr=smooth(ATransY,10);
transZErr=smooth(ATransZ,10);
tiltErr=smooth(ATilt,10);
errFreq=Af;

% loglog(Af,transXErr)
% plot(ETMYX_out)
end