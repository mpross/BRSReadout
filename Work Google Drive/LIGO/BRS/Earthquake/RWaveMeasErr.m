function [errFreq,transErr,tiltErr]=RWaveMeasErr(trans,tilt,sampf)

[ATilt, ~] = asd2(tilt,1/sampf, 3, 1, @hann);
[ATrans, Af] = asd2(trans,1/sampf, 3, 1, @hann);

transErr=smooth(ATrans,10);
tiltErr=smooth(ATilt,10);
errFreq=Af;
end