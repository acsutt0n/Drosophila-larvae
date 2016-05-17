% important notes: 1) make sure your moving window is at least 2x as long as
% the slowest your interested in (i.e. for a 10 hz oscillation, your window
% should be at least 200 ms) 2) make sure your sampling frequency is at
% least 2x as fast as the fastest event your interested in (i.e. for a 100
% hz osciilation, make sure your Fs is at least 200 Hz).

% name your file data.txt, make sure is 1 column of numbers
function moving_PSD(trace)

%load AK.txt
%x = AK(:,1);

% YOUR INPUT: moving time spectrogram, units = ms, x = data, 100 = moving window width,
% 50 = step width, [0:100] = output bandwidth, 1000 = sampling frequency
% (Hz)
% FOR INFO ON OUTPUT: type "help spectrogram" into Matlab command window

[S,F,T,P] = spectrogram(trace,100000,50000,[0.001:60],50000);

% graph specifics
surf(T,F,log(P))
shading interp
view(0,90)