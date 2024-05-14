Ts = 1/10000;
t = linspace(0,0.01,10000); %range in time domain 

mt=5*cos(200*pi*t)+10*cos(400*pi*t); % defining message signal

ct=2*cos(2000*pi*t); % defining carrier signal

xt=mt+ct; %defining modulated signal

zt=60*xt+xt.^2; %defining output signal of non-linear device

%plot(t,zt); %plotting each function plot(t,mt);plot(t,xt)
%xlabel('time (s)')
%ylabel('Amplitude of z(t)')
%title('Plot of the non-linear device signal')
 
%//////////////////// B ////////////////
mf = fft(mt);
xf = fft(xt);
zf = fft(zt);

fs = 1/Ts;

N = length(t);

fshift = (-N/2:N/2-1)*(1/(t(2)-t(1)))/N;

%ff1shift = fftshift(mf);
%ff2shift = fftshift(xf);
ff3shift = fftshift(zf);

%stem(fshift,abs(ff1shift))
%stem(fshift,abs(ff2shift))
%stem(fshift,abs(ff3shift))

%xlabel('Frequency (Hz)')
%ylabel('Magnitude')
%title('Plot of the Z(f)')
%xlim([-1000,1000]);

%//////////////////// C ////////////////

center_frequency = 1000; % Center frequency in Hz
bandwidth = 401; % Bandwidth in Hz

% Find indices corresponding to the bandpass filter
lower_cutoff = center_frequency - bandwidth/2;
upper_cutoff = center_frequency + bandwidth/2;

% Apply bandpass filter
ff3shift_filtered = zeros(size(ff3shift));
ff3shift_filtered(abs(fshift) >= lower_cutoff & abs(fshift) <= upper_cutoff) = ff3shift(abs(fshift) >= lower_cutoff & abs(fshift) <= upper_cutoff);


% Inverse Fourier transform to obtain time-domain signal y(t)
yt = ifft(ifftshift(ff3shift_filtered));

% Plot the time-domain signal
plot(t, yt);
xlabel('Time (s)');
ylabel('Amplitude');
title('Filtered signal y(t)');


% Envelope detection
env = envelope(yt);
mtt = fft(env);


ff4shift = fftshift(mtt);
stem(fshift,abs(ff4shift))
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Plot of the demodulated signal m(f)')
xlim([-1000,1000]);

% Plot the envelope in the time domain
%plot(t,env);
%hold on
%xlabel('Time (s)');
%plot(t,yt);
%hold off
%ylabel('Amplitude');
%title('demodulated signal m(𝑡)');









