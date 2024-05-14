Fs= 10000;
t = -2:1/Fs:6;
N = length(t);
f = (-N/2:N/2-1)*(1/(t(2)-t(1)))/N;

%defining message signal
mt = zeros(1,length(t));
mt(0<=t & t<=2) = 1;
mt(2<t & t<=4) = -1;

%fourier transform of m(t)
m = fft(mt);
mf = fftshift(m)./Fs;

%ploting message signal
figure;
subplot(2,1,1);
plot(t,mt);

xlabel('time (s)')
ylabel('|m(t)|')
title('Plot of the message signal')

subplot(2,1,2);
plot(f,abs(mf));

xlabel('f(Hz)')
ylabel('|M(f)|')
title('Magnitude Spectrum of the message signal')

%defining carrier signal
ct=5*cos(2*pi*250*t);

% defining ∅(𝑡)
kf=50;
integral_mt = cumtrapz(t, mt);
fi=2.*pi.*kf.*integral_mt;

% plotting ∅(𝑡)
figure;
subplot(2,1,1);
plot(t,fi);

xlabel('time (s)')
ylabel('|∅(t)|')
title('Plot of the ∅(t)')

% ∅(f)
fi_fourier=fft(fi);
fi_fourier_shifted=fftshift(fi_fourier);

%plotting ∅(f)
subplot(2,1,2);
plot(f,abs(fi_fourier_shifted)./Fs);

xlabel('f(Hz)')
ylabel('|∅(f)|')
title('Magnitude Spectrum of the∅(f)')


% Defining y(t)
yt=5.*cos((2.*pi.*250.*t)+fi);

% Fourier Transform of y(f)
y_fourier=fft(yt);
y_fourier_shifted=fftshift(y_fourier)./Fs;

figure;
subplot(2,1,1);
plot(t,yt);

xlabel('time (s)')
ylabel('|y(t)|')
title('Plot of the y(t)')

subplot(2,1,2);
plot(f,abs(y_fourier_shifted));

xlabel('f(Hz)')
ylabel('|y(f)|')
title('Magnitude Spectrum of the y(f)')

%demodulasyon
y_diff= diff(yt)./(t(2)-t(1));

% differantiator frequency domain
y_diff_f=fft(y_diff);
y_diff_ff=fftshift(y_diff_f)./Fs;

%diyottan geçirme
y_diyot = max(0, y_diff);

%y_diyot frekans domaininde
y_diyot_f= fft(y_diyot);
y_diyot_ff= fftshift(y_diyot_f)./Fs;

%graph
figure;
subplot(2,2,1);
plot(t(1:end-1),y_diff);

xlabel('time (s)')
ylabel('|ydiff(t)|')
title('Plot of the signal after Differentiator')


subplot(2,2,2);
plot(t(1:end-1),y_diyot);
xlabel('time (s)')
ylabel('|ydiyot(t)|')
title('Plot of the signal after Diyot')

% Assuming length(y_diyot_ff) is the correct length for the frequency domain
f_diff = (-length(y_diff_ff)/2:length(y_diff_ff)/2-1) * (1/(t(2)-t(1)))/length(y_diff_ff);


% Assuming length(y_diyot_ff) is the correct length for the frequency domain
f_diyot = (-length(y_diyot_ff)/2:length(y_diyot_ff)/2-1) * (1/(t(2)-t(1)))/length(y_diyot_ff);

% Plotting frequency domain of y_diyot_ff
subplot(2,2,4);
plot(f_diyot, abs(y_diyot_ff));

xlabel('f(Hz)')
ylabel('|ydiyot(f)|')
title('Magnitude Plot of the signal after Diyot ')

% Plotting frequency domain of y_diff_ff
subplot(2,2,3);
plot(f_diff, abs(y_diff_ff));

xlabel('f(Hz)')
ylabel('|ydiff(f)|')
title('Magnitude Plot of the signal after Differantiator ')


%RC Filter
R=200000;
C=1e-6;
tao=R.*C;
drop_off = exp(-((t(2)-t(1))./tao));

y_envelope=ones(size(y_diyot));

for i= 2: length(y_diyot)
    if y_envelope(i-1) < y_diyot(i)

        y_envelope(i) = y_diyot(i);

    else 
        y_envelope(i)=drop_off.*y_envelope(i-1);
    end
end

%y(t)türevinde eğer envelopdaki değer y(tin türevinden(y_diyot) küçükse envelopa
%eşit değilse ve büyükse exponansiyel formul kadar azalcak (mantığım)

% plot y_envelope
figure;
plot(t(1:end-1), y_envelope);
xlabel('time (s)')
ylabel('|yenvelope(t)|')
title('Plot of the Envelope of y(t)')

%y_envelope fourier

y_envelope_f=fft(y_envelope);
y_envelope_fft= fftshift(y_envelope_f)./Fs;

figure
plot(f_diyot, abs(y_envelope_fft));
xlabel('f(Hz)')
ylabel('|yenvelope(f)|')
title('Magnitude Spectrum Plot of the y envelope')


%Low pass filter

center_frequency = 0; % Center frequency in Hz
bandwidth = 101; % Bandwidth in Hz

% Find indices corresponding to the bandpass filter
lower_cutoff = center_frequency - bandwidth/2;
upper_cutoff = center_frequency + bandwidth/2;

% Apply bandpass filter
y_envelope_filtered = zeros(size(y_envelope_fft));
y_envelope_filtered(abs(f) >= lower_cutoff & abs(f) <= upper_cutoff) = y_envelope_fft(abs(f) >= lower_cutoff & abs(f) <= upper_cutoff);

% inverse fourier
y_envelope_filtered_time=ifft(ifftshift(y_envelope_filtered));

%plotting
figure
subplot(2,1,1)
plot(t(1:end-1),y_envelope_filtered_time);
xlabel('time (s)')
ylabel('|y envelope- filtered(t)|')
title('Plot of the Filtered Envelope of y(t)')

subplot(2,1,2)
plot(f_diyot,y_envelope_filtered);
xlabel('f(Hz)')
ylabel('|y envelope filtered(f)|')
title('Magnitude Spectrum Plot of the y envelope filtered')

% DC Blocking
y_envelope_filtered_dcblocked = y_envelope_filtered;
dc_index = find(f_diyot == 0); % Find the index corresponding to DC component

% Set the DC component to zero
y_envelope_filtered_dcblocked(dc_index) = 0;

% inverse fourier after DC blocking
y_envelope_filtered_dcblocked_time = ifft(ifftshift(y_envelope_filtered_dcblocked));

% Plotting
figure
subplot(2,1,1)
plot(t(1:end-1), y_envelope_filtered_dcblocked_time);
xlabel('time (s)')
ylabel('|y-DC blocked(t)|')
title('Plot of the Filtered and DC Blocked Envelope of y(t)')

subplot(2,1,2)
plot(f_diyot, abs(y_envelope_filtered_dcblocked));
xlabel('f(Hz)')
ylabel('|y-DC blocked(f)|')
title('Magnitude Spectrum Plot of the Filtered and DC Blocked y envelope')

