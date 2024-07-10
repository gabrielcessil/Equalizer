close all;
clear sound;
load handel.mat

% Free copy right song found in:
% https://pixabay.com/music/search/music%20with%20lyrics%20contest/
% ** INPUTS **

% Choose version played: (-1) Noised; (0) Original; (1) Filtered
version = 0;  

% Update the path for your specific case
filename = 'C:\Users\gabri\Desktop\Posts\Audio processing - Wordpress/rum-and-finesse-220404.mp3';
SNR = 20; % Added signal-to-noise ratio 

% Design relevant info 
Bass = [50 500]; % (Hz) Lower frequencies
LowMedium = [500 2000]; % (Hz) Low-Medium frequencies

% Gains for frequency window
GainBass = 1; % SubBass gain (0~1)
GainMedium = 1; % Bass gain (0~1)

% Filters specification 
Amax = 8; % Maximum ripple in Passband(dB)
Amin = 10; % Minimum attenuaion in Stopband (dB)
tran = 0.1; % Transition band fraction (0~1): How much of each band is transitional 


 
% ** COMPUTING **

% Lambda functions for conversions
rads_to_Hz = @(w) w/(2 * pi);

% Reading the audio file
[audio_signal, Fs] = audioread(filename); 

% Converting stereo audio to single channel
audio_signal = (audio_signal(:,1) + audio_signal(:,2))/ max([audio_signal(:,1).'  audio_signal(:,2).'] ); 

% Add noise for didatic porpouse
Signal_To_Noise_Ratio = 20*log10(SNR);
noisy_audio_signal = awgn(audio_signal,Signal_To_Noise_Ratio);

N = length(noisy_audio_signal); % Audio noisy_audio_signal lenght
Ts = 1/Fs; % Audio noisy_audio_signal sampling time (s)
Ws = 2*pi*Fs; % Audio noisy_audio_signal sampling frequency (rad/s)

% Designing the Analog Passband Filter 
fprintf('\nBass Frequencies:\n');
[T1,Ts_filtros1] = PassBand_Butter(Bass,Amax,Amin,tran,GainBass);
fprintf('\nMedium Frequencies:\n');
[T2,Ts_filtros2]= PassBand_Butter(LowMedium,Amax,Amin,tran,GainMedium);

% Estimating the new sampling time for the filtering as multiple of Ts
Ts_filtros = Ts*floor(min([Ts_filtros1 Ts_filtros2])/Ts);
Fs_filtros = 1/Ts_filtros;
fprintf('\nOriginal signal sampling time: %.3f us\n', Ts*10^6);
fprintf('\nFiltering sampling time: %.3f us\n', Ts_filtros*10^6);

% Discretinzing the analog filters
T1z = c2d(T1, Ts_filtros, 'tustin');
T2z = c2d(T2, Ts_filtros, 'tustin');

% Resampling the original signal to match the filters
audio_signal_resampled = ResambleSignal(noisy_audio_signal,Fs,Fs_filtros);
% Filter the signal
filter_response_1 = filter(T1z.Numerator{1},T1z.Denominator{1},audio_signal_resampled);
filter_response_2 = filter(T2z.Numerator{1},T2z.Denominator{1},audio_signal_resampled);

% Combining filtered channels
signal_out = filter_response_1 + filter_response_2;
% Returning to original audio level
signal_out = (signal_out / max(signal_out)) * max(noisy_audio_signal);

% Resampling the output to match the original sampling
signal_out_resampled = ResambleSignal(signal_out,Fs_filtros,Fs);

% Play audios
if version>0
    sound(signal_out_resampled, Fs);
elseif version<0
    sound(noisy_audio_signal, Fs);
else
    sound(audio_signal, Fs);
end

% Saving audios
audiowrite('noisy_audio.ogg',noisy_audio_signal,Fs)
audiowrite('filtered_audio.ogg',signal_out_resampled,Fs)

% ** PLOTTING **
figure(1);
PlotMagnitudes(audio_signal,signal_out_resampled,noisy_audio_signal, Ts, SNR)
figure(2);
PlotBode(T1z, T2z, Bass, LowMedium)
figure(3);
PlotFFT(audio_signal,signal_out_resampled,noisy_audio_signal, Bass, LowMedium, Ts, Ts_filtros)


% Fourier transform of original noisy_audio_signal

function PlotFFT(audio_signal,signal_out_resampled,noisy_audio_signal,Band1, Band2, Ts, Ts_filtros)
    add_vertical_line = @(x, color) line([x x], ylim, 'Color', color, 'LineStyle', '--', 'LineWidth', 1.5, 'HandleVisibility','off');
    
    signal_FFT = fft(audio_signal); % Shift in frequency
    signal_FFT = fftshift(signal_FFT); % Shift in frequency
    mag_signal_FFT = mag2db(abs(signal_FFT)); % Convert to dB
    
    noisy_signal_FFT = fft(noisy_audio_signal); % Output audio Fourier transform
    noisy_signal_FFT = fftshift(noisy_signal_FFT); % Output audio Fourier transform
    noisy_mag_signal_FFT = mag2db(abs(noisy_signal_FFT)); % Convert to dB
    
    signal_out_FFT = fft(signal_out_resampled); % Output audio Fourier transform
    signal_out_FFT = fftshift(signal_out_FFT); % Shift in frequency
    mag_signal_out_FFT = mag2db(abs(signal_out_FFT)); % Convert to dB
    
    Fs =  1/Ts;
    Fs_filtros = 1/Ts_filtros;
    N = length(audio_signal); % Audio lenght
    N_output = length(signal_out_resampled); % Output audio lenght
    
    fshift_output = (-N_output/2:N_output/2-1)*(Fs/N_output); % Shifted frequency
    fshift = (-N/2:N/2-1)*(Fs/N); % Shifted frequency
    ylimFFT_max = real(max([mag_signal_FFT.' noisy_mag_signal_FFT.' mag_signal_out_FFT.']));
    ylimFFT_min = real(min([mag_signal_FFT.' noisy_mag_signal_FFT.' mag_signal_out_FFT.']));
    
    subplot(3,1,1);
    hold on
    plot(fshift, mag2db(abs(signal_FFT)));
    add_vertical_line(Band1(1), 'r');
    add_vertical_line(Band1(2), 'r');
    add_vertical_line(Band2(1), 'm');
    add_vertical_line(Band2(2), 'm');
    add_vertical_line(Fs_filtros, 'k');
    title('Original Audio Spectrum');
    xlabel('Hz');
    ylabel('Amplitude');
    ylim([ylimFFT_min*1.1 ylimFFT_max*1.1]);
    xlim([-24000 24000]);
    hold off
    subplot(3,1,2);
    hold on
    plot(fshift, mag2db(abs(noisy_signal_FFT)));
    add_vertical_line(Band1(1), 'r');
    add_vertical_line(Band1(2), 'r');
    add_vertical_line(Band2(1), 'm');
    add_vertical_line(Band2(2), 'm');
    add_vertical_line(Fs_filtros, 'k');
    title('Noisy Audio Spectrum');
    xlabel('Hz');
    ylabel('Amplitude');
    ylim([ylimFFT_min*1.1 ylimFFT_max*1.1]);
    xlim([-24000 24000]);
    hold off
    subplot(3,1,3);
    hold on
    plot(fshift_output, mag2db(abs(signal_out_FFT)))
    add_vertical_line(Band1(1), 'r');
    add_vertical_line(Band1(2), 'r');
    add_vertical_line(Band2(1), 'm');
    add_vertical_line(Band2(2), 'm');
    add_vertical_line(Fs_filtros, 'k');
    title('Filtered Audio Spectrum');
    xlabel('Hz');
    ylabel('Amplitude');
    ylim([ylimFFT_min*1.1 ylimFFT_max*1.1]);
    xlim([-24000 24000]);
    hold off
end
function PlotBode(T1z, T2z, Band1, Band2)
    
    % Lambda function to add vertical lines
    add_vertical_line = @(x, color) line([x x], ylim, 'Color', color, 'LineStyle', '--', 'LineWidth', 1.5, 'HandleVisibility','off');
    [mag1, phase1, wout1] = bode(T1z); % Get the Bode plot data
    [mag2, phase2, wout2] = bode(T2z); % Get the Bode plot data
    freq1 = wout1 / (2 * pi); % Convert angular frequency (rad/s) to Hz
    freq2 = wout2 / (2 * pi); % Convert angular frequency (rad/s) to Hz
    mag1_dB = mag2db(squeeze(mag1)); % Convert magnitude to dB
    mag2_dB = mag2db(squeeze(mag2)); % Convert magnitude to dB
    phase1_deg = squeeze(phase1); % Convert phase to degrees
    phase2_deg = squeeze(phase2); % Convert phase to degrees

    subplot(2,1,1);
    hold on
    semilogx(freq1, mag1_dB, 'b', 'LineWidth', 1.5);
    semilogx(freq2, mag2_dB, 'r', 'LineWidth', 1.5);
    add_vertical_line(Band1(1), 'b');
    add_vertical_line(Band2(1), 'b');
    add_vertical_line(Band1(2), 'r');
    add_vertical_line(Band2(2), 'r');
    title('Bode Diagram');
    ylabel('Magnitude (dB)');
    legend('Bass Channel Filter', 'LowMedium Channel Filter');
    grid on;
    hold off


    subplot(2,1,2);
    hold on
    semilogx(freq1, phase1_deg, 'b', 'LineWidth', 1.5);
    semilogx(freq2, phase2_deg, 'r', 'LineWidth', 1.5);
    xlabel('Frequency (Hz)');
    ylabel('Phase (degrees)');
    grid on;
    add_vertical_line(Band1(1), 'b');
    add_vertical_line(Band1(2), 'b');
    add_vertical_line(Band2(1), 'r');
    add_vertical_line(Band2(2), 'r');
    hold off
end
function PlotMagnitudes(audio_signal,signal_out_resampled,noisy_audio_signal,Ts, SNR)
    
    N = length(audio_signal); % Audio lenght
    N_output = length(signal_out_resampled); % Output audio lenght
    
    
    % Plot the audio signals magnitudes
    time_signal = (0:(N-1))*Ts; % Time noisy_audio_signal (s)
    subplot(4,1,1);
    plot(time_signal, audio_signal)
    title('Original Audio Signal')
    ylabel('Amplitude')
    xlim([0 max(time_signal)]);
    subplot(4,1,2);
    plot(time_signal, audio_signal-noisy_audio_signal, 'r')
    title("Added Noise (SNR="+SNR+")")
    ylabel('Amplitude')
    ylim([-1 1]);
    xlim([0 max(time_signal)]);
    subplot(4,1,3);
    plot(time_signal, noisy_audio_signal)
    title('Noisy Audio Signal')
    ylabel('Amplitude')
    xlim([0 max(time_signal)]);
    subplot(4,1,4);
    time_signal_Filters = (0:(N_output-1))*Ts; % Time noisy_audio_signal (s)
    plot(time_signal_Filters, signal_out_resampled)
    title('Filtered Audio Signal')
    xlabel('Time (s)')  
    ylabel('Amplitude')
    xlim([0 max(time_signal)]);

end
function resambled = ResambleSignal(signal, originalFs,desiredFs)
    [p,q] = rat(desiredFs / originalFs);
    resambled = resample(signal,p,q);
end 
function [T,Ts_filtros] = PassBand_Butter(Band,Amax,Amin,tran, Gain)
    [ws1,wp1,wp2,ws2] = GetBandFrequencySpec(Band(1),Band(2),tran);
    [ws, wp, w0] =  PassBand_to_NormalizedPassLow(ws1, wp1, wp2, ws2);
    E = sqrt(10^(Amax/10)-1);
    n = ceil(log10( (10^(Amin/10)-1)/(E.^2) )/(log10(ws)*2)); 
    [num,den] = butter(n,[ws1,ws2],'s');
    T = tf(num,den)*Gain;
    Ts_filtros = (2*pi)/(ws2*3);
    fprintf('Passband: %.1f < f < %.1f Hz\n' ,(wp1)/(2*pi), (wp2)/(2*pi));
    fprintf('Stopbands: %.1f < f, f > %.1f Hz\n', (ws1)/(2*pi), (ws2)/(2*pi));
end 

function [ws1,wp1,wp2,ws2] = GetBandFrequencySpec(f1,f2,tran)
    Hz_to_rads = @(f) 2 * pi * f;
    fs1 = min(f1,f2); % Beginnig of transition band (Hz)
    fs2 = max(f1,f2); % End of transition band (Hz)
    fp1 = fs1+ (fs2-fs1)*tran/2; % Lowest frequency on band (Hz)
    fp2 = fs2- (fs2-fs1)*tran/2; % Highest frequency on band (Hz)      

    
    wp1 = Hz_to_rads(fp1); % Lowest frequency on band (rad/s)
    wp2 = Hz_to_rads(fp2); % Highest frequency on band (rad/s)  
    ws1 = Hz_to_rads(fs1); % Beginnig of transition band (rad/s)
    ws2 = Hz_to_rads(fs2); % End of transition band (rad/s)
end

function [ws_, wp_, w0] = PassBand_to_NormalizedPassLow(ws1, wp1, wp2, ws2)
    w0 = sqrt(wp1*wp2);
    ws_ = (ws2-ws1)/(wp2-wp1);
    wp_ = 1;
end

function T = Butterworth(Amax, Amin, ws, wp, s_)

    % Butterworth filter design
    E = sqrt(10^(Amax/10)-1);
    %n = ceil(  log10(10^(Amin/10)-1)/ (2/log(ws))  );
    n = ceil(log10( (10^(Amin/10)-1)/(E.^2) )/(log10(ws)*2));
    
    poles = zeros(1, n);
    
    for i = 1:n
        sigma = E^(-1/n) * sin( ((2*i - 1) * pi) / (2*n) );
        omega = E^(-1/n) * cos( ((2*i - 1) * pi) / (2*n) );
        
        poles(i) = sigma + omega*1j; % Use 1j instead of j for imaginary unit
    end
    
    A = 1;
    for i = 1:n
        if mod(n,2) ~= 0 && i == ceil(n/2)
            A = A * (s_ + real(poles(i)));
        elseif i > (n/2)
            break;
        else
            A = A * (s_^2 + s_*real(poles(i) + poles(n-(i-1))) + real(poles(i)*poles(n-(i-1))));
        end
    end
    
    G0 = prod(poles);
    T = G0 / A;
    
    fprintf('Filter order: %.1f\n', n);

end

