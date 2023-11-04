%% *Homework2_2*
%% Programmers
% Mohammad Mahdi Elyasi - 9823007
%
% Moein Nasiri - 9823093
%% Clear Workspace
close all;
clear ;
clc;
%% Homework_1

[x_t, fs] = audioread('Audio01.wav');
x_t = x_t';
t_axis = linspace(0, length(x_t) / fs, length(x_t));
f0 = 10000;
% plot
plot(t_axis, x_t);

plot(fftshift(abs(fft(x_t))));

% sound(x_t, fs);
% pause(length(x_t) / fs);

s = 2 * cos((2 * pi * f0) * t_axis);
FD_FIR = filter_FIR;
y = FD_FIR.filter(x_t);
plot(t_axis, y);
y1 = y .* s;
plot(t_axis, y1);
y2 = FD_FIR.filter(y1);
imp_filter = load('filter.mat').Num;
reso_freq = 10000;
f_axis = linspace(-fs / 2, fs / 2, reso_freq);
% f_axis1 = linspace(-fs / 2, fs / 2, size(imp_filter));

plot(f_axis, fftshift(abs(fft(imp_filter, reso_freq))) / fs);
% sound(y2, fs);
% pause(length(y2) / fs);

figure('Name', 'Signals');
subplot(4, 2, 1)
plot(t_axis, x_t);
xlabel('Samples');
ylabel('Amplitude');
title('Input Signal');
grid on;

subplot(4, 2, 2)
plot(fftshift(abs(fft(x_t))) / fs);
xlabel('Samples');
ylabel('Amplitude');
grid on;
title('fft');

subplot(4, 2, 3)
plot(t_axis, y);
xlabel('Samples');
ylabel('Amplitude');
title('Filtered signal');
grid on;

subplot(4, 2, 4)
plot(fftshift(abs(fft(y))) / fs);
xlabel('Samples');
ylabel('Amplitude');
grid on;
title('fft');

subplot(4, 2, 5)
plot(t_axis, y1);
xlabel('Samples');
ylabel('Amplitude');
title('Product of carrier');
grid on;

subplot(4, 2, 6)
plot(fftshift(abs(fft(y1))) / fs);
xlabel('Samples');
ylabel('Amplitude');
grid on;
title('fft');

subplot(4, 2, 7)
plot(t_axis, y2);
xlabel('Samples');
ylabel('Amplitude');
title('Filtered signal');
grid on;

subplot(4, 2, 8)
plot(fftshift(abs(fft(y2))) / fs);
xlabel('Samples');
ylabel('Amplitude');
grid on;
title('fft');

figure('Name', 'Signals');
subplot(2, 1, 1);
plot(f_axis, 20 * log10(abs(fftshift(fft(imp_filter, reso_freq)))));
xlabel('Samples');
ylabel('Amplitude');
title('fft');
grid on;

subplot(2, 1, 2);
plot(f_axis, unwrap(angle(fftshift(fft(imp_filter, reso_freq)))));
xlabel('Samples');
ylabel('Amplitude');
% xlim([0 25000]);
grid on;
title('phase');
% 2.4
x2 = y2 .* s;
x4 = filter(imp_filter, 1, x2);
% sound(x4, fs);
% pause(length(x4) / fs);

%% Homework_2

f0 = 100;
fs = 500;
t = 0:1 / fs:2;
x_1 = sin(2 * pi * f0 * t);
x_2 = chirp(t, 400, 1, 200);
x_3 = 50 * [zeros(1, 249) 1 zeros(1, 751)];
x3 = x_1 + x_2 + x_3;
figure('Name', 'FFT');
plot(abs(fft(x)));
xlabel('freq(Hz)');
ylabel('amp');
title('fft');

w1 = hamming(256);
figure('Name', 'Hamming');
spectrogram(x3, w1, 64, 128, 1000, 'yaxis');
title('Spectogram 256');

w2 = hamming(512);
figure('Name', 'Hamming');
spectrogram(x3, w2, 64, 128, 1000, 'yaxis');
title('Spectogram 512');
