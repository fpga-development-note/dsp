clear all;close all;clc;

% 参数设置
Fs = 2500e6;        % 采样率 2500 MHz
f1=10e6;
f2=100e6;
f3=400e6;
N = 10240;           % 采样点数
t = (0:N-1)/Fs;     % 时间向量


% Signal_DDC_N=0:1:N-1;
% Signal_DDC_T=Signal_DDC_N/fs;
Signal_DDC1=cos(2*pi*f1*t);
Signal_DDC2=cos(2*pi*f2*t);
Signal_DDC3=cos(2*pi*f3*t);
% 添加一些噪声使信号更真实
noise_power = 0.01;
noise_signal = sqrt(noise_power)*randn(size(Signal_DDC1));

%Signal = Signal_DDC1+Signal_DDC2+Signal_DDC3+noise_signal;
Signal = Signal_DDC1+Signal_DDC2+Signal_DDC3;

% 计算频谱
X = fft(Signal, N);      % FFT变换
X_abs = abs(X);     % 取幅度
X_abs = fftshift(X_abs); % 将零频移到中心


% 频率向量 (MHz)
f = (-N/2:N/2-1) * (Fs/N) / 1e6; % 转换为MHz单位


%滤波器设计
n0 = 127;               % Filter order
f0 = [0 0.03125 0.0625 1];    % Frequency band edges
A0 = [1  1  0 0];      % Amplitudes
bPM0 = firpm(n0,f0,A0);

% 计算滤波器频率响应
[H_proto, f_proto] = freqz(bPM0, 1, 1024, Fs);
%y=filter(bPM0,1,Signal);

%滤波计算
y=Filter_FFT(Signal,bPM0);%数字滤波

After_Signal_DDC_fft=fftshift(abs(fft(y,N)));

figure;
% 子图1: 原型滤波器频率响应
plot(f_proto/1e6, 20*log10(abs(H_proto)), 'LineWidth', 2);
xlabel('频率 (MHz)');
ylabel('幅度 (dB)');
title('原型滤波器频率响应');
grid on;
xlim([0, 200]);
hold on;


figure;
subplot(4,1,1);
plot(t,real(Signal));
%plot(t_axis.*1e6,real(Signal_DDC));
xlabel('时间');
title('滤波前输入信号');

subplot(4,1,2);
plot(t,real(y));
%plot(t_axis.*1e6,real(Signal_DDC));
xlabel('时间');
title('滤波后输入信号');


subplot(4,1,3);
plot(f,X_abs);
xlabel('频率');
title('输入信号的频谱');

subplot(4,1,4);
plot(f,After_Signal_DDC_fft);
xlabel('频率');
title('输入信号的频谱');




