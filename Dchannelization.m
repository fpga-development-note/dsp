clear all;close all;clc;

fs=2.5e+9;%采样率2.4GHz
f=1.0e+8;%采样频率100Mhz
IF=80e+7;%中频350MHz
MF=3.6e+9;%高频3.6GHz

N=1024;
Signal_DDC_N=0:1:N-1;
Signal_DDC_T=Signal_DDC_N./fs;
Signal_DDC=sin(2*pi*f*Signal_DDC_T);%取f=100Hz的1000个点作为离散数据
 
 
N=length(Signal_DDC);
nfft=2^nextpow2(N);%求取FFT点数，fft点数1024，nextpow2是求取2的指数位宽，如nextpow2(1001)=10
t_axis=(0:N-1)./fs;
f_axis=(0:nfft-1)./nfft*fs-fs/2;

figure;
subplot(211)
plot(t_axis,real(Signal_DDC))
%plot(t_axis.*1e6,real(Signal_DDC));
xlabel('时间')
title('输入信号')
subplot(212)
Signal_DDC_fft=(fftshift(fft(Signal_DDC,nfft)));%fft运算，采样点数为1000，输出点数为1024，
plot(f_axis./1e6,(abs(Signal_DDC_fft)))
xlabel('频率')
title('输入信号的频谱')

%%原型滤波器设计
fI=IF;
K=16;%划分16个信道
Channel_Freq_Range=[((0:K-1)-(K-1)/2).*fs/K-fs/K/2;((0:K-1)-(K-1)/2).*fs/K+fs/K/2]./1e6;

n0 = 127;               % Filter order
f0 = [0 0.03125 0.0625 1];    % Frequency band edges
A0 = [1  1  0 0];      % Amplitudes
h_LP = firpm(n0,f0,A0);

%h_LP=fir1(1023,1/K,'low');%fir1(n,Wn)：n代表设计滤波器阶数，Wn代表归一化截止频率，截止频率/奈奎斯特频率。计算方式为2*pi*f/Fs 如截止频率f=300Hz，采样率Fs=1000Hz，则Wn=0.6
%fvtool(h_LP,1);
freqz(h_LP,1,512);

%%滤波系数通道提取
M=length(h_LP);%滤波器阶数
Q=fix(M/K);%每个通道滤波器的阶数，128/16=8
H=zeros(K,Q);
for d=1:K
    H(d,:)=h_LP(d:K:(Q-1)*K+d);%提取每个通道滤波器系数
end

%%补0
tic;%调用计时器，显示单位为秒
temp=mod(length(Signal_DDC),K);%取余
if temp~=0
    Signal_DDC0=[Signal_DDC(1:end),zeros(1,K-temp)];%添0
    X=reshape(Signal_DDC0,K,length(Signal_DDC0)/K);%数据16倍抽取，将1024点变成16X64矩阵
else
    X=reshape(Signal_DDC(1:end),K,length(Signal_DDC)/K);
end
%X=flipud(X);%将矩阵上下翻转，为了进行滤波器卷积运算

%与(-1)^m相乘
[rx,L]=size(X);
if mod(K,2)==0
    X=X.*repmat((-1).^(0:L-1),K,1);%repmat为重复矩阵，将(-1)^(0:L-1)重复K行1列，用于每个通道滤波之前的选择
end
%多相滤波  
Y=zeros(K,L);
for d=1:K
    Y(d,:)=Filter_FFT(X(d,:),H(d,:));%数字滤波
end
for ll=1:L
    temp=Y(:,ll).*(-1).^(0:K-1).';
    temp=temp.*exp(j*(0:K-1).'*pi/K);
    Y(:,ll)=ifft(temp,K).*K;
end
toc;%计时器结束，信道化结束

%信道化结果观察
figure;
range=1:16;%观察1-16通道的数据
T_range=length(range);
for d=1:length(range)
    if mod(length(range),4)==0
        subplot(length(range)/2,4,d);
    else
        subplot((length(range)+1)/2,2,d);
    end
    t_axis=((0:L-1))./fs.*K+(range(d)-1)./fs;
    plot(t_axis.*1e6,abs(Y(range(d),:)))
    ylim([0 max(max(abs(Y)))+1])
end
for m=1:460
    f       = 150+m*5;
    sig     = floor(127*sin(2*pi*f*t));

    sig_vec = reshape(sig,32,N/32)';
    y       = fft(sig_vec,[],2);
    y       = y(:,1:16);
 
 
    [max_value,max_idx] = max(abs(y(1,:)));
    phase_ori = angle(y(:,max_idx))';
 

    

    sub_pse_ori = phase_ori(1,2:end) - phase_ori(1,1:end-1);
    for n = 1:length(sub_pse_ori)
       if(sub_pse_ori(1,n)>=pi) 
           sub_pse_hand(1,n) = sub_pse_ori(1,n)-2*pi;
       elseif(sub_pse_ori(1,n)<=-pi) 
           sub_pse_hand(1,n) = sub_pse_ori(1,n)+2*pi;
       else
           sub_pse_hand(1,n) = sub_pse_ori(1,n);
       end
    end
    
    pse_unwrap     = unwrap(phase_ori);
    sub_pse_unwrap = pse_unwrap(1,2:end) - pse_unwrap(1,1:end-1);
    diff = sub_pse_hand-sub_pse_unwrap;
    
  
    sub_pse_0 = sub_pse_hand;
    ratio_0 = sum(sub_pse_0(1,1:8))/8/2/pi;
    sub_pse_1 = sub_pse_unwrap;
    ratio_1 = sum(sub_pse_1(1,1:8))/8/2/pi;
    a_result(1,m) = max_idx-1+ratio_0;
    b_result(1,m) = max_idx-1+ratio_1;
    c_result(1,m) = f/FS*32;
end

