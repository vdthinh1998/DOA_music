%Author: Huy Vuong 
%School: HaNoi University Of Science and Technology
%Project: MUSIC DOA
% aaaaaa
%=============================================================
% Hang so 
%=============================================================
clear all
clc
f = 12000;
v = 1500;
lamda=v/f;
M = 2; %so phan tu mang thu
gain = 0;
k=2*pi/lamda;
D=1; %so nguon tin hieu
d = lamda/2;
%=============================================================
% Filter
%=============================================================
%[U,fs] = audioread('D:/THUAT_TOAN_MUSIC/file_wav/28-11/test_12k_2811_3.wav');
[U,fs] = audioread('data_2kenh_radar_12k.wav');
[b,a]=butter(5,[11000,13000]/(fs/2),'bandpass');
filtsig=filter(b,a,U);  %filtered signal
maxkenh1 = max(filtsig(:,1));
maxkenh2 = max(filtsig(:,2));
heso = maxkenh2/maxkenh1; %Tim do lech bien do giua 2 kenh
y1 = heso*filtsig(:,1);
y2 = filtsig(:,2);
y = [ y1 y2 ];%Ghep 2 kenh 
yy = hilbert(y);% Chuyen tu real => phuc
Nb = length(y); %so mau
UU = yy';
dt = 1/fs;
t = 0:dt:(length(y)*dt)-dt;
subplot(3,1,1);
plot(t,y1);
hold on
plot(t,y2,'g');
xlabel('Seconds'); 
ylabel('Amplitude');
subplot(3,1,2);
plot(psd(spectrum.periodogram,y,'Fs',fs,'NFFT',length(y)));
%AP DUNG THUAT TOAN MUSIC DE TIM DOA
%=============================================================
%Tinh covarian cua tin hieu vao
%=============================================================

Ruu=UU*UU'/Nb;

%=============================================================
%Xac dinh gia tri rieng va vector rieng cua covarian cua tin hieu loi vao
%=============================================================
[eigVector,eigValue]=eig(Ruu);

eigValueMax=max(max(eigValue));

%disp(eigValue);

%Xac dinh so nguon tin hieu den
%signals=length(find(diag(eigValue)>eigValueMax/1000000));
%disp(signals);
signals = 1;
%signals = 2;
%signals = 4;

%Xac dinh cac vector rieng cua khong gian nhieu
eigVectorNoise=eigVector(:,1:M-signals);

%Pho khong gian cua tin hieu
i=1; 
for theta=1:.1:180
    i=i+1; 
    A0_tmp=10^(gain/10)*exp(j*k*(0:M-1)*d*(cos(theta*pi/180))); 
    A0=A0_tmp.';
    P(i)=10*log((A0'*A0)/(A0'*eigVectorNoise*eigVectorNoise'*A0))/10;
    %P(i)=((A0'*A0)/(A0'*eigVectorNoise*eigVectorNoise'*A0))
    %[maxP,index] = max(real(P));
end

for index=1:length(P)
    if (index == 1) && (P(1) > P(2))
        disp(index/10);
    elseif (index == length(P)) && (P(length(P)) > P(length(P)-1))
        disp(index/10);
    elseif (index~= length(P))&& (index ~= 1) && (P(index) > P(index-1)) && (P(index) > P(index+1))
        disp(index/10);
    end
end 
    
    
%result_DOA = index/10;
%display(result_DOA);

%=============================================================
%Bieu dien
%=============================================================
subplot(3,1,3);
theta=1:.1:180;
plot(theta,real(P(1,1:1791)),'k','linewidth',2); 
xlabel('DOA(do)');
ylabel('Pho khong gian MUSIC(dB)'); 
hold on;