clear all 
%close all 
clc 
f_sonar = 12000;
fs = 192000;
c_sound = 1500;
%THAM SO HE ANTEN ULA
%Buoc song cua tin hieu(m)
lamda=c_sound/f_sonar;

%So phan tu cua mang anten
M=4;                 

%Khoang cach giua cap anten lien tiep so voi buoc song
d=0.5*lamda;

%GAIN
gain=0;

%THAM SO NGUON TIN HIEU DEN
%So nguon tin hieu 
D=2;

%GOC TOI CUA CAC NGUON TIN HIEU
angles=[90 60]*(pi/180);   

%SNR cua cac nguon tin hieu (dB)-ung voi moi kenh I va Q
SNRdB=15;
SNRdBs=SNRdB*[1 1];

%THAM SO CHUNG
%He so song
k=2*pi/lamda;         

%So mau tin hieu thu
Nb=960; 


%Tao ma tran vecto dau vao tin hieu ban dau S[D,Nb] va ma tran vecto lai A_tmp(D,M)
for i=1:D 
    %S(i,:)=(20^(SNRdBs(i)/10))*1*(randn(1,Nb)+j*randn(1,Nb));
    %S(i,:)=(20^(SNRdBs(i)/10))*1*(randn(1,Nb));
    
    t = 0:1/fs:0.005-1/fs; % phat sin trong 0.0005s
    S(i,:) = 10*sin(2*pi*f_sonar*t);
    A_tmp(i,:)=10^(gain/10)*exp(j*k*(0:M-1)*d*(cos(angles(i)))); 
end
subplot(2,1,1);
plot(S(1,:));
hold on
plot(S(2,:));
%Tao ma tran nhieu N[M,Nb] bien do moi kenh bang 1.
N=1*(randn(M,Nb)+j*randn(M,Nb));

%Tao ma tran du lieu thu duoc boi mang anten U[M,Nb]
A=A_tmp.';
%U=A*S+N;
U=A*S+N;

%AP DUNG THUAT TOAN MUSIC DE TIM DOA

%Tinh covarian cua tin hieu vao
Ruu=U*U'/Nb;

%Xac dinh gia tri rieng va vector rieng cua covarian cua tin hieu loi vao
[eigVector,eigValue]=eig(Ruu);

eigValueMax=max(max(eigValue));

%disp(eigValue);

%Xac dinh so nguon tin hieu den
%signals=length(find(diag(eigValue)>eigValueMax/1000000));
%disp(signals);
%signals=1;
signals=2;
%signals=4;
%Xac dinh cac vector rieng cua khong gian nhieu
eigVectorNoise=eigVector(:,1:M-signals);

%Pho khong gian cua tin hieu
i=0; 
for theta=0:.1:180
    i=i+1; 
    A0_tmp=10^(gain/10)*exp(j*k*(0:M-1)*d*(cos(theta*pi/180))); 
    A0=A0_tmp.';
    P(i)=10*log((A0'*A0)/(A0'*eigVectorNoise*eigVectorNoise'*A0))/10;
    %P(i)=((A0'*A0)/(A0'*eigVectorNoise*eigVectorNoise'*A0))
end

%Bieu dien
theta=0:.1:180;
subplot(2,1,2);
plot(theta,real(P),'k','linewidth',2); 
xlabel('DOA(do)');
ylabel('Pho khong gian MUSIC(dB)'); 
hold on;