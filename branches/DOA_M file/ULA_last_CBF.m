clear all 
close all 
clc 

%THAM SO HE ANTEN UCA
%Buoc song cua tin hieu(m)
lamda=0.328;

%So phan tu cua mang anten
M=6;                 

%Khoang cach giua cap anten lien tiep so voi buoc song
d=0.5*lamda;

%GAIN
gain=0;

%THAM SO NGUON TIN HIEU DEN
%So nguon tin hieu 
D=2;

%GOC TOI CUA CAC NGUON TIN HIEU
angles=[90 92]*(pi/180);   

%SNR cua cac nguon tin hieu (dB)-ung voi moi kenh I va Q
SNRdB=15;
SNRdBs=SNRdB*[1 1];

%THAM SO CHUNG
%He so song
k=2*pi/lamda;         

%So mau tin hieu thu
Nb=1000; 


%Tao ma tran vecto dau vao tin hieu ban dau S[D,Nb] va ma tran vecto lai A_tmp(D,M)
for i=1:D 
    S(i,:)=(20^(SNRdBs(i)/10))*1*(randn(1,Nb)+j*randn(1,Nb));
    A_tmp(i,:)=10^(gain/10)*exp(j*k*(0:M-1)*d*(cos(angles(i)))); 
end

%Tao ma tran nhieu N[M,Nb] bien do moi kenh bang 1.
N=1*(randn(M,Nb)+j*randn(M,Nb));

%Tao ma tran du lieu thu duoc boi mang anten U[M,Nb]
A=A_tmp.';
U=A*S+N;

%AP DUNG THUAT TOAN MUSIC DE TIM DOA

%Tinh covarian cua tin hieu vao
Ruu=U*U'/Nb;

%Pho khong gian cua tin hieu
i=0; 
for theta=0:.1:180
    i=i+1; 
    A0_tmp=10^(gain/10)*exp(j*k*(0:M-1)*d*(cos(theta*pi/180))); 
    A0=A0_tmp.';
    %P(i)=10*log((A0'*A0)/(A0'*eigVectorNoise*eigVectorNoise'*A0))/10;
    %P(i)=10*log((A0'*A0)/(A0'*Ruu^(-1)*A0))/10;
    P(i)=10*log((A0'*Ruu*A0)/(A0'*A0))/10;
end

%Bieu dien
theta=0:.1:180;
plot(theta,real(P),'k','linewidth',2); 
xlabel('DOA(do)');
ylabel('Pho khong gian CAPON(dB)'); 
hold on;
%Ghi chu:
%Gia thiet:
%-Tin hieu gom I va Q ??u co phan bo Gaussian trang -> s=I+jQ va I va Q deu
%dung ham randn -normal distribution co var=1.
%-Nhieu la AWGN va tac dong len ca hai kenh I va Q nen n=nI+jnQ va cung
%dung ham randn.
%SNR ung voi moi kenh I va Q cua tin hieu.