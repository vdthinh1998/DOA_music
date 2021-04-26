clear all 
%close all 
clc 
f_sonar = 12000;
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
U_out = ihtrans(U);