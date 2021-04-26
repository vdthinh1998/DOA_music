tic
clear all 
close all 
clc 
M=5;                 %So phan tu cua mang anten 
Nb=1000;               %So mau tin hieu thu duoc
lamda=0.328;            %Buoc song cua tin hieu (m)
d=0.5*lamda;               %Khoang cach giua cac phan tu anten trong mang ULA (m)
D=2;                   %So nguon tin hieu 
angles=[90 92]*(pi/180); %Mang cac goc toi cua cac nguon tin hieu tuong ung
k=2*pi/lamda;     
gain=0;
SNRdB=15;
SNRdBs=SNRdB*[1 1];

%-------------------------------------------------
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
%-------------------------------------------------


%Tim ma tran cac vecto rieng cua U*UH

%Xac dinh gia tri rieng va vector rieng cua covarian cua tin hieu loi vao
[eigVector,eigValue]=eig(Ruu);
disp(eigValue);
eigValueMax=max(max(eigValue));
disp(eigValueMax);

%Xac dinh so nguon tin hieu den
%signals=length(find(diag(eigValue)>eigValueMax/1000));
%disp(signals);
signals=2;

%Xac dinh cac vector rieng cua khong gian nhieu
disp(eigVector);
eigVectorNoise=eigVector(:,1:M-signals);

%Xac dinh cac vector rieng cua khong gian tin hieu
eigVectorSignal=eigVector(:,M-signals+1:M);
%disp (eigVectorSignal);

%Xac dinh Q0 va Q1
Q0=eigVectorSignal(1:M-1,:);
Q1=eigVectorSignal(2:M,:);
disp(Q0);
disp(Q1);

%Xac dinh V
[U,D,V]=svd([Q0 Q1]);
disp (V);

%Xac dinh V12 va V22
V12=V(1:signals,signals+1:2*signals);
V22=V(signals+1:2*signals,signals+1:2*signals);
disp(V12);
disp(V22);

% Tinh DOA
[eigVectorSi,eigValueSi]=eig(-V12*inv(V22));
disp(eigValueSi);
for i=1:signals
    eigValueSI(i)=eigValueSi(i,i);
end    
%phi=(pi-acos((angle(eigValueSI))/(2*pi*0.5)))*180/pi;
phi=(acos((angle(eigValueSI))/(2*pi*d/lamda)))*180/pi;
disp(phi);

toc
