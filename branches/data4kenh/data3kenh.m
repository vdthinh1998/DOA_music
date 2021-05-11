[y, fs]=audioread('CV_TN_80_165_200kHz_2.wav');
frame=length(y)/50;
subplot(2,1,1);
plot(y);
[Bf1 Af1]=butter(5,[5/12 1/2]);%40kHz----->50kHz
[Bf2 Af2]=butter(5,[55/96 65/96]);%55kHz----->65kHz
[Bf3 Af3]=butter(5,[75/96 85/96]);%75kHz----->85kHz
y1=filter(Bf1,Af1,y);% Loc y lay y1=45k
y2=filter(Bf2,Af2,y);%Loc y lay y2=60k
y3=filter(Bf3,Af3,y);%Loc y lay y2=80k
y1=y1/(max(y1));
y2=y2/(max(y2));
y3=y3/(max(y3));
y1 = y1(1:576000);
y2 = y2(1:576000);
y3 = y3(1:576000);
%subplot(2,1,2);
%plot(y1);
dt = 1/fs;                   % seconds per sample
StopTime = 3;             % seconds
t = (0:dt:StopTime-dt)';     % seconds
%%Sine wave:
Fc1 = 33000;                     % hertz
x1 = cos(2*pi*Fc1*t);
y_nhan1 = y1.*x1;
%plot(y_nhan);
[b,a]=butter(5,[11000,13000]/(fs/2),'bandpass');
y_12k_1 = filter(b,a,y_nhan1);  %filtered signal
figure(2)
plot(y_12k_1);
title('tin hieu 12k');
%================== 
%Xu li kenh 1
%==================
y_12k_1 = abs(y_12k_1);
  y1_enveloped=envelope(y_12k_1,150,'rms'); %lay mau 150 mau/lan
  y1_enveloped=y1_enveloped/(max(y1_enveloped));
  y1_enveloped = y1_enveloped(1:80000);
  figure(3)
  plot(y1_enveloped);
  xlabel('Time(s)');
  ylabel('Amplitude(V)');
  title('ENVELOPE SIGNAL');

 i = 1;
 stt = 0;
 D1=[];
while i<=(frame)
    switch stt
        case 0
            if(y1(i)>= 0.4) 
                index1 = i;
                i = i+400;
                stt = 1;
            else
                i = i+1;
                stt = 0;
            end
        case 1
            if(y1(i)>= max(y1(i:i+16000)))
                index2 = i;
                i = i+400;
                time = index2-index1;
                distance=time*1500/2;
                disp('Distance 1(m)');
                distance=distance*1/fs;
                disp(distance);
                stt = 0;
            else
                i = i+1;
                stt = 1;
            end
        otherwise
            i = i+1;
            stt = 0;
    end
   %D1=[D1 distance];
end
