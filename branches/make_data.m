% Make data
%huyvq0
[U,fs] = audioread('D:/THUAT_TOAN_MUSIC/file_wav/28-11/test_12k_2811_3.wav');
U1 = U(:,1);
random = -6e-4 + (9e-4)*rand(76704,1);
data1 = [ zeros(10000,1); U1(30000:30095); random; U1(30096:30191);-6e-4 + (9e-4)*rand(76704,1);U1(30286:30381);
          -6e-4 + (9e-4)*rand(76704,1);U1(30000:30095); -6e-4 + (9e-4)*rand(76704,1); U1(30096:30191);zeros(76704,1);U1(30286:30381)
    ];

subplot(2,1,1);
plot(data1);

%Loc 12k
[b,a]=butter(5,[11000,13000]/(fs/2),'bandpass');
filtsig=filter(b,a,U);  %filtered signal
maxkenh1 = max(filtsig(:,1));
maxkenh2 = max(filtsig(:,2));
heso = maxkenh2/maxkenh1; %Tim do lech bien do giua 2 kenh
y1 = heso*filtsig(:,1);
y2 = filtsig(:,2);
y = [ y1 y2 ];%Ghep 2 kenh 
subplot(2,1,2);
plot(y);