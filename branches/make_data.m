% Make data
%huyvq0
[U,fs] = audioread('D:/THUAT_TOAN_MUSIC/file_wav/28-11/test_12k_2811_3.wav');
U1 = U(:,1);
U2 = U(:,2);
%random = -6e-4 + (9e-4)*rand(76704,1);
nhieu_1_part1 = 0.05*U1(30000:33840);
nhieu_2_part1 = 0.05*U2(30000:33840);
nhieu_1_part2 = 0.05*U1(30000:102864);
nhieu_2_part2 = 0.05*U2(30000:102864);
reflex1 = 0.3*U1(30000:30095)-0.001;
reflex2 = 0.3*U2(30000:30095)-0.001;
data1 = [ zeros(10000,1); U1(30000:30095); nhieu_1_part1;reflex1;nhieu_1_part2; U1(30096:30191);nhieu_1_part1;reflex1;nhieu_1_part2;U1(30286:30381);
          nhieu_1_part1;reflex1;nhieu_1_part2;U1(30000:30095); nhieu_1_part1;reflex1;nhieu_1_part2; U1(30096:30191);nhieu_1_part1;reflex1;nhieu_1_part2;
          U1(30286:30381);nhieu_1_part1;reflex1;nhieu_1_part2 ];
data2 = [ zeros(10000,1); U2(30000:30095); nhieu_2_part1;reflex2;nhieu_2_part2; U2(30096:30191);nhieu_2_part1;reflex2;nhieu_2_part2;U2(30286:30381);
          nhieu_2_part1;reflex2;nhieu_2_part2;U2(30000:30095); nhieu_2_part1;reflex2;nhieu_2_part2; U2(30096:30191);nhieu_2_part1;reflex2;nhieu_2_part2;
          U2(30286:30381);nhieu_2_part1;reflex2;nhieu_2_part2 ];
data = [data1 data2];
subplot(2,1,1);
plot(data);
audiowrite('data_2kenh_radar_12k.wav',data,192000);
%Loc 12k
[b,a]=butter(5,[11000,13000]/(fs/2),'bandpass');
filtsig=filter(b,a,data);  %filtered signal
maxkenh1 = max(filtsig(:,1));
maxkenh2 = max(filtsig(:,2));
heso = maxkenh2/maxkenh1; %Tim do lech bien do giua 2 kenh
y1 = heso*filtsig(:,1);
y2 = filtsig(:,2);
y = [ y1 y2 ];%Ghep 2 kenh 
subplot(2,1,2);
plot(y);
