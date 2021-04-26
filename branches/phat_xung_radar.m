function phat_xung_radar(tanso,time_len)
% Phien ban 1 tan so co gan tin hieu ma hoa thoi gian va thu tu phat xung

chu_ky=0.2;
do_rong_xung=0.0002;
flm=192000;
tim=clock
thoigian_phat=tim(4)*3600+tim(5)*60+round(tim(6));
tg_bin=de2bi(thoigian_phat,18,'left-msb');
TG=[];
num_of_pulse=floor(time_len/chu_ky);
t2=0:2*pi/20:1000*pi;
f2=sin(t2);
do_dai_header=10000;
f2=f2(1:do_dai_header);
for ii=1:num_of_pulse
    xx=[1 1 tg_bin de2bi(ii,10,'left-msb') zeros(1,do_dai_header/20-30) ];
    xx=kron(ones(20,1),xx);
    xx=xx(:);
    xx=xx.*f2';
    TG=[TG; xx'];
end
  TG=TG/5;
DRX=round(do_rong_xung*flm);
CK=round(flm*chu_ky)-do_dai_header       % 1000 is Header len
x=[ones(1,DRX) zeros(1,CK-DRX)];
t=0:2*pi/(flm/tanso):50000*pi;

f0=sin(t);
f0=f0(1:CK);
y=x.*f0;
  
yy=[];
for i=1:num_of_pulse
    yg=[TG(i,:) y];
    yy=[yy yg];
end
   plot(yy)
  %  wavwrite(yy,192000,16,'xung_Phat_gan_thoi_gian.wav');
% plot(abs(fft(yy)))
 sound([yy' yy'],192000)