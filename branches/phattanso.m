function phattanso(tanso,thoigian)
%t=0:2*pi/(192000/tanso):3000000*pi;
%x=sin(t);
%x=x(1:192000*thoigian);
%
somau=1:1:192000*thoigian;
x=sin(2*pi*(tanso/192000)*somau);
x=x/2+0.5;
size(x)
plot(x) 
sound([x' x'],192000);