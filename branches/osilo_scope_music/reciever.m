function varargout = reciever(varargin)

% Author Nguyen Quoc Khuong
% Modify by Huy Vuong
% Ha noi University of Sience and Technology
% Radio electronic and Teclecommunication School
% mail:  khuong.nguyenquoc@mail.hust.vn


% This version  using downsample technique


% RECIEVER MATLAB code for reciever.fig
%      RECIEVER, by itself, creates a new RECIEVER or raises the existing
%      singleton*.
%
%      H = RECIEVER returns the handle to a new RECIEVER or the handle to
%      the existing singleton*.
%
%      RECIEVER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RECIEVER.M with the given input arguments.
%
%      RECIEVER('Property','Value',...) creates a new RECIEVER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before reciever_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      reset.  All inputs are passed to reciever_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help reciever

% Last Modified by GUIDE v2.5 09-Apr-2021 11:00:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @reciever_OpeningFcn, ...
                   'gui_OutputFcn',  @reciever_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);


if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT


% --- Executes just before reciever is made visible.
function reciever_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to reciever (see VARARGIN)

% Choose default command line output for reciever
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes reciever wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = reciever_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
 
% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
 
cla(handles.axes3);
cla(handles.axes13);
cla(handles.axes14);
clc;
set(handles.disp_txt,'String','Show text');
receiver_parametter;
 
 set(handles.mod_method,'String',  R_method);
set(handles.edit_FS,'String',num2str(192000));
set(handles.edit_F0,'String',num2str(12400));
set(handles.edit_df,'String',num2str(400));
set(handles.edit_bitrates,'String',num2str(100));


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3



function edit_F0_Callback(hObject, eventdata, handles)
% hObject    handle to edit_F0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_F0 as text
%        str2double(get(hObject,'String')) returns contents of edit_F0 as a double


% --- Executes during object creation, after setting all properties.
function edit_F0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_F0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in mod_method.
function mod_method_Callback(hObject, eventdata, handles)
% hObject    handle to mod_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns mod_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mod_method

warning off
val = get(handles.mod_method,'Value');
Bang_tan = get(handles.popupmenu3,'Value');
Chon_so_kenh=get(handles.popupmenu4,'Value');
Filter_Option = get(handles.popupmenu2,'Value');

flm=get(handles.edit_FS,'String');       % Number of QASK modulation level
flm=str2double(flm); 

do_rong_xung=get(handles.edit46,'String');       % Number of QASK modulation level
do_rong_xung=str2double(do_rong_xung); 

chu_ky_xung=get(handles.edit47,'String');       % Number of QASK modulation level
chu_ky_xung=str2double(chu_ky_xung); 

Somau_chuky=round(flm*chu_ky_xung);

Hien_thi_nhieu_chu_ky=get(handles.checkbox28,'Value');

bac_bo_loc=get(handles.edit43,'String');       
bac_bo_loc=str2double(bac_bo_loc); 

Tan_so_duoi=get(handles.edit44,'String');       
Tan_so_duoi=str2double(Tan_so_duoi); 
Tan_so_tren=get(handles.edit45,'String');       
Tan_so_tren=str2double(Tan_so_tren); 
 
if Bang_tan==1
switch Filter_Option
    case 2
        [Bf Af]=butter(bac_bo_loc,2*Tan_so_tren/flm);
 
    case 3
        [Bf Af]=butter(bac_bo_loc,[2*Tan_so_duoi/flm 2*Tan_so_tren/flm]);
     
    case 4
        [Bf Af]=butter(bac_bo_loc,2*Tan_so_tren/flm,'high');
    
    case 5
         Bf=fir1(bac_bo_loc,2*Tan_so_tren/flm);
   
end
else
   switch Bang_tan
    case 2
        [Bf Af]=butter(bac_bo_loc,[5/48 5/24]);
    case 3
        [Bf Af]=butter(bac_bo_loc,[5/24 5/16]);
     case 4
        [Bf Af]=butter(bac_bo_loc,[5/16 5/12]);
     case 5
        [Bf Af]=butter(bac_bo_loc,[5/12 1/2]);
     case 6
        [Bf Af]=butter(bac_bo_loc,[45/96 55/96]);
     case 7
        [Bf Af]=butter(bac_bo_loc,[55/96 65/96]);
     case 8
        [Bf Af]=butter(bac_bo_loc,[65/96 75/96]);
      case 9
        [Bf Af]=butter(bac_bo_loc,[75/96 85/96]);
       case 10
        [Bf Af]=butter(bac_bo_loc,[3/8 11/24]);   
end

end

Tin_hieu_goc_1=get(handles.channel1,'Value');
Tin_hieu_goc_2=get(handles.channel2,'Value');
Pho_Tin_hieu_goc_1=get(handles.checkbox19,'Value');
Pho_Tin_hieu_goc_2=get(handles.checkbox20,'Value');

Tin_hieu_loc_1=get(handles.checkbox21,'Value');
Tin_hieu_loc_2=get(handles.checkbox22,'Value');
Pho_Tin_hieu_loc_1=get(handles.checkbox25,'Value');
Pho_Tin_hieu_loc_2=get(handles.checkbox26,'Value');

 

system_parametter;
global r; 
r = audiorecorder(flm, 16, 1);
stop(r);
 rt=get(handles.edit_rectime,'String');
rt=str2double(rt);
Count_down=rt;  
switch val
        case 2
             r = audiorecorder(flm, 16, 2);
            record(r);
            tttt = clock;    
            while etime(clock, tttt)<rt
                set(handles.edit50,'String',num2str(Count_down-etime(clock, tttt)));
                pause(0.01);
            end
             
            stop(r);
            y =  getaudiodata(r,'double'); 
             switch Chon_so_kenh
                 case 1 
                axes(handles.axes3);
                plot(-y);  
                set(handles.axes3,'YLim',[-2 2]);
                 Y=abs(fft(y));
                [mm nn]=size(y)
                axes(handles.axes13);
                plot(Y);  
                 case 2
                     y=y(:,1);
                     axes(handles.axes3);
                    plot(-y);  
                    set(handles.axes3,'YLim',[-2 2]);
                     Y=abs(fft(y));
                    [mm nn]=size(y)
                    axes(handles.axes13);
                    plot(Y);  
                 case 3
                     y=y(:,2);
                     axes(handles.axes3);
                    plot(-y);  
                    set(handles.axes3,'YLim',[-2 2]);
                     Y=abs(fft(y));
                    [mm nn]=size(y)
                    axes(handles.axes13);
                    plot(Y); 
                  
            end
    end            
%
while val==3   %Osilo 1 kenh
         
    Magnitude=get(handles.osilo_mag,'String');
    Magnitude=str2double(Magnitude);
    if rt<0.05 rt=0.05;
    end
        r = audiorecorder(flm, 16, 1);
    record(r);
 
%record(handles.r2);
 
    while val==3
            pause(rt);
            stop(r);
            y =  getaudiodata(r,'double' ); 
             
            record(r); 
                mm=length(y);
           
            axes(handles.axes3);
        
            plot(-Magnitude*y);  
         set(handles.axes3,'YLim',[-2 2]);
       
       Y=abs(fft(y));
       
       axes(handles.axes13);
 
        plot(Y);  
       [ii dd]=max(y(mm/40:mm/2));
        set(handles.disp_txt,'String',['Main Frequency:  ' num2str((mm/40+dd-2)*192000/mm)]);
        
        if Filter_Option>1
         yf=filter(Bf,Af,y);
         YF=abs(fft(yf));
         axes(handles.axes14);
         plot(-yf)
          axes(handles.axes15);
         plot(YF)
        
        
        if Hien_thi_nhieu_chu_ky==1
            N=floor(mm/Somau_chuky)
            yy=yf(1:N*Somau_chuky);
            size(yy)
            yy=reshape(yy,Somau_chuky,N);
            yy=sum(yy,2);
            axes(handles.axes21);
            plot(yy)

        end
        end
        
        val = get(handles.mod_method,'Value');
      
    end
   stop(r);

end

%%%%%%%%%% End of 3

 while val==4   %Osilo 2 kenh
         
    Magnitude=get(handles.osilo_mag,'String');
    Magnitude=str2double(Magnitude);
  
    if rt<0.05 rt=0.05;
    end
    
        r = audiorecorder(flm, 16, 2);
   
record(r);
 
%record(handles.r2);
    
    while val==4
 
           pause(rt);
            stop(r);
            y =  getaudiodata(r,'double' ); 
             
            record(r); 
             switch Chon_so_kenh
                 case 1 
                axes(handles.axes3);
                yy=y;
                yy(:,1)=yy(:,1)+1;
                yy(:,2)=yy(:,2)-1;
                plot(-Magnitude*yy);  
                set(handles.axes3,'YLim',[-2 2]);

               Y=abs(fft(y));
               [mm, nn]=size(y);
               axes(handles.axes13);

                plot(Y);  
               [ii dd]=max(y(mm/40:mm/2));
                set(handles.disp_txt,'String',['Main Frequency:  ' num2str((mm/40+dd-2)*192000/mm)]);

                if Filter_Option>1
                 yf=filter(Bf,Af,y);
                 YF=abs(fft(yf));
                 axes(handles.axes14);
                 yy=yf;
                    yy(:,1)=yy(:,1)+max(yy(:,1));;
                    yy(:,2)=yy(:,2)-max(yy(:,2));;
                 plot(-yy)
                  axes(handles.axes15);
                 plot(YF)
                end
        
                if Hien_thi_nhieu_chu_ky==1
                    N=floor(mm/Somau_chuky);
                    yy=yf(1:N*Somau_chuky,:);

                        yy1=reshape(yy(:,1),Somau_chuky,N);
                          yy1=sum(yy1,2);
                        yy2=reshape(yy(:,2),Somau_chuky,N);
                          yy2=sum(yy2,2);
                          yy=[yy1 yy2];

                    axes(handles.axes21);
                    plot(yy)
        

                end
                 case 2
                     y=y(:,1);
                           mm=length(y);
           
            axes(handles.axes3);
        
            plot(-Magnitude*y);  
         set(handles.axes3,'YLim',[-2 2]);
       
       Y=abs(fft(y));
       
       axes(handles.axes13);
 
        plot(Y);  
       [ii dd]=max(y(mm/40:mm/2));
        set(handles.disp_txt,'String',['Main Frequency:  ' num2str((mm/40+dd-2)*192000/mm)]);
        
        if Filter_Option>1
         yf=filter(Bf,Af,y);
         YF=abs(fft(yf));
         axes(handles.axes14);
         plot(-yf)
          axes(handles.axes15);
         plot(YF)
        
        
        if Hien_thi_nhieu_chu_ky==1
            N=floor(mm/Somau_chuky)
            yy=yf(1:N*Somau_chuky);
            size(yy)
            yy=reshape(yy,Somau_chuky,N);
            yy=sum(yy,2);
            axes(handles.axes21);
            plot(yy)

        end
        end
                     
                 case 3
                     y=y(:,1);
                           mm=length(y);
           
            axes(handles.axes3);
        
            plot(-Magnitude*y);  
         set(handles.axes3,'YLim',[-2 2]);
       
       Y=abs(fft(y));
       
       axes(handles.axes13);
 
        plot(Y);  
       [ii dd]=max(y(mm/40:mm/2));
        set(handles.disp_txt,'String',['Main Frequency:  ' num2str((mm/40+dd-2)*192000/mm)]);
        
        if Filter_Option>1
         yf=filter(Bf,Af,y);
         YF=abs(fft(yf));
         axes(handles.axes14);
         plot(-yf)
          axes(handles.axes15);
         plot(YF)
        
        
        if Hien_thi_nhieu_chu_ky==1
            N=floor(mm/Somau_chuky)
            yy=yf(1:N*Somau_chuky);
            size(yy)
            yy=reshape(yy,Somau_chuky,N);
            yy=sum(yy,2);
            axes(handles.axes21);
            plot(yy)

        end
        end
                     
             end
        val = get(handles.mod_method,'Value');
      
    end
   stop(r);

end

%  End of 5
%  Start of Video rec
 

% --- Executes during object creation, after setting all properties.
function mod_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mod_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_rectime_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rectime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rectime as text
%        str2double(get(hObject,'String')) returns contents of edit_rectime as a double


% --- Executes during object creation, after setting all properties.
function edit_rectime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rectime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_data.
function save_data_Callback(hObject, eventdata, handles)
% hObject    handle to save_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
flm=get(handles.edit_FS,'String');       % Number of QASK modulation level
flm=str2double(flm); 
 
 
l= get(handles.axes3, 'Children'); 
y=get(l,'Ydata');
[mm nn]=size(y);
    if mm==2   
        y=cell2mat(y)';
    end
 File_name=get(handles.File_name,'String');
     
        
val = get(handles.mod_method,'Value');

cc=clock;
dt=[ num2str(cc(3))  num2str(cc(4)) num2str(cc(5)) num2str(round(cc(6)))];
 
name=strcat(File_name,dt,'.wav'); 
size(y)
audiowrite(y,flm,16,name);
  
% --- Executes on button press in load_data.
function load_data_Callback(hObject, eventdata, handles)
% hObject    handle to load_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname]=uigetfile({'*wav'});
path=[pathname,filename];

 
if (filename==0)
    set(handles.load_data,'String','Cannot load');
    pause(1);
    set(handles.load_data,'String','Load wav file');
 %    uicontrol(hObject)
else
    %y=wavread(path);
    [y,flm1]=audioread(path); 

%====================== fix in 13:13 24/03
% huy_f = 12000;
% huy_v = 1500;
% huy_lamda=huy_v/huy_f;
% huy_M = 2; %so phan tu mang thu
% huy_gain = 0;
% huy_k=2*pi/huy_lamda;
% %huy_D=1; %so nguon tin hieu
% huy_d = huy_lamda/2;
% [huy_b,huy_a]=butter(5,[11000,13000]/(flm1/2),'bandpass');
% huy_filtsig=filter(huy_b,huy_a,y);  %filtered signal
% maxkenh1 = max(huy_filtsig(:,1));
% maxkenh2 = max(huy_filtsig(:,2));
% heso2kenh = maxkenh2/maxkenh1; %Tim do lech bien do giua 2 kenh
% data_kenh1 = heso2kenh*huy_filtsig(:,1);
% data_kenh2 = huy_filtsig(:,2);
% ghep2kenh = [ data_kenh1 data_kenh2 ];%Ghep 2 kenh 
% hilbert2kenh = hilbert(ghep2kenh);% Chuyen tu real => phuc
% huy_somau = length(ghep2kenh); %so mau
% dao_hilbert = hilbert2kenh';
% time_lm = 1/flm1;
% Ruu=dao_hilbert*dao_hilbert'/huy_somau;
% [eigVector,eigValue]=eig(Ruu);
% huy_signals = 1;
% %Xac dinh cac vector rieng cua khong gian nhieu
% eigVectorNoise=eigVector(:,1:huy_M-huy_signals);
% 
% %Pho khong gian cua tin hieu
% huy_i=0; 
% for theta=0:.1:180
%     huy_i=huy_i+1; 
%     A0_tmp=10^(huy_gain/10)*exp(j*huy_k*(0:huy_M-1)*huy_d*(cos(theta*pi/180))); 
%     A0=A0_tmp.';
%     P(huy_i)=10*log((A0'*A0)/(A0'*eigVectorNoise*eigVectorNoise'*A0))/10;
%     %P(i)=((A0'*A0)/(A0'*eigVectorNoise*eigVectorNoise'*A0))
%     [maxP,huy_index] = max(real(P));
% end
% result_DOA = huy_index/10;
% display(result_DOA);
%======================
axes(handles.axes3);
plot(y);  
axes(handles.axes13);
plot(abs(fft(y)));  

 
flm=get(handles.edit_FS,'String');       % Number of QASK modulation level
flm=str2double(flm); 
 
 
val = get(handles.mod_method,'Value')

%set(handles.disp_txt,'String',['Goc toi    :  ' num2str(result_DOA)]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
 
  
 % --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
 %get(handles.checkbox1,'Value')

 


 
 
   
function edit22_Callback(hObject, eventdata, handles)
 
 function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fc_Callback(hObject, eventdata, handles)
 function edit_fc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


 


function edit_mod_Callback(hObject, eventdata, handles)
 function edit_mod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function No_datasubs_Callback(hObject, eventdata, handles)
 function No_datasubs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to No_datasubs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_pro_Callback(hObject, eventdata, handles)
 function edit_pro_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


 
function osilo_mag_Callback(hObject, eventdata, handles)
 function osilo_mag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to osilo_mag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

  
  
    
 
% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit43_Callback(hObject, eventdata, handles)
% hObject    handle to edit43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit43 as text
%        str2double(get(hObject,'String')) returns contents of edit43 as a double


% --- Executes during object creation, after setting all properties.
function edit43_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox19.
function checkbox19_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox19


% --- Executes on button press in checkbox20.
function checkbox20_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox20


% --- Executes on button press in checkbox21.
function checkbox21_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox21


% --- Executes on button press in checkbox22.
function checkbox22_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox22


% --- Executes on button press in checkbox25.
function checkbox25_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox25


% --- Executes on button press in checkbox26.
function checkbox26_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox26


% --- Executes on button press in checkbox27.



function edit44_Callback(hObject, eventdata, handles)
% hObject    handle to edit44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit44 as text
%        str2double(get(hObject,'String')) returns contents of edit44 as a double


% --- Executes during object creation, after setting all properties.
function edit44_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit45_Callback(hObject, eventdata, handles)
% hObject    handle to edit45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit45 as text
%        str2double(get(hObject,'String')) returns contents of edit45 as a double


% --- Executes during object creation, after setting all properties.
function edit45_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function File_name_Callback(hObject, eventdata, handles)
% hObject    handle to File_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of File_name as text
%        str2double(get(hObject,'String')) returns contents of File_name as a double


% --- Executes during object creation, after setting all properties.
function File_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to File_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Filter.
function Filter_Callback(hObject, eventdata, handles)
% hObject    handle to Filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
flm=get(handles.edit_FS,'String');       % Number of QASK modulation level
flm=str2double(flm); 

Bang_tan = get(handles.popupmenu3,'Value');

do_rong_xung=get(handles.edit46,'String');       % Number of QASK modulation level
do_rong_xung=str2double(do_rong_xung); 

chu_ky_xung=get(handles.edit47,'String');       % Number of QASK modulation level
chu_ky_xung=str2double(chu_ky_xung); 

X= get(handles.axes3, 'XLim')
X=round(X)
 


Hien_thi_nhieu_chu_ky=get(handles.checkbox28,'Value');
 
l= get(handles.axes3, 'Children'); 
y=get(l,'Ydata');
[mm nn]=size(y)
    if mm==2   
        y=cell2mat(y)';
        [nn mm]=size(y)
    end
if X(2)>length(y) X(2)=length(y); 
 end
  
 if X(1)<=0   X(1)=1;
 end
  
 if mm==2
     y=y(X(1):X(2),:);
 else
     y=y(X(1):X(2))';
 end
 
    Filter_Option = get(handles.popupmenu2,'Value');

flm=get(handles.edit_FS,'String');       % Number of QASK modulation level
flm=str2double(flm); 
bac_bo_loc=get(handles.edit43,'String');       
bac_bo_loc=str2double(bac_bo_loc); 

Tan_so_duoi=get(handles.edit44,'String');       
Tan_so_duoi=str2double(Tan_so_duoi); 
Tan_so_tren=get(handles.edit45,'String');       
Tan_so_tren=str2double(Tan_so_tren); 
if Bang_tan==1
 
switch Filter_Option
    case 2
        [Bf Af]=butter(bac_bo_loc,2*Tan_so_tren/flm);
        axes(handles.axes20);
        [F P]=freqz(Bf);
        plot(abs(F))
    case 3
        [Bf Af]=butter(bac_bo_loc,[2*Tan_so_duoi/flm 2*Tan_so_tren/flm]);
        axes(handles.axes20);
        [F P]=freqz(Bf);
        plot(abs(F))
    case 4
        [Bf Af]=butter(bac_bo_loc,2*Tan_so_duoi/flm,'high');
        axes(handles.axes20);
        [F P]=freqz(Bf);
        plot(abs(F))
    case 5
         Bf=fir1(bac_bo_loc,2*Tan_so_tren/flm);
         axes(handles.axes20);
        [F P]=freqz(Bf);
        plot(abs(F))
         case 6
         Bf=fir1(bac_bo_loc,[2*Tan_so_duoi/flm 2*Tan_so_tren/flm]);
         axes(handles.axes20);
        [F P]=freqz(Bf);
        plot(abs(F))
         case 7
         Bf=fir1(bac_bo_loc,2*Tan_so_tren/flm,'high');
         axes(handles.axes20);
        [F P]=freqz(Bf);
        plot(abs(F))
end
 

else
   switch Bang_tan
   
   case 2
        [Bf Af]=butter(bac_bo_loc,[5/48 5/24]);
    case 3
        [Bf Af]=butter(bac_bo_loc,[5/24 5/16]);
     case 4
        [Bf Af]=butter(bac_bo_loc,[5/16 5/12]);
     case 5
        [Bf Af]=butter(bac_bo_loc,[5/12 1/2]);
     case 6
        [Bf Af]=butter(bac_bo_loc,[45/96 55/96]);
     case 7
        [Bf Af]=butter(bac_bo_loc,[55/96 65/96]);
     case 8
        [Bf Af]=butter(bac_bo_loc,[65/96 75/96]);
      case 9
        [Bf Af]=butter(bac_bo_loc,[75/96 85/96]);
       case 10
        [Bf Af]=butter(bac_bo_loc,[3/8 11/24]);   
end

end
if Filter_Option>1
    if Filter_Option>5
        if mm==1 y=conv(y,Bf);
        else
           y1=conv(Bf,y(:,1)); 
           y2=conv(Bf,y(:,2));
           
           y=[y1 y2];
        end
    else
        y=filter(Bf,Af,y);
    end
end
%Chon_so_kenh=get(handles.checkbox27,'Value');
apdungthuattoan=get(handles.checkbox30,'Value');  % Thuat toan cong tin hieu tren 1 kenh cho truong hop khong dong bo phat thu
SMCK=round(flm*chu_ky_xung);
SMDRX=2*round(flm*do_rong_xung)
RC=30;  % Khoang quet
axes(handles.axes14);
if mm==1
         plot(-y)
else
    yy=y;
    yy(:,1)=yy(:,1)+max(yy(:,1));
    yy(:,2)=yy(:,2)-max(yy(:,2)); 
    plot(-yy)
end
 axes(handles.axes15);
         plot(abs(fft(y)))
  % y=y(15000:end);      nn=length(y);
  [nn mm]=size(y)
  if Hien_thi_nhieu_chu_ky==1 
      N=floor(nn/SMCK)
           
            if mm==1
                if apdungthuattoan==0
                    yy=y(1:N*SMCK);
                    yy=reshape(yy,SMCK,N);
                      yy=sum(yy,2);
                else
                 yy=y(1:N*SMCK);
                yy=reshape(yy,SMCK,N);
                [MX  MD]=max(yy(:,1))
                Ketqua=[];
                GTChuan=yy(MD-SMDRX:MD+SMDRX,1);
                for ii=2:N
                    SS=[];
                    for kk=-RC:1:RC
                        gt=sum(GTChuan-yy(MD-SMDRX+kk:MD+SMDRX+kk,ii));
                        SS=[SS; gt];
                    end
                    Ketqua=[Ketqua  SS];
                end
                [NX ND]=min(Ketqua)
                SHIFT=ND-RC
                yyy=yy(RC+2:end-RC-2,1);
                L=length(yyy);
                for ii=2:N
                        yyy=[yyy yy(RC+SHIFT(ii-1):L+RC+SHIFT(ii-1)-1,ii)];
                end
                  yy=sum(yyy,2);
                end
              
                
            else
                if apdungthuattoan==0
                yy=y(1:N*SMCK,:);
                yy1=reshape(yy(:,1),SMCK,N);
                   yy1=sum(yy1,2);
                yy2=reshape(yy(:,2),SMCK,N);
                   yy2=sum(yy2,2);
                  yy=[yy1 yy2];
                else
                      yy=y(1:N*SMCK,:);
                yy1=reshape(yy(:,1),SMCK,N);
                 yy2=reshape(yy(:,2),SMCK,N);
                 %%%%%%%%
                 size(y)
                 [MX  MD]=max(yy1(:,1));
                 [MX2  MD2]=max(yy2(:,1));
                Ketqua=[];
                Ketqua2=[];
                GTChuan=yy1(MD-SMDRX:MD+SMDRX,1);
                GTChuan2=yy2(MD-SMDRX:MD+SMDRX,1);
                for ii=2:N
                    SS=[]; SS2=[];
                    for kk=-RC:1:RC
                        gt=sum(GTChuan-yy1(MD-SMDRX+kk:MD+SMDRX+kk,ii));
                        SS=[SS; gt];
                        gt2=sum(GTChuan2-yy2(MD-SMDRX+kk:MD+SMDRX+kk,ii));
                        SS2=[SS2; gt2];
                    end
                    Ketqua=[Ketqua  SS];
                    Ketqua2=[Ketqua2  SS2];
                end
                [NX ND]=min(Ketqua);
                [NX2 ND2]=min(Ketqua2);
                SHIFT=ND-RC;
                SHIFT2=ND2-RC;
                yyy1=yy1(RC+2:end-RC-2,1);
                yyy2=yy2(RC+2:end-RC-2,1);
                L=length(yyy1);
                for ii=2:N
                        yyy1=[yyy1 yy1(RC+SHIFT(ii-1):L+RC+SHIFT(ii-1)-1,ii)];
                        yyy2=[yyy2 yy2(RC+SHIFT2(ii-1):L+RC+SHIFT2(ii-1)-1,ii)];
                end
                %%%%%%%%%%%%%%
                 
                  yy=[sum(yyy1,2) sum(yyy2,2)];
                  
                end
            end
            size(yy)
axes(handles.axes21);
if mm==1
           plot(yy)
else
    yyy=yy;
    yyy(:,1)=yyy(:,1)+max(yyy(:,1));;
    yyy(:,2)=yyy(:,2)-max(yyy(:,2));; 
    plot(yyy)
end
         % plot(yy(:,1:3))
       %   axes(handles.axes22);
        %  plot(sum(yyy,2))  
  end
Combination2Channel=get(handles.checkbox29,'Value');
if Combination2Channel==1
        yy=sum(yy,2);
    axes(handles.axes22);
         plot(yy)

end

function edit46_Callback(hObject, eventdata, handles)
% hObject    handle to edit46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit46 as text
%        str2double(get(hObject,'String')) returns contents of edit46 as a double


% --- Executes during object creation, after setting all properties.
function edit46_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit47_Callback(hObject, eventdata, handles)
% hObject    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit47 as text
%        str2double(get(hObject,'String')) returns contents of edit47 as a double


% --- Executes during object creation, after setting all properties.
function edit47_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox28.
function checkbox28_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox28


% --- Executes on button press in checkbox29.
function checkbox29_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox29



function edit48_Callback(hObject, eventdata, handles)
% hObject    handle to edit48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit48 as text
%        str2double(get(hObject,'String')) returns contents of edit48 as a double


% --- Executes during object creation, after setting all properties.
function edit48_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit49_Callback(hObject, eventdata, handles)
% hObject    handle to edit49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit49 as text
%        str2double(get(hObject,'String')) returns contents of edit49 as a double


% --- Executes during object creation, after setting all properties.
function edit49_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox30.
function checkbox30_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox30

 
% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
flm=get(handles.edit_FS,'String');       % Number of QASK modulation level
flm= (flm); 
X= get(handles.axes21, 'XLim');
X=round(X);
KC=(X(2)-X(1))*1500/flm/2;
set(handles.disp_txt,'String',['Khoang cach:  ' num2str(KC)]);
%set(handles.disp_txt,'String',['Goc toi    :  ' num2str(result_DOA)]);
%%=========fix this

% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
flm=get(handles.edit_FS,'String');       % Number of QASK modulation level
flm=str2double(flm); 
X= get(handles.axes14, 'XLim');
X=round(X);
KC=(X(2)-X(1))*1500/flm/2
set(handles.disp_txt,'String',['Khoang cach:  ' num2str(KC)]);
%set(handles.disp_txt,'String',['Goc toi: 39']);
%===== fix this

% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
flm=get(handles.edit_FS,'String');       % Number of QASK modulation level
flm=str2double(flm); 

Bang_tan = get(handles.popupmenu3,'Value');

do_rong_xung=get(handles.edit46,'String');       % Number of QASK modulation level
do_rong_xung=str2double(do_rong_xung); 

chu_ky_xung=get(handles.edit47,'String');       % Number of QASK modulation level
chu_ky_xung=str2double(chu_ky_xung); 

X= get(handles.axes3, 'XLim')
X=round(X)
 


Hien_thi_nhieu_chu_ky=get(handles.checkbox28,'Value');
 
l= get(handles.axes3, 'Children'); 
y=get(l,'Ydata');
[mm nn]=size(y)
    if mm==2   
        y=cell2mat(y)';
        [nn mm]=size(y)
    end
if X(2)>length(y) X(2)=length(y); 
 end
  
 if X(1)<=0   X(1)=1;
 end
  
 if mm==2
     y=y(X(1):X(2),:);
 else
     y=y(X(1):X(2))';
 end
 
 axes(handles.axes14);
         plot(y)
         Y=abs(fft(y));

         axes(handles.axes15);
         plot(Y)
          [mm, nn]=size(Y)
       axes(handles.axes13);
 
        plot(Y);  
       [ii dd]=max(Y(200:mm/2))
        set(handles.disp_txt,'String',['Main Frequency:  ' num2str((200+dd-2)*192000/mm)]);
        %set(handles.disp_txt,'String',['Goc :  ' num2str(result_DOA)]); % alo alo

% --- Executes on button press in checkbox33.
 
% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4
Chon_so_kenh=get(handles.popupmenu4,'Value');
switch Chon_so_kenh
    case 1
    set(handles.channel1,'Value',1,'Enable','On','Visible','On');
     set(handles.checkbox19,'Value',1,'Enable','On','Visible','On');
      set(handles.checkbox21,'Value',1,'Enable','On','Visible','On');
       set(handles.checkbox25,'Value',1,'Enable','On','Visible','On');
       set(handles.channel2,'Value',1,'Enable','On','Visible','On');
     set(handles.checkbox20,'Value',1,'Enable','On','Visible','On');
      set(handles.checkbox22,'Value',1,'Enable','On','Visible','On');
       set(handles.checkbox26,'Value',1,'Enable','On','Visible','On');
    case 2
      set(handles.channel1,'Value',1,'Enable','Off','Visible','On');
     set(handles.checkbox19,'Value',1,'Enable','Off','Visible','On');
      set(handles.checkbox21,'Value',1,'Enable','Off','Visible','On');
       set(handles.checkbox25,'Value',1,'Enable','Off','Visible','On');
       set(handles.channel2,'Value',0,'Enable','On','Visible','Off');
     set(handles.checkbox20,'Value',0,'Enable','On','Visible','Off');
      set(handles.checkbox22,'Value',0,'Enable','On','Visible','Off');
       set(handles.checkbox26,'Value',0,'Enable','On','Visible','Off');
    case 3
       set(handles.channel1,'Value',0,'Enable','On','Visible','Off');
     set(handles.checkbox19,'Value',0,'Enable','On','Visible','Off');
      set(handles.checkbox21,'Value',0,'Enable','On','Visible','Off');
       set(handles.checkbox25,'Value',0,'Enable','On','Visible','Off');
       set(handles.channel2,'Value',1,'Enable','Off','Visible','On');
     set(handles.checkbox20,'Value',1,'Enable','Off','Visible','On');
      set(handles.checkbox22,'Value',1,'Enable','Off','Visible','On');
       set(handles.checkbox26,'Value',1,'Enable','Off','Visible','On');
end

% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit50_Callback(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit50 as text
%        str2double(get(hObject,'String')) returns contents of edit50 as a double


% --- Executes during object creation, after setting all properties.
function edit50_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nut_tinh_goc.
function nut_tinh_goc_Callback(hObject, eventdata, handles)
% hObject    handle to nut_tinh_goc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
flm=get(handles.edit_FS,'String');       % Number of QASK modulation level
flm=str2double(flm); 

Bang_tan = get(handles.popupmenu3,'Value');

do_rong_xung=get(handles.edit46,'String');       % Number of QASK modulation level
do_rong_xung=str2double(do_rong_xung); 

chu_ky_xung=get(handles.edit47,'String');       % Number of QASK modulation level
chu_ky_xung=str2double(chu_ky_xung); 

X= get(handles.axes3, 'XLim')
X=round(X)
 


Hien_thi_nhieu_chu_ky=get(handles.checkbox28,'Value');
 
l= get(handles.axes3, 'Children'); 
y=get(l,'Ydata');
[mm nn]=size(y)
    if mm==2   
        y=cell2mat(y)';
        [nn mm]=size(y)
    end
if X(2)>length(y) X(2)=length(y); 
 end
  
 if X(1)<=0   X(1)=1;
 end
  
 if mm==2
     y=y(X(1):X(2),:);
 else
     y=y(X(1):X(2))';
 end
 
huy_f = 12000;
huy_v = 1500;
huy_lamda=huy_v/huy_f;
huy_M = 2; %so phan tu mang thu
huy_gain = 0;
huy_k=2*pi/huy_lamda;
huy_D=1; %so nguon tin hieu
huy_d = huy_lamda/2;
flm1 = 48000;
[huy_b,huy_a]=butter(5,[11000,13000]/(flm1/2),'bandpass');
huy_filtsig=filter(huy_b,huy_a,y);  %filtered signal
maxkenh1 = max(huy_filtsig(:,1));
maxkenh2 = max(huy_filtsig(:,2));
heso2kenh = maxkenh2/maxkenh1; %Tim do lech bien do giua 2 kenh
data_kenh1 = heso2kenh*huy_filtsig(:,1);
data_kenh2 = huy_filtsig(:,2);
ghep2kenh = [ data_kenh1 data_kenh2 ];%Ghep 2 kenh 
hilbert2kenh = hilbert(ghep2kenh);% Chuyen tu real => phuc
huy_somau = length(ghep2kenh); %so mau
dao_hilbert = hilbert2kenh';
time_lm = 1/flm1;
huy_t = 0:time_lm:(length(ghep2kenh)*time_lm)-time_lm;
Ruu=dao_hilbert*dao_hilbert'/huy_somau;
[eigVector,eigValue]=eig(Ruu);
eigValueMax=max(max(eigValue));
%Xac dinh so nguon tin hieu den
%huy_signals=length(find(diag(eigValue)>eigValueMax/1000000));
huy_signals = 1;

%Xac dinh cac vector rieng cua khong gian nhieu
eigVectorNoise=eigVector(:,1:huy_M-huy_signals);

%Pho khong gian cua tin hieu
huy_i=0; 
 for theta=0:.1:180
    huy_i=huy_i+1; 
    A0_tmp=10^(huy_gain/10)*exp(j*huy_k*(0:huy_M-1)*huy_d*(cos(theta*pi/180))); 
    A0=A0_tmp.';
    P(huy_i)=10*log((A0'*A0)/(A0'*eigVectorNoise*eigVectorNoise'*A0))/10;
    %P(i)=((A0'*A0)/(A0'*eigVectorNoise*eigVectorNoise'*A0))
    [maxP,huy_index] = max(real(P));
 end
 result_DOA = huy_index/10;
 display(result_DOA);
set(handles.disp_txt,'String',['Goc : ',num2str(result_DOA)]);
