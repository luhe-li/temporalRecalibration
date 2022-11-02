function varargout = guitest_MCD(varargin)
% GUITEST_MCD MATLAB code for guitest_MCD.fig
%      GUITEST_MCD, by itself, creates a new GUITEST_MCD or raises the existing
%      singleton*.
%
%      H = GUITEST_MCD returns the handle to a new GUITEST_MCD or the handle to
%      the existing singleton*.
%
%      GUITEST_MCD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUITEST_MCD.M with the given input arguments.
%
%      GUITEST_MCD('Property','Value',...) creates a new GUITEST_MCD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guitest_MCD_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guitest_MCD_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guitest_MCD

% Last Modified by GUIDE v2.5 20-Jul-2022 23:03:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guitest_MCD_OpeningFcn, ...
                   'gui_OutputFcn',  @guitest_MCD_OutputFcn, ...
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


% --- Executes just before guitest_MCD is made visible.
function guitest_MCD_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guitest_MCD (see VARARGIN)

% Choose default command line output for guitest_MCD
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

axes(handles.axes1);
%default figure values
%set(gcf,'color','w');
set(gca,'gridcolor','w');
set(gca,'xminorgrid','off','yminorgrid','off')
set(gca,'gridalpha',1);
set(gca,'LineWidth',1); 
set(gca,'box','off'); % was on
set(gca,'color',[234/256,234/256,238/256]);
set(gca,'TickLength',[0 0]);
set(gca,'fontsize',16); % was 26
xticks([]); yticks([]);

% UIWAIT makes guitest_MCD wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = guitest_MCD_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function text_bias_Callback(hObject, eventdata, handles)
% hObject    handle to text_bias (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_bias as text
%        str2double(get(hObject,'String')) returns contents of text_bias as a double


% --- Executes during object creation, after setting all properties.
function text_bias_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_bias (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function resetPlots(handles,plot_idx)
%eval(['axes(handles.axes',num2str(plot_idx),')']);  
cla reset; 
set(gca,'gridcolor','w');
set(gca,'xminorgrid','off','yminorgrid','off')
set(gca,'gridalpha',1);
set(gca,'LineWidth',1); 
set(gca,'box','off'); % was on
set(gca,'color',[234/256,234/256,238/256]);
set(gca,'TickLength',[0 0]);
set(gca,'fontsize',16); % was 26
xticks([]); yticks([]);
text(0.45,0.5,'Loading...','fontSize',20,'Color',[0.47,0.67,0.19]);


function text_t_V_Callback(hObject, eventdata, handles)
% hObject    handle to text_t_V (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_t_V as text
%        str2double(get(hObject,'String')) returns contents of text_t_V as a double


% --- Executes during object creation, after setting all properties.
function text_t_V_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_t_V (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_t_A_Callback(hObject, eventdata, handles)
% hObject    handle to text_t_A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_t_A as text
%        str2double(get(hObject,'String')) returns contents of text_t_A as a double


% --- Executes during object creation, after setting all properties.
function text_t_A_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_t_A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_t_AV_Callback(hObject, eventdata, handles)
% hObject    handle to text_t_AV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_t_AV as text
%        str2double(get(hObject,'String')) returns contents of text_t_AV as a double


% --- Executes during object creation, after setting all properties.
function text_t_AV_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_t_AV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to text_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_alpha as text
%        str2double(get(hObject,'String')) returns contents of text_alpha as a double


% --- Executes during object creation, after setting all properties.
function text_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_sigma_SOA_Callback(hObject, eventdata, handles)
% hObject    handle to text_sigma_SOA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_sigma_SOA as text
%        str2double(get(hObject,'String')) returns contents of text_sigma_SOA as a double


% --- Executes during object creation, after setting all properties.
function text_sigma_SOA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_sigma_SOA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function addProgressBar(handles, idx)
eval(['axes(handles.axes',num2str(10 + idx),')']);
patch([0, 1, 1, 0], [0, 0, 1, 1], [0.47,0.67,0.19],...
    'EdgeColor',[0.47,0.67,0.19],'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);

function clearProgressBar(handles)
for idx = 1:9
    eval(['axes(handles.axes',num2str(10 + idx),')']);
    patch([0, 1, 1, 0], [0, 0, 1, 1], [0.94,0.94,0.94],...
        'EdgeColor',[0.94,0.94,0.94]);
end


function text_num_simulations_Callback(hObject, eventdata, handles)
% hObject    handle to text_num_simulations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_num_simulations as text
%        str2double(get(hObject,'String')) returns contents of text_num_simulations as a double


% --- Executes during object creation, after setting all properties.
function text_num_simulations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_num_simulations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_ts_duration_Callback(hObject, eventdata, handles)
% hObject    handle to text_ts_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_ts_duration as text
%        str2double(get(hObject,'String')) returns contents of text_ts_duration as a double


% --- Executes during object creation, after setting all properties.
function text_ts_duration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_ts_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_fs_Callback(hObject, eventdata, handles)
% hObject    handle to text_fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_fs as text
%        str2double(get(hObject,'String')) returns contents of text_fs as a double


% --- Executes during object creation, after setting all properties.
function text_fs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_exp_trials_Callback(hObject, eventdata, handles)
% hObject    handle to text_exp_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_exp_trials as text
%        str2double(get(hObject,'String')) returns contents of text_exp_trials as a double


% --- Executes during object creation, after setting all properties.
function text_exp_trials_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_exp_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function radiobutton_boxes_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_boxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value') == 1
    handles.radiobutton_sliders.Value = 0;
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes on button press in radiobutton_sliders.
function radiobutton_sliders_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_sliders (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value') == 1
    handles.radiobutton_boxes.Value = 0;
    % Update handles structure
    guidata(hObject, handles);
end

% --- Executes on slider movement.
function slider_bias_Callback(hObject, eventdata, handles)
% hObject    handle to slider_bias (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_bias_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_bias (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slider_t_V_Callback(hObject, eventdata, handles)
% hObject    handle to slider_t_V (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_t_V_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_t_V (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_t_A_Callback(hObject, eventdata, handles)
% hObject    handle to slider_t_A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_t_A_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_t_A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_t_AV_Callback(hObject, eventdata, handles)
% hObject    handle to slider_t_AV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_t_AV_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_t_AV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to slider_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_sigma_SOA_Callback(hObject, eventdata, handles)
% hObject    handle to slider_sigma_SOA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_sigma_SOA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_sigma_SOA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in radiobutton_boxes.
global setup param saveResults
if get(hObject,'Value')
    SimResults = {setup, param, saveResults};
    save(['SimResults_', datestr(datetime('now')), '.mat'],'SimResults');
    handles.pushbutton_save.Value = 0;
end

% --- Executes on button press in pushbutton_start.
function pushbutton_start_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for i = 1:10; resetPlots(handles,i); end
clearProgressBar(handles)
global setup param saveResults

% set parameters
adaptor_soas = [-700, -300:100:300, 700]./1000; % s
n_soas = length(adaptor_soas);

% set session parameters
setup.sim_trial      = str2double(get(handles.text_num_simulations,'String'));%1;
setup.exposure_trial = str2double(get(handles.text_exp_trials,'String'));%250;

% define fixed parameters to create time series;
% see parameter details in the function below
setup.fs          = str2double(get(handles.text_fs,'String'));%1e3; % hz
stim_dura   = 0.033; % in sec
setup.ts_duration = str2double(get(handles.text_ts_duration,'String')); %16; % in sec, 7 x 2 padding + 2 stimulus duration

% free parameters to be fitted
if get(handles.radiobutton_sliders,'Value') == 1 %read off from the sliders
    param.bias          = handles.slider_bias.Value;
    param.tv            = handles.slider_t_V.Value;
    param.ta            = handles.slider_t_A.Value;
    param.tav           = handles.slider_t_AV.Value;
    param.learning_rate = handles.slider_alpha.Value;
    param.sigma_soa     = handles.slider_sigma_SOA.Value;
elseif get(handles.radiobutton_boxes,'Value') == 1 %read off from the text boxes
    param.bias          = str2double(get(handles.text_bias,'String')); % in sec
    param.tv            = str2double(get(handles.text_t_V,'String'));%0.0873;
    param.ta            = str2double(get(handles.text_t_A,'String')); %0.0684;
    param.tav           = str2double(get(handles.text_t_AV,'String'));%0.7859;
    param.learning_rate = str2double(get(handles.text_alpha,'String'));%0.001;
    param.sigma_soa     = str2double(get(handles.text_sigma_SOA,'String'));%0.2;
else
    errordlg('You need to choose between sliders and text boxes!'); 
    return
end

% initiate recalibration effect for all soas

soa_recal = cell(1, n_soas);
last_recal = NaN(n_soas, setup.sim_trial); % summarize the last recalibration effect

% take MCD parts out from the loop

% initiate time series for the filter
% t should be the same length as the signal time series
nsample = setup.ts_duration * setup.fs + 1; 
t = linspace(0, setup.ts_duration, nsample);

% low-pass filters (Equation 1)
fv=fft(t.*exp(-t./param.tv));
fa=fft(t.*exp(-t./param.ta));
fav=fft(t.*exp(-t./param.tav));

% run simulation

for i = 1:n_soas
    soa = adaptor_soas(i);% in s, adaptor_soa, fixed in session
    %load the progress bar
    addProgressBar(handles, i);
    for t = 1:setup.sim_trial

        soa_recal{i}(t,:) = update_recal_gaussian(setup.exposure_trial, soa, ...
            setup.ts_duration, setup.fs, stim_dura, fa, fv, fav, param.bias, ...
            param.sigma_soa, param.learning_rate);

        last_recal(i, t) = soa_recal{i}(t,end);

    end
    
    %plot the histogram
    eval(['axes(handles.axes',num2str(1+i),');']);
    histogram(last_recal(i,:),20,'FaceColor','k','FaceAlpha',0.3,'EdgeColor','w');
    yticks([]); xticks(round(mean(last_recal(i,:)),4));
end


last_recal_iSOA_error = std(last_recal, [], 2);
last_recal_iSOA = mean(last_recal, 2);
saveResults = {soa_recal, last_recal, last_recal_iSOA_error, last_recal_iSOA};

axes(handles.axes1);
set(gca,'FontSize',15,'linewidth',2); hold on; 
errorbar(adaptor_soas, last_recal_iSOA, last_recal_iSOA_error,'.','LineWidth',2); hold off;
yline(0)
xticks(adaptor_soas)
xticklabels(adaptor_soas)
xlabel('Adaptor SOA (s)')
ylabel('recalibration (s)')
