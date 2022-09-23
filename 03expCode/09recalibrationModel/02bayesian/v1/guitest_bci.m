function varargout = guitest_bci(varargin)
% GUITEST_bci MATLAB code for guitest_bci.fig
%      GUITEST_bci, by itself, creates a new GUITEST_bci or raises the existing
%      singleton*.
%
%      H = GUITEST_bci returns the handle to a new GUITEST_bci or the handle to
%      the existing singleton*.
%
%      GUITEST_bci('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUITEST_bci.M with the given input arguments.
%
%      GUITEST_bci('Property','Value',...) creates a new GUITEST_bci or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guitest_bci_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guitest_bci_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guitest_bci

% Last Modified by GUIDE v2.5 23-Sep-2022 14:26:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guitest_bci_OpeningFcn, ...
                   'gui_OutputFcn',  @guitest_bci_OutputFcn, ...
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


% --- Executes just before guitest_bci is made visible.
function guitest_bci_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guitest_bci (see VARARGIN)

% Choose default command line output for guitest_bci
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

% UIWAIT makes guitest_bci wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = guitest_bci_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function text_p_c1_Callback(hObject, eventdata, handles)
% hObject    handle to text_p_c1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_p_c1 as text
%        str2double(get(hObject,'String')) returns contents of text_p_c1 as a double


% --- Executes during object creation, after setting all properties.
function text_p_c1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_p_c1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function resetPlots(handles,plot_idx)
% eval(['axes(handles.axes',num2str(plot_idx),')']);  
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


function text_sigma_c1_Callback(hObject, eventdata, handles)
% hObject    handle to text_sigma_c1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_sigma_c1 as text
%        str2double(get(hObject,'String')) returns contents of text_sigma_c1 as a double


% --- Executes during object creation, after setting all properties.
function text_sigma_c1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_sigma_c1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_sigma_c2_Callback(hObject, eventdata, handles)
% hObject    handle to text_sigma_c2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_sigma_c2 as text
%        str2double(get(hObject,'String')) returns contents of text_sigma_c2 as a double


% --- Executes during object creation, after setting all properties.
function text_sigma_c2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_sigma_c2 (see GCBO)
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



function text_sigma_soa_Callback(hObject, eventdata, handles)
% hObject    handle to text_sigma_soa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_sigma_soa as text
%        str2double(get(hObject,'String')) returns contents of text_sigma_soa as a double


% --- Executes during object creation, after setting all properties.
function text_sigma_soa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_sigma_soa (see GCBO)
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
        'EdgeColor',[0.47,0.67,0.19],'FaceAlpha', 0.5, 'EdgeAlpha', 0.5)

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
function slider_p_c1_Callback(hObject, eventdata, handles)
% hObject    handle to slider_p_c1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_p_c1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_p_c1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slider_sigma_c1_Callback(hObject, eventdata, handles)
% hObject    handle to slider_sigma_c1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_sigma_c1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_sigma_c1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_sigma_c2_Callback(hObject, eventdata, handles)
% hObject    handle to slider_sigma_c2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_sigma_c2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_sigma_c2 (see GCBO)
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
function slider_sigma_soa_Callback(hObject, eventdata, handles)
% hObject    handle to slider_sigma_soa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_sigma_soa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_sigma_soa (see GCBO)
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

% free parameters to be fitted
if get(handles.radiobutton_sliders,'Value') == 1 %read off from the sliders
    param.p_c1          = handles.slider_p_c1.Value;
    param.sigma_c1            = handles.slider_sigma_c1.Value;
    param.sigma_c2            = handles.slider_sigma_c2.Value;
    param.sigma_soa     = handles.slider_sigma_soa.Value;
    param.alpha = handles.slider_alpha.Value;
elseif get(handles.radiobutton_boxes,'Value') == 1 %read off from the text boxes
    param.p_c1          = str2double(get(handles.text_p_c1,'String')); % in sec
    param.sigma_c1             = str2double(get(handles.text_sigma_c1,'String'));
    param.sigma_c2           = str2double(get(handles.text_sigma_c2,'String')); 
    param.sigma_soa     = str2double(get(handles.text_sigma_soa,'String'));
    param.alpha         = str2double(get(handles.text_alpha,'String'));

else
    errordlg('You need to choose between sliders and text boxes!'); 
    return
end

%% initiate recalibration effect for all soas

delta_s               = NaN(n_soas, setup.sim_trial, setup.exposure_trial + 1);
last_recal            = NaN(n_soas, setup.sim_trial); % summarize the last recalibration effect

%% running simulation

for i                 = 1:n_soas

    adaptor_soa                   = adaptor_soas(i);% in s, adaptor_soa, fixed in session

    for t                 = 1:setup.sim_trial

        delta_s = update_recal_bayesian_MA(setup.exposure_trial, adaptor_soa, ...
            param.p_c1, param.sigma_soa, param.sigma_c1, param.sigma_c2, param.alpha);

        last_recal(i,t)     =  - delta_s(:, end);

    end

    %plot the histogram
    eval(['axes(handles.axes',num2str(1+i),');']);
    histogram(last_recal(i,:),20,'FaceColor','k','FaceAlpha',0.3,'EdgeColor','w');
    yticks([]); xticks(round(mean(last_recal(i,:)),4));
end

%% summarize and plot
last_recal_iSOA_error = std(last_recal, [], 2);
last_recal_iSOA = mean(last_recal, 2);
saveResults = {delta_s, last_recal, last_recal_iSOA_error, last_recal_iSOA};

axes(handles.axes1);
set(gca,'FontSize',15,'linewidth',2); hold on; 
errorbar(adaptor_soas, last_recal_iSOA, last_recal_iSOA_error,'.','LineWidth',2); hold off;
yline(0)
xticks(adaptor_soas)
xticklabels(adaptor_soas)
xlabel('adaptor SOA (s)')
ylabel('recalibration (s)')
