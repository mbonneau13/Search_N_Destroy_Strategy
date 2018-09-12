function varargout = GUI_melaleuca(varargin)
% GUI_MELALEUCA MATLAB code for GUI_melaleuca.fig
%      GUI_MELALEUCA, by itself, creates a new GUI_MELALEUCA or raises the existing
%      singleton*.
%
%      H = GUI_MELALEUCA returns the handle to a new GUI_MELALEUCA or the handle to
%      the existing singleton*.
%
%      GUI_MELALEUCA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_MELALEUCA.M with the given input arguments.
%
%      GUI_MELALEUCA('Property','Value',...) creates a new GUI_MELALEUCA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_melaleuca_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_melaleuca_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_melaleuca

% Last Modified by GUIDE v2.5 06-Sep-2018 10:49:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_melaleuca_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_melaleuca_OutputFcn, ...
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


% --- Executes just before GUI_melaleuca is made visible.
function GUI_melaleuca_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_melaleuca (see VARARGIN)

% Choose default command line output for GUI_melaleuca
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
%%% Plot the UF logo 
axes(handles.axes1)
X = imread('UF.jpg');
imshow(X)

%%% Plot the USGS logo 
axes(handles.axes2)
X = imread('USGS.jpg');
imshow(X)

%% Background

% UIWAIT makes GUI_melaleuca wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_melaleuca_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in RUN.
function RUN_Callback(hObject, eventdata, handles)
% hObject    handle to RUN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
time = tic;
h = waitbar(0,'Please wait...');
oldpointer = get(handles.figure1, 'pointer');
set(handles.figure1, 'pointer', 'watch')
drawnow;

if exist('NbRow.mat','file') == 2
    
    load('NbRow.mat')
else
    
    errordlg('Don''t forget to enter the number of rows')
    return
end

if exist('NbCol.mat','file') == 2
    
    load('NbCol.mat')
else
    
    errordlg('Don''t forget to enter the number of columns')
    return
end

ZN = zeros(ligne,colonne);

if exist('Data.mat','file') == 2
    
    load('Data.mat')
else
    
    errordlg('Don''t forget to enter the input data')
    return
end

Inside = Data(:,1);
ZN(Inside) = ones;

Search = -ones(ligne,colonne);
Search(Inside) = Data(:,3);

O = zeros(ligne,colonne); 
O(Inside) = Data(:,2);

Destroy = -ones(ligne,colonne);
Destroy(Inside) = Data(:,end);

load('Budget.mat')

if exist('alpha1.mat','file') == 2
    
    load('alpha1.mat')
    delete('alpha1.mat')
else
    
    alpha1 = [-0.5 , -0.5];
end

if exist('beta1.mat','file') == 2
    
    load('beta1.mat')
    delete('beta1.mat')
else
    
    beta1 = 0.5;
end

if exist('Theta1.mat','file') == 2
    
    load('Theta1.mat')
    delete('Theta1.mat')
else
    
    Theta1 = [0.1 , 0.3];
end

if exist('nTest.mat','file') == 2
    
    load('nTest.mat')
    delete('nTest.mat')
else
    
    nTest = 5;
end

if exist('max_iter.mat','file') == 2
    
    load('max_iter.mat')
    delete('max_iter.mat')
else
    
    max_iter = 100;
end

if exist('Nit.mat','file') == 2
    
    load('Nit.mat')
    delete('Nit.mat')
else
    
    Nit = 10;
end

Run_FunctionGUI(O,ZN,Search,Destroy,Budget,alpha1,beta1,Theta1,nTest,max_iter,Nit,ligne,colonne)
delete('NbCol.mat')
delete('NbRow.mat')
delete('Data.mat')
delete('Budget.mat')
disp(['Job was done in ' num2str(toc(time)) ' secondes !'])


function OutputName_Callback(hObject, eventdata, handles)
% hObject    handle to OutputName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OutputName as text
%        str2double(get(hObject,'String')) returns contents of OutputName as a double

Out = hObject.String;
save('Out.mat','Out')


% --- Executes during object creation, after setting all properties.
function OutputName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OutputName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Input_Data_Callback(hObject, eventdata, handles)
% hObject    handle to Input_Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Input_Data as text
%        str2double(get(hObject,'String')) returns contents of Input_Data as a double
try
    Name = hObject.String;
    i = find(Name == '.');
    if isempty(i)
        errordlg('Please, write the file extension (csv or xlsx)')
    else
        if isequal(Name(i:end),'.xlsx')
            
            Data = xlsread(Name);
        else
            Data = csvread(Name);
        end
        if ~isequal(size(Data,2),4)
            errordlg('Data should have 4 columns (Inside sites number, Observation value, Search cost, Destroy cost)')
        else
            save('Data.mat','Data')
        end
    end
catch
    errordlg(['No file ' hObject.String ' is present in the folder, type another name.' ])
end


% --- Executes during object creation, after setting all properties.
function Input_Data_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Input_Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function E_Destroy_Callback(hObject, eventdata, handles)
% hObject    handle to E_Destroy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of E_Destroy as text
%        str2double(get(hObject,'String')) returns contents of E_Destroy as a double
try
    Name = hObject.String;
    i = find(Name == '.');
    if isempty(i)
        errordlg('Please, write the file extension (csv or xlsx)')
    else
        if isequal(Name(i:end),'.xlsx')
            
            Destroy = xlsread(Name);
        else
            Destroy = csvread(Name);
        end
        [l , c] = size(Destroy);
        if c > l
            Destroy = Destroy';
        end
        save('Destroy.mat','Destroy')
    end
catch
    errordlg(['No file ' hObject.String 'is present in the folder, type another name.' ])
end


% --- Executes during object creation, after setting all properties.
function E_Destroy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to E_Destroy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function E_Obs_Callback(hObject, eventdata, handles)
% hObject    handle to E_Obs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of E_Obs as text
%        str2double(get(hObject,'String')) returns contents of E_Obs as a double
try
    Name = hObject.String;
    i = find(Name == '.');
    if isempty(i)
        errordlg('Please, write the file extension (csv or xlsx)')
    else
        if isequal(Name(i:end),'.xlsx')
            
            Obs = xlsread(Name);
        else
            Obs = csvread(Name);
        end
        [l , c] = size(Obs);
        if c > l
            Obs = Obs';
        end
        save('Obs.mat','Obs')
    end
catch
    errordlg(['No file ' hObject.String 'is present in the folder, type another name.' ])
end


% --- Executes during object creation, after setting all properties.
function E_Obs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to E_Obs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function E_Search_Callback(hObject, eventdata, handles)
% hObject    handle to E_Search (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of E_Search as text
%        str2double(get(hObject,'String')) returns contents of E_Search as a double
try
    Name = hObject.String;
    i = find(Name == '.');
    if isempty(i)
        errordlg('Please, write the file extension (csv or xlsx)')
    else
        if isequal(Name(i:end),'.xlsx')
            
            Search = xlsread(Name);
        else
            Search = csvread(Name);
        end
        [l , c] = size(Search);
        if c > l
            Search = Search';
        end
        save('Search.mat','Search')
    end
catch
    errordlg(['No file ' hObject.String 'is present in the folder, type another name.' ])
end

% --- Executes during object creation, after setting all properties.
function E_Search_CreateFcn(hObject, eventdata, handles)
% hObject    handle to E_Search (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nRow_Callback(hObject, eventdata, handles)
% hObject    handle to nRow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nRow as text
%        str2double(get(hObject,'String')) returns contents of nRow as a double
try
    ligne = str2double(hObject.String);
    save('NbRow.mat','ligne')
catch
    errordlg('The number of row should be an integer')
end

% --- Executes during object creation, after setting all properties.
function nRow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nRow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nCol_Callback(hObject, eventdata, handles)
% hObject    handle to nCol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nCol as text
%        str2double(get(hObject,'String')) returns contents of nCol as a double
try
    colonne = str2double(hObject.String);
    save('NbCol.mat','colonne')
catch
    errordlg('The number of column should be an integer')
end

% --- Executes during object creation, after setting all properties.
function nCol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nCol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function budget_Callback(hObject, eventdata, handles)
% hObject    handle to budget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of budget as text
%        str2double(get(hObject,'String')) returns contents of budget as a double
try
    Budget = str2double(hObject.String);
    save('Budget.mat','Budget')
catch
    errordlg('Problem witht he budget.')
end

% --- Executes during object creation, after setting all properties.
function budget_CreateFcn(hObject, eventdata, handles)
% hObject    handle to budget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function Alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Beta_Callback(hObject, eventdata, handles)
% hObject    handle to Beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Beta as text
%        str2double(get(hObject,'String')) returns contents of Beta as a double
try
    beta1 = str2double(hObject.String);
    save('beta1.mat','beta1')
catch
    errordlg('Problem witht the beta value')
end


% --- Executes during object creation, after setting all properties.
function Beta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function theta_Callback(hObject, eventdata, handles)
% hObject    handle to theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of theta as text
%        str2double(get(hObject,'String')) returns contents of theta as a double
try
    V = hObject.String;
    
    Theta1 = str2double(hObject.String);
    save('Theta1.mat','Theta1')
catch
    errordlg('Problem witht the Theta value')
end

% --- Executes during object creation, after setting all properties.
function theta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Alpha_Callback(hObject, eventdata, handles)
% hObject    handle to Alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Alpha as text
%        str2double(get(hObject,'String')) returns contents of Alpha as a double
try
    alpha1 = str2double(hObject.String);
    save('alpha1.mat','alpha1')
catch
    errordlg('Problem witht the alpha value')
end


function nTest_Callback(hObject, eventdata, handles)
% hObject    handle to nTest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nTest as text
%        str2double(get(hObject,'String')) returns contents of nTest as a double
try
    nTest = str2double(hObject.String);
    save('nTest.mat','nTest')
catch
    errordlg('Problem witht the nTest value')
end


% --- Executes during object creation, after setting all properties.
function nTest_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nTest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max_iter_Callback(hObject, eventdata, handles)
% hObject    handle to max_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_iter as text
%        str2double(get(hObject,'String')) returns contents of max_iter as a double
try
    max_iter = str2double(hObject.String);
    save('max_iter.mat','max_iter')
catch
    errordlg('Problem witht the max_iter value')
end


% --- Executes during object creation, after setting all properties.
function max_iter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Nit_Callback(hObject, eventdata, handles)
% hObject    handle to Nit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nit as text
%        str2double(get(hObject,'String')) returns contents of Nit as a double
try
    Nit = str2double(hObject.String);
    save('Nit.mat','Nit')
catch
    errordlg('Problem witht the Nit value')
end

% --- Executes during object creation, after setting all properties.
function Nit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in HelpNbRow.
function HelpNbRow_Callback(hObject, eventdata, handles)
% hObject    handle to HelpNbRow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('You have to define the number of rows in the study area.','Help','help')

% --- Executes on button press in HelpColumn.
function HelpColumn_Callback(hObject, eventdata, handles)
% hObject    handle to HelpColumn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('You have to define the number the columns in the study area.','Help','help')


% --- Executes on button press in HelpBudget.
function HelpBudget_Callback(hObject, eventdata, handles)
% hObject    handle to HelpBudget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('You have to define the Search and Destroy budget here. It should be an integer number','Help','help')

% --- Executes on button press in HelpBeta.
function HelpBeta_Callback(hObject, eventdata, handles)
% hObject    handle to HelpBeta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('You have to define the initial value of the beta parameter. This intial value is used by the EM algorithl to estimate the model parameters. If you have no idea of what can be the initial value, just leave it as it is.','Help','help')


% --- Executes on button press in HelpAlpha.
function HelpAlpha_Callback(hObject, eventdata, handles)
% hObject    handle to HelpAlpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('You have to define the initial value of the alpha parameter. This intial value is used by the EM algorithl to estimate the model parameters. Be carefull the values have to be written between square brackets and separated by a comma. The first value is alpha_0 and the second alpha_1. If you have no idea of what can be the initial value, just leave it as it is.','Help','help')


% --- Executes on button press in HelpTheta.
function HelpTheta_Callback(hObject, eventdata, handles)
% hObject    handle to HelpTheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('You have to define the initial value of the theta parameter. This intial value is used by the EM algorithl to estimate the model parameters. Be carefull the values have to be written between square brackets and separated by a comma. The first value is the False positive rate and is the false negative rate. If you have no idea of what can be the initial value, just leave it as it is.','Help','help')


% --- Executes on button press in HelpnTest.
function HelpnTest_Callback(hObject, eventdata, handles)
% hObject    handle to HelpnTest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('You have to define the initial value of the nTest parameter. This is the number of run used by the EM algorithm. If you have no idea of what can be the initial value, just leave it as it is.','Help','help')

% --- Executes on button press in HelpMaxIter.
function HelpMaxIter_Callback(hObject, eventdata, handles)
% hObject    handle to HelpMaxIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('You have to define the initial value of the max_iter parameter. This is the maximal number of iteration of the EM algorithm for a specific run. If you have no idea of what can be the initial value, just leave it as it is.','Help','help')

% --- Executes on button press in HelpNit.
function HelpNit_Callback(hObject, eventdata, handles)
% hObject    handle to HelpNit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('You have to define the initial value of the Nit parameter. This is the number of steps for the Gibbs Sampling algorithm. If you have no idea of what can be the initial value, just leave it as it is.','Help','help')

% --- Executes on button press in HelpOutput.
function HelpOutput_Callback(hObject, eventdata, handles)
% hObject    handle to HelpOutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('Just put the filename where you want to have the results. The results will be a list of site numbers to be searched.','Help','help')


% --- Executes on button press in HelpInput.
function HelpInput_Callback(hObject, eventdata, handles)
% hObject    handle to HelpInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('Specify the filename and extension (ex data.csv) of the data input. The file should be an xlsx or csv file. The first column should contain the site number that are in the study area. The second column should have the observation value. The third column should have the Seach cost and the last column the destroy cost.','Help','help')
