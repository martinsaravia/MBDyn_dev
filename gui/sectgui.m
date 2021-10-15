function varargout = sectgui(varargin)
% MBDYNhandles M-file for sectgui.fig
%      MBDYNhandles, by itself, creates a new MBDYNhandles or raises the existing
%      singleton*.
%
%      H = MBDYNhandles returns the handle to a new MBDYNhandles or the handle to
%      the existing singleton*.
%
%      MBDYNhandles('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MBDYNhandles.M with the given input arguments.
%
%      MBDYNhandles('Property','Value',...) creates a new MBDYNhandles or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the handles before sectgui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sectgui_OpeningFcn via varargin.
%
%      *See handles Options on handlesDE's Tools menu.  Choose "handles allows only one
%      instance to run (singleton)".
%
% See also: handlesDE, handlesDATA, handlesHANDLES

% Edit the above text to modify the response to help sectgui

% Last Modified by handlesDE v2.5 04-Apr-2013 18:30:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sectgui_OpeningFcn, ...
                   'gui_OutputFcn',  @sectgui_OutputFcn, ...
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






% --- Executes just before sectgui is made visible.
function sectgui_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for sectgui
handles.output = hObject;
handles.guiname = get(handles.figure1,'Name');

handles.ploth(3) = plot(handles.axes3,0,0);   handles.axesh(3) = handles.axes3;
handles.ploth(4) = plot(handles.axes4,0,0);   handles.axesh(4) = handles.axes4;

% Update handles structure
guidata(hObject, handles);








% --- Outputs from this function are returned to the command line.
function varargout = sectgui_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles;  % Ouputs the handles to the command line



% --- Executes on button press in runbutton.
function runbutton_Callback(hObject, eventdata, handles)
set(handles.sectable,'Data',zeros(9,9))

sectkernel

set(handles.sectable,'Data',D)



function fileedit_Callback(hObject, eventdata, handles)

handles.filename = get(hObject,'String');

set (handles.figure1, 'Name', [  'mbDyn V2.17 - ' handles.filename]); % Change the name of the figure to add filename 

guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of fileedit as text
%        str2double(get(hObject,'String')) returns contents of fileedit as a double


% --- Executes during object creation, after setting all properties.
function fileedit_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in stopbutton.
function stopbutton_Callback(hObject, eventdata, handles)

disp ('----------   Cancelling Run   -------------')

set(hObject,'UserData',1);  % Set userdata to 1 in order to set the flago for cancelling the run (note that the GUI interacts with the model directly)

guidata(hObject,handles); % Update the object properties



% --- Executes on button press in pausebutton.
function pausebutton_Callback(hObject, eventdata, handles)


if get(hObject,'UserData')==1
    set(hObject,'UserData',0);
    set(hObject,'String','Pause');
    disp ('----------   Resuming Run...  -------------')
else
    set(hObject,'UserData',1);
    disp ('----------   Pausing Run...  -------------')
    set(hObject,'String','Cont.');
    pause on
    pause
end

guidata(hObject,handles); % Update the object properties






% --- Executes on button press in openfile.
function openfile_Callback(hObject, eventdata, handles)

[name dire] = uigetfile ({'*.inps';'*.sinp'},'Open File'); % Open Dialog

[pathstr, name, ext, versn] = fileparts(name);

if name == 0; return; end

set(handles.fileedit,'String',[dire name])
set (handles.figure1, 'Name', [ 'mbDyn V2.17' ' - ' dire name]); % Change the name of the figure to add filename 


handles.dire = dire;
handles.exte = ext;
handles.file = name;

guidata(hObject, handles); % Updates Handles
