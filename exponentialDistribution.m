function varargout = exponentialDistribution(varargin)
% EXPONENTIALDISTRIBUTION MATLAB code for exponentialDistribution.fig
%      EXPONENTIALDISTRIBUTION, by itself, creates a new EXPONENTIALDISTRIBUTION or raises the existing
%      singleton*.
%
%      H = EXPONENTIALDISTRIBUTION returns the handle to a new EXPONENTIALDISTRIBUTION or the handle to
%      the existing singleton*.
%
%      EXPONENTIALDISTRIBUTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EXPONENTIALDISTRIBUTION.M with the given input arguments.
%
%      EXPONENTIALDISTRIBUTION('Property','Value',...) creates a new EXPONENTIALDISTRIBUTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before exponentialDistribution_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to exponentialDistribution_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help exponentialDistribution

% Last Modified by GUIDE v2.5 16-Oct-2015 00:12:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @exponentialDistribution_OpeningFcn, ...
                   'gui_OutputFcn',  @exponentialDistribution_OutputFcn, ...
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
end


% --- Executes just before exponentialDistribution is made visible.
function exponentialDistribution_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to exponentialDistribution (see VARARGIN)

% Choose default command line output for exponentialDistribution
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
end
% UIWAIT makes exponentialDistribution wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = exponentialDistribution_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
x1 = str2double(get(handles.edit1, 'String'));
setappdata(0, 'x1', x1);
end
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
x2 = str2double(get(handles.edit2, 'String'));
setappdata(0, 'x2', x2);
end
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
l1 = str2double(get(handles.edit3, 'String'));
setappdata(0, 'l1', l1);
end
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x1 = getappdata(0, 'x1');
x2 = getappdata(0, 'x2');
l1 = getappdata(0, 'l1');
y1 = x1:0.1:x2;

axes(handles.axes1);
if l1 > 0
    e1 = l1 * exp(-abs(l1 * y1));
    plot(y1, e1);
    setappdata(0, 'exponentialDistribution1', e1);
end
end

function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
x3 = str2double(get(handles.edit4, 'String'));
setappdata(0, 'x3', x3);
end
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
x4 = str2double(get(handles.edit5, 'String'));
setappdata(0, 'x4', x4);
end
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
l2 = str2double(get(handles.edit6, 'String'));
setappdata(0, 'l2', l2);
end
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x3 = getappdata(0, 'x3');
x4 = getappdata(0, 'x4');
l2 = getappdata(0, 'l2');
y2 = x3:0.1:x4;

axes(handles.axes2);
if l2 > 0
    e2 = l2 * exp(-abs(l2 * y2));
    plot(y2, e2);
    setappdata(0, 'exponentialDistribution2', e2);
end
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x1 = getappdata(0, 'exponentialDistribution1');
x2 = getappdata(0, 'exponentialDistribution2');

%Input Error handling
if ~isvector(x1) || ~isvector(x2)
   set(handles.edit7, 'String', sprintf('MWWTEST requires vector rather than matrix data.'));
end 
if ~all(isfinite(x1)) || ~all(isnumeric(x1)) || ~all(isfinite(x2)) || ~all(isnumeric(x2))
    set(handles.edit7, 'String', sprintf('Warning: all X1 and X2 values must be numeric and finite'));
end

%set the basic parameter
n1=length(x1); n2=length(x2); NP=n1*n2; N=n1+n2; N1=N+1; k=min([n1 n2]);

[A,B]=tiedrank([x1(:); x2(:)]); %compute the ranks and the ties
R1=A(1:n1); R2=A(n1+1:end); 
T1=sum(R1); T2=sum(R2);
U1=NP+(n1*(n1+1))/2-T1; U2=NP-U1;    
if round(exp(gammaln(N1)-gammaln(k+1)-gammaln(N1-k))) > 20000
    mU=NP/2;
    if B==0
        sU=realsqrt(NP*N1/12);
    else
        sU=realsqrt((NP/(N^2-N))*((N^3-N-2*B)/12));
    end
    Z1=(abs(U1-mU)-0.5)/sU; Z2=(abs(U2-mU)-0.5)/sU; 
    [p,h] = ranksum(x1, x2);
    if B == 0
        set(handles.edit7, 'String', sprintf('Group1 Group2\n\nnumerosity %i %i\nSum of Ranks(W) %0.1f %0.1f\nMean rank %0.1f %0.1f\nTest variable(U) %0.1f %0.1f\nZ corrected for continuity %0.4f %0.4f\nStandard deviation %0.4f\np-value(1-tailed) %0.5f\np-value(2-tailed) %0.5f\nMean %0.1f\nh = %i\nif h = 1(rejection of the null hypothesis)\nif h = 0 (failure to reject the null hypothesis)\nat the 0.5 significance level\nSample size is large enough to use the normal distribution approximation', n1, n2, T1, T2, T1/n1, T2/n2, U1, U2, Z1, Z2, sU, p, 2*p, mU, h));
    else
        set(handles.edit7, 'String', sprintf('Group1 Group2\n\nnumerosity %i %i\nSum of Ranks(W) %0.1f %0.1f\nMean rank %0.1f %0.1f\nTest variable(U) %0.1f %0.1f\nZ corrected for continuity %0.4f %0.4f\nStandard deviation corrected for ties %0.4f\np-value(1-tailed) %0.5f\np-value(2-tailed) %0.5f\nMean %0.1f\nh = %i\nif h = 1(rejection of the null hypothesis)\nif h = 0 (failure to reject the null hypothesis)\nat the 0.5 significance level\nSample size is large enough to use the normal distribution approximation', n1, n2, T1, T2, T1/n1, T2/n2, U1, U2, Z1, Z2, sU, p, 2*p, mU, h));
    end
else
    if n1<=n2
        w=T1;
    else
        w=T2;
    end
    pdf=sum(nchoosek(A,k),2);
    P = [sum(pdf<=w) sum(pdf>=w)]./length(pdf);
    [p,h] = min(P);
    set(handles.edit7, 'String', sprintf('Group1 Group2\n\nnumerosity %i %i\nSum of Ranks(W) %0.1f %0.1f\nMean rank %0.1f %0.1f\nTest variable(U) %0.1f %0.1f\nTest variable(W) %0.1f\np-value(1-tailed) %0.5f\np-value(2-tailed) %0.5f\nh = %i\nif h = 1(rejection of the null hypothesis)\nif h = 0 (failure to reject the null hypothesis)\nat the 0.5 significance level\nSample size is small enough to use the exact Mann-Whitney-Wilcoxon distribution', n1, n2, T1, T2, T1/n1, T2/n2, U1, U2, w, p, 2*p, h));
end
U = [U1 U2];
axes(handles.axes3);
normplot(U);
end


function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
end
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
c = str2double(get(handles.edit8, 'String'));
setappdata(0, 'c', c);
end
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
end
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x1 = getappdata(0, 'exponentialDistribution1');
x2 = getappdata(0, 'exponentialDistribution2');
tst = getappdata(0, 'c');
tail = getappdata(0, 'tail');

global n v alpha
if ~isvector(x1) || ~isvector(x2)
   set(handles.edit9, 'String', sprintf('TESTT requires vector rather than matrix data.'));
end 
if ~all(isfinite(x1)) || ~all(isnumeric(x1)) || ~all(isfinite(x2)) || ~all(isnumeric(x2))
    set(handles.edit9, 'String', sprintf('Warning: all X1 and X2 values must be numeric and finite'));
end

switch tst
    case 0 
        n=[length(x1) length(x2)]; 
        m=[mean(x1) mean(x2)]; 
        v=[var(x1) var(x2)]; 
        if v(2)>v(1) 
            v=fliplr(v);
            m=fliplr(m);
            n=fliplr(n);
        end
        F=v(1)/v(2); 
        DF=n-1;
        p = fcdf(F,DF(1),DF(2)); 
        p = 2*min(p,1-p);
        if p<alpha
            a=v./n; b=sum(a);
            denom=sqrt(b);
            gl=b^2/sum(a.^2./(n-1));
        else
            gl=sum(n)-2; 
            s=sum((n-1).*v)/(sum(n)-2); 
            denom=sqrt(sum(s./n));
        end
        dm=diff(m);
        clear H n m v a b s 
        t=abs(dm)/denom; 
        p=(1-tcdf(t,gl))*tail; 
        h = ttest2(x1,x2);
        [x,str] = powerStudent(t,gl,tail,alpha);
        set(handles.edit9, 'String', sprintf('STUDENT''S T-TEST FOR UNPAIRED SAMPLES\nt-test: %0.5f\nDF: %0.4f\ntail: %d\np-value: %0.4f\nhypothesis: %i\nPower is: %0.4f\n%s', t, gl, tail, p, h, x, str));
    case 1 
        n=length(x1); 
        gl=n-1; 
        d=x1-x2; 
        dm=mean(d); 
        vc=tinv(1-alpha/tail,gl); 
        ic=[abs(dm)-vc abs(dm)+vc]; 
        denom=sqrt((sum((d-dm).^2))/(n*(n-1))); 
        clear n d 
        str=[num2str((1-alpha)*100) '%% C.I.\n'];
        t=abs(dm)/denom; 
        p=(1-tcdf(t,gl))*tail; 
        h = ttest(x1,x2);
        [y,str1] = powerStudent(t,gl,tail,alpha);
        set(handles.edit9, 'String', sprintf('STUDENT''S T-TEST FOR PAIRED SAMPLES\nt-test: %0.5f\ndegrees of freedom: %0.4f\ntail: %d\np-value: %0.4f\nMean of difference: %0.4f\nConfidence interval: %0.4f %0.4f\nhypothesis: %i\nPower is: %0.4f\n%s\n%s', t, gl, tail, p, abs(dm), ic, h, y, str1, str));
end

axes(handles.axes4);
x = [x1 x2];
boxplot(x, 1);
h1 = gca;
h1.XTick = [1 2];
h1.XTickLabel = {'A','B'};
xlabel('X')
ylabel('Y')
end

function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
tail = str2double(get(handles.edit12, 'String'));
setappdata(0, 'tail', tail);
end
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function [x,str] = powerStudent(t,df,c,alpha)
if nargin < 4, 
    alpha = 0.05; %default
end 

if ~isscalar(alpha)
   str = sprintf('POWERSTUDENT requires a scalar ALPHA value.');
end

if ~isnumeric(alpha) || isnan(alpha) || (alpha <= 0) || (alpha >= 1)
   str = sprintf('POWERSTUDENT requires 0 < ALPHA < 1.');
end

if nargin < 3,
    c = 1;
end

if nargin < 2, 
    str = sprintf('Requires at least two input arguments.'); 
end 

t = abs(t);

if c == 1;
   a = alpha;
   P = 1-tcdf(t,df);
   if P >= a;
      tp = tinv(1-a,df) - t;  
      x = 1-tcdf(tp,df);
      str = sprintf('one-tailed hypothesis test\n(The null hypothesis was not statistically significative)\nPower is: %0.4f\nThe associated probability for the Student''s t test is equal or larger than %3.2f\nSo, the assumption of sample means are equal was met.', x, alpha);
   else
      tb = t - tinv(1-a,df); 
      x = tcdf(tb,df);
      str = sprintf('one-tailed hypothesis test\n(The null hypothesis was statistically significative)\nPower is: %0.4f\nThe associated probability for the Student''s t test is smaller than %3.2f\nSo, the assumption of sample means are equal was not met.', x, alpha);
   end
elseif c == 2;
   a = alpha/2;
   P = 1-tcdf(t,df);
   if P >= a;
      tp1 = tinv(1-a,df) - t;  %Power estimation.
      Power1 = 1-tcdf(tp1,df);
      tp2 = t + tinv(1-a,df);
      Power2 = 1-tcdf(tp2,df);
      x = Power1 + Power2;
      str = sprintf('two-tailed hypothesis test\n(The null hypothesis was not statistically significative)\nPower is: %0.4f\nThe associated probability for the Student''s t test is equal or larger than %3.2f\nSo, the assumption of sample means are equal was met.', x, alpha);      
   else
      tb1 = t - tinv(1-a,df); 
      b1 = 1-tcdf(tb1,df);
      tb2 = t + tinv(1-a,df);
      b2 = 1-tcdf(tb2,df);
      x = 1 - (b1 - b2);
      str = sprintf('two-tailed hypothesis test\n(The null hypothesis was statistically significative)\nPower is: %0.4f\nThe associated probability for the Student''s t test is smaller than %3.2f\nSo, the assumption of sample means are equal was not met.', x, alpha);
   end
end
return,
end
