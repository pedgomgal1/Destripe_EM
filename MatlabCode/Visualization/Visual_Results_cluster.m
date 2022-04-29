clc; clear; close all;

build = true;
Results_path = ['Results',filesep,'V2_TVL_proj'];

if build
tic
files = dir([Results_path,filesep,'R*.mat']);

for file_number = 1:length(files)
	load([Results_path,filesep,files(file_number).name]);
	disp(files(file_number).name);
	lam_X(file_number) = lambda_X;
	lam_Z(file_number) = lambda_Ztild;
	Res{file_number} = Z;
end
Res{length(lam_X)+1} = double(slices(1:512,1:512,:))/255;

clearvars('-except','Res','lam_Z','lam_X','Results_path');
toc
% save([Results_path,filesep,'Res_Visualization.mat']);
else 
	load([Results_path,filesep,'Res_Visualization.mat']);
end

for i = 1:length(lam_X)
lambda_string(i) =  string(['TV: ',num2str(lam_X(i)),' L: ',num2str(lam_Z(i))]); 
end
lambda_string = [lambda_string,"original"];

global result_number_1
global result_number_2
global slice
result_number_1 = length(lam_X)+1;
result_number_2 = 1;
slice = 1;
[m,n,z] = size(Res{end});


fig1 = figure;
ax(1)=subplot(1,2,1); im1 = imagesc(Res{result_number_1}(:,:,slice)); colormap gray; axis image; title('image1');%caxis([0,255]);
ax(2)=subplot(1,2,2); im2 = imagesc(real(Res{result_number_2}(:,:,slice))); colormap gray; axis image; title('image2');%caxis([0,255]);
linkaxes(ax,'xy');

drawnow;

txt = uicontrol(fig1,...
        'Style','text',...
        'String','slice 1',...
        'Units', 'Normalized',...
        'Position', [0.2 0.05 0.6 0.05]);



sld = uicontrol(fig1,...
        'Style', 'slider',...
        'Min',1,'Max',z,'Value',1,...
        'Units', 'Normalized',...
        'Position', [0.2 0.0 0.6 0.05],...
        'Callback', {@fun2, Res, txt, im1, im2},...
        'SliderStep', [1/(z-1) 1]);

pop_up_1 = uicontrol(fig1,...
		'Style', 'popup',...
		'Value', result_number_1,...
   		'String', lambda_string,...
   		'Units', 'Normalized',...
   		'Callback', {@fun3, im1, Res },...
   		'Position', [0.18 0.1 0.25 0.05]);
txt_l = uicontrol(fig1,...
        	'Style','text',...
        	'String','image1',...
        	'Units', 'Normalized',...
        	'Position', [0.08 0.1 0.1 0.05]);

pop_up_2 = uicontrol(fig1,...
		'Style', 'popup',...
   		'String', lambda_string,...
   		'Units', 'Normalized',...
   		'Callback', {@fun4, im2, Res },...
   		'Position', [0.6 0.1 0.25 0.05]);
txt_2 = uicontrol(fig1,...
        	'Style','text',...
        	'String','image2',...
        	'Units', 'Normalized',...
        	'Position', [0.5 0.1 0.1 0.05]);

function fun2(hObject, eventdata, Res, txt, im1, im2)
global result_number_1
global result_number_2
global slice
slice = floor(get(hObject,'Value'));
im1.CData = Res{result_number_1}(:,:,slice);
im2.CData = Res{result_number_2}(:,:,slice);
txt.String = ['slice ',num2str(slice)];
end

function fun3(hObject, eventdata, im1, Res)
global result_number_1
global result_number_2
global slice
result_number_1 = get(hObject,'Value');
im1.CData = Res{result_number_1}(:,:,slice);
end

function fun4(hObject, eventdata, im2, Res)
global result_number_1
global result_number_2
global slice
result_number_2 = get(hObject,'Value');
im2.CData = Res{result_number_2}(:,:,slice);
end

