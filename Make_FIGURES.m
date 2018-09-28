function Make_FIGURES (S_map1, S_map2)

%  Purpose: To map somatomaps or statstical maps on the manikin

% Input Parameters: 
%       S_map1 and S_map2: Matrices with dimensions of 344 by 400. You can 
% input an empty matrix or a matrix of zeros instead of S_map2 in case you
% do not want to map two different activations on one manikin
%        

% Optional Parameters (OP): 
%   First OP : name for the tif file to save 
% 	Second OP : is a number that controls the max value of hot maps
% 	Third OP: you can input a struct file that contains information about
% 	the color maps 


% Example 
% c.map1 = hot(100) ; c.name = 'hot';
% Make_FIGURES (S_map1, [], 'Proportional Map', 15, c)


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%        Author: Mahlega Hassanpour
%        Date : Sat September 29 15:47:14 EDT 2018
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% set defaults
nPx=344;
nPy =400;

if nargin > 3
    file_name = varargin{1};
else 
    file_name = 'Somatomap';
end

if nargin > 4
    max_val = varargin{2};
else 
    max_val = round (max (S_map1(:)));
end

if nargin > 5
    c= varargin{3};
else
    c.name ={'hot', 'winter'};
    c.map1 = hot(100);
    c.map2 = winter(100);
end

f1=figure('units','pixels','Position',[10 100 900 900]);



%% Load the manikin
dir_name = 'C:\Users\User\Desktop\New folder\'; % change this to where the files are savedmask_400,344.mat'
mk_s=imread([dir_name,'white_body_mobile_resized.png']);
load([dir_name,'mask_400,344.mat'],'mask')
 
mk=mk_s;
msk=find(mask==1);
mk2=reshape(mk_s,[],3);
mk2(msk,:)=mk2(msk,:)-50;
mk=reshape(mk2,size(mk,1),size(mk,2),3);
xx=find(mask==0);

%% Make figures


if ~isempty(S_map2) % for mapping Group1>Group2 and vise versa in one map
    S_map1 =S_map1'.*mask;
    S_map2 =S_map2'.*mask;
    max_val = max(round (max (S_map1(:))), round (max (S_map2(:))));
    
    S_map1(1)=max_val;
    S_map2(1)=max_val;
    
    dd=[-max_val  -max_val/2 0 max_val/2 max_val];
    FIGrgb=ActivationOverlay_2Acts(sum(mk,3),'bone',S_map1,S_map2,c.name{1},c.name{2},0.01);
    c.map = [flip(c.map2,1);c.map1];
else
    
    
    S_map1 =S_map1'.*mask;
    S_map1(1)=max_val;
    c.map = [c.map1]; dd=[0 max_val/4 2*max_val/4 3*max_val/4 max_val];
    FIGrgb=ActivationOverlay(sum(mk,3),'bone',S_map1,c.name{1},0.01);
end

% +++ condition 1
image(FIGrgb),title(file_name,'FontSize', 20);


colormap(c.map)
cp3=colorbar ; set (cp3,'Ticks',[0 0.25 0.5 0.75 1 ],'TickLabels',dd);
axis off
axis image;%axis off;
saveas(gcf,[file_name,'.tif'])

end
