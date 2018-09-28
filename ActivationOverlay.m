function FIGrgb=ActivationOverlay(bknd,bkgndMap,act,overlayMap,cutoff)

% Purpose: This function overlays functional activations on anatomical backgrounds.
% Separate colormaps are used for the background and the activation.



% Input Parameters: 
%       bknd = background image (typically an anitomical image)
%       bkgndMap = colormap for background, typically bone
%       act = activation image
%       overlayMap = colormap for activations, typically hot
%       cutoff = threshold for activation defined as a percent of maximum
%           activation [0 1.00]
%


%
%
% Output Parameters: 

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%        Author: Mahlega Hassanpour
%        Date : Sat September 08 15:47:14 EDT 2018
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%% 0) sanity check to see if background and activation images are the same size
if (size(bknd) ~= size(act))
    disp('Error in ActivationOverlay: images are not the same size')
    status = -1;
    return
end

%% 1) Set parameters
res=100; % resolution of colormaps
cMap1=eval([bkgndMap,'(',num2str(res),');']);
cMap2=eval([overlayMap,'(',num2str(res),');']);
cMap3=flipdim(cool(100),1); % default for negative activations is cool (inverted)
step=1/res;
cutoffP=round(res*cutoff);
PLOT=0;

%% 2) Normalize so that doubles can define the colormap
bknd=bknd./max(abs(bknd(:)));
act=act./max(act(:));
actN=act./min(act(:));
cutoffN=round(0.5*res);

%% 3) Set new values from colormaps
FIGrgb=zeros(numel(bknd),3);

% Background
for i=1:res
    Bcells=intersect(find(bknd>=step*(i-1)),find(bknd<=step*(i)));
    for j=1:3
        FIGrgb(Bcells,j)=cMap1(i,j);
    end
end

% Activation
for i=1:res
    AcellsPos=[];AcellsNeg=[];
    if i>cutoffP
        AcellsPos=intersect(find(act>=step*(i-1)),find(act<=step*(i)));
    end
%     if i>cutoffN
%         AcellsNeg=intersect(find(actN>=step*(i-1)),find(actN<=step*(i)));
%     end
    for j=1:3
        FIGrgb(AcellsPos,j)=cMap2(i,j);
%         FIGrgb(AcellsNeg,j)=cMap3(i,j);
    end
end
FIGrgb=reshape(FIGrgb,size(bknd,1),size(bknd,2),3);

%% Plot
if PLOT==1
    figure;
    image(FIGrgb);
    axis image
    axis off
end