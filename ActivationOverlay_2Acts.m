function FIGrgb=ActivationOverlay_2Acts(bknd,bkgndMap,actA,actB,overlayMapA,overlayMapB,cutoff)

% Purpose: This function overlays functional activations on anatomical backgrounds.
% Separate colormaps are used for the background and the activation.



% Input Parameters: 
%       bknd = background image (typically an anitomical image)
%       bkgndMap = colormap for background, typically bone
%       actA and actB = two activation images
%       overlayMapA and overlayMapB= colormaps for activations A and B
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
if (size(bknd) ~= size(actA))
    disp('Error in ActivationOverlay: images are not the same size')
    status = -1;
    return
end

%% 1) Set parameters
res=100; % resolution of colormaps
cMap1=eval([bkgndMap,'(',num2str(res),');']);

colors={'red','green','blue','cyan','magenta','yellow'};
if any(strcmp(overlayMapA,colors))
    cMap2=GenerateHotMap(overlayMapA,res);
else
    cMap2=eval([overlayMapA,'(',num2str(res),');']);
end
if any(strcmp(overlayMapB,colors))
    cMap3=GenerateHotMap(overlayMapB,res);
else
    cMap3=eval([overlayMapB,'(',num2str(res),');']);
end

step=1/res;
cutoffP=round(res*cutoff);
PLOT=0;

%% 2) Normalize so that doubles can define the colormap
bknd=bknd./max(abs(bknd(:)));
actA=actA./max(actA(:));
actB=actB./max(actB(:));

%% 3) Set new values from colormaps
FIGrgb=zeros(numel(bknd),3);

% Background
for i=1:res
    Bcells=intersect(find(bknd>=step*(i-1)),find(bknd<=step*(i)));
    for j=1:3
        FIGrgb(Bcells,j)=cMap1(i,j);
    end
end

% Darken background anatomy
beta=-0.25;
FIGrgb=FIGrgb.^(1/(1+beta));


% Activation
for i=1:res
    AcellsPosA=[];
    AcellsPosB=[];
    if i>cutoffP
        AcellsPosA=intersect(find(actA>=step*(i-1)),find(actA<=step*(i)));
        AcellsPosB=intersect(find(actB>=step*(i-1)),find(actB<=step*(i)));
    end
    for j=1:3
        FIGrgb(AcellsPosA,j)=cMap2(i,j);
        FIGrgb(AcellsPosB,j)=cMap3(i,j);
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