function [C_th_map]= Cluster_ANALYSIS(z_map,C_S_th,varargin)

% Purpose: To cluster size threshold a given map at the given size 


% Input Parameters: 
%       z_map : statistical map, it is a vector 1 by number of voxels
%       C_S_th : is the threshold for the cluster size, the defult can be 
% pi *(FWHM)^2, FWHM is the smoothness size of the map.


% Optional Parameters (OP): 
% 	First OP : Z_th: threshold for z-map
% 	Second OP: p_val : is the p_value of interest, the default is 0.05
% 	Third OP: tail: you can 'one' (default) for one sided comparison, or 
% ‘two’ for two sid


% Output Parameters: 
%   C_th_map: is cluster thresholded map that has same dimentions as z_map
 
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%        Author: Mahlega Hassanpour
%        Date : Sat September 08 15:47:14 EDT 2018
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% set defaults
nPx=344;
nPy =400;

C_th_map=zeros(size(z_map));
%% Intensity Thresholding
if nargin > 2 
    Z_th = varargin{1};
    z_map(z_map<Z_th)=0;
else
    Z_th = min(z_map);
end

%% Cluster Thresholding
foo=reshape(z_map,nPx,nPy);
bw_data = foo; % Binary map
bw_data(foo<0) = 0;
bw_data(foo>0) = 1;
cc_Th = bwconncomp(bw_data , 8); % find clusters
numPixels = cellfun(@numel,cc_Th.PixelIdxList);
listPixels = cc_Th.PixelIdxList;
xx2=find(numPixels>C_S_th);
for cc=1:length(xx2)
    C_th_map(listPixels{xx2(cc)})=  z_map(listPixels{xx2(cc)});
end

%% Saving
if nargin > 3
    file_name=varargin{2};
    
    save([pwd,'/',file_name, '_cluster_thresholded'], 'C_th_map', 'Z_th','C_S_th' )
    
 end
end
