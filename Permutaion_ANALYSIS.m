function [z_map,z_map_th,all_perms,p_val_map,z_map_th_neg,p_val_map_neg,Keep_Pix]= Permutaion_ANALYSIS(S_map1,S_map2,N_Perm,file_name,varargin)
% Permutation analysis for group comparison



% Purpose: To perform permutation analysis for group comparison

% Input Parameters: 
%       S_map1 and S_map2: Matrices with dimensions of N_s1 by N_Vox and N_s2 by
% N_Pix, respectively. N_s1 and N_s2 are number of subjects in group one 
% and  two, respectively. N_Pix is the number of pixels in the somatomap image.
%       N_Perm : Number of permutations
%       file_name: Name of the output file

% Optional Parameters (OP): 
% 	First OP : all_perms: all the possible permutations given N_s1, N_s2 
% and N_Perm
% 	Second OP: p_val : is the p_value of interest, the default is 0.05
% 	Third OP: tail: you can 'one' (default) for one sided comparison, or 
% ‘two’ for two sided comparison, default is ‘one’


% Output Parameters: 
%       z_map: non-thresholded z-map
%       z_map_th: thresholded z- map for group 1> group 2
%       all_perms: see above
%       p_val_map : p-value map for group1> group2 comparison
%       z_map_th_neg: thresholded z- map for group 2> group 1
%       Keep_Pix : all the pixels that were hit by at least one participant

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%        Author: Mahlega Hassanpour
%        Date : Sat September 08 15:47:14 EDT 2018
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% set defaults
nPx=344; 
nPy =400;
FWHM=6; % 6 pixel smoothing

% Condition data
S_map1=S_map1./max(S_map1(:));
S_map2=S_map2./max(S_map2(:));

N_s1=size(S_map1,1);
N_s2=size(S_map2,1);

%% Proportional maps
Prop1=sum(S_map1)./N_s1;
Prop2=sum(S_map2)./N_s2;

%% Spatial smoothing

HSize = round(1.5*FWHM); % smoothing extent
Sigma= FWHM/sqrt(8*log(2));
H = fspecial('gaussian',HSize,Sigma);

Prop1 = Smooth_SMaps(Prop1,H,[nPx,nPy],[]);
Prop2 = Smooth_SMaps(Prop2,H,[nPx,nPy],[]);

%% Calculate z statistics for difference between two group
% Contrast map
Prop= (N_s1.*Prop1+N_s2.*Prop2)./(N_s1+N_s2);
diff_Prop=Prop1-Prop2;
Var_Prop=Prop.*(1-Prop).*(1/N_s1+1/N_s2);

Keep_Pix=find(Var_Prop>0);
% Z_map
z_map=zeros(1,size(S_map1,2));
z_map(Keep_Pix)=diff_Prop(Keep_Pix) ./ sqrt(Var_Prop(Keep_Pix));

clear Prop Prop1 Prop2 diff_Prop Var_Prop

nargin
%% Set parameters and orders for permutation test
N_m=max(N_s1,N_s2);
if nargin > 4
    all_perms=varargin{1};
    if isempty(all_perms)
        all_perms=single(nchoosek(1:N_s1+N_s2,N_m));
        perm_orders=randsample(length(all_perms),N_Perm);
        all_perms=all_perms(perm_orders,:);
        save(['All_',num2str(N_m),'selections_outof_',num2str(N_s1+N_s2)],'all_perms')
    end
else
    all_perms=single(nchoosek(1:N_s1+N_s2,N_m));
    perm_orders=randsample(length(all_perms),N_Perm);
    all_perms=all_perms(perm_orders,:);
    save(['All_',num2str(N_m),'selections_outof_',num2str(N_s1+N_s2)],'all_perms')
    
end

if nargin > 5
    p_val=varargin{2};
else
    p_val=0.05;
    
end

if nargin >6
    tail=varargin{3};
else
    tail='one';
    
end
%% Start Permutation
S_map_all=[S_map1;S_map2]; % combine all data
foo=1:N_s1+N_s2;
z_map_random=zeros(length(Keep_Pix),1);
S_map_all=S_map_all(:,Keep_Pix);
parfor or=1:N_Perm
    if N_s1-N_s2 >0
        perms1{or}=all_perms(or,:);
        perms2{or}=foo(~ismember(foo,perms1{or}));
    else
        perms2{or}=all_perms(or,:);
        perms1{or}=foo(~ismember(foo,perms2{or}));
    end
        
    Prop1(:,or)=sum(squeeze(S_map_all(perms1{or},:)))./N_s1;
    Prop2(:,or)=sum(squeeze(S_map_all(perms2{or},:)))./N_s2;
    
     % Apply smoothing
     Prop1(:,or)= Smooth_SMaps_2(Prop1(:,or),H,[nPx,nPy],Keep_Pix);
     Prop2(:,or)= Smooth_SMaps_2(Prop2(:,or),H,[nPx,nPy],Keep_Pix);
    
    % Calculate z_map
        
    Prop(:,or)= (N_s1.*Prop1(:,or)+N_s2.*Prop2(:,or))./(N_s1+N_s2);
    
    diff_Prop(:,or)=Prop1(:,or)-Prop2(:,or);
    Var_Prop(:,or)=Prop(:,or).*(1-Prop(:,or)).*(1/N_s1+1/N_s2);
    
    z_map_random(:,or)=diff_Prop(:,or) ./ sqrt(Var_Prop(:,or));
    
    
end 
clear Prop Prop1 Prop2 diff_Prop Var_Prop
%% Distribution evaluation
B=unique(z_map(Keep_Pix));
p_value_b=zeros(size(B));
p_value_b_neg=p_value_b;

p_val_map=zeros(size(z_map));
p_val_map_neg=zeros(size(z_map));

z_map_th= z_map;
z_map_th_neg= -z_map;

z_rand=z_map_random(:);

switch tail
    case 'one'
        for bb=1:length(B)
           
            p_value_b(bb)=sum(z_rand>B(bb))/length(z_rand);
           
            p_value_b_neg(bb)=sum(-z_rand>-B(bb))/length(z_rand);
            
        end
        % find positive threshold
        B_sub= B( p_value_b<=p_val);
        Z_th_P=min(B_sub);
        if ~isempty(Z_th_P)
            z_map_th(z_map<Z_th_P)=0;
        end
        % find negative threshold
        B_sub= B( p_value_b_neg<=p_val);
        Z_th_N=-max(B_sub);
        if ~isempty(Z_th_N)
            z_map_th_neg(z_map_th_neg<Z_th_N)=0;
                  
        end
    case 'two'
        parfor bb=1:length(B)
            p_value_b(bb)=sum(abs(z_rand)>abs(B(bb)))/length(z_rand);
            
        end
        B_sub= B( p_value_b<=p_val);
        Z_th=min(B_sub);
        z_map_th=abs(z_map);
        z_map_th(z_map<Z_th)=0;
        
end
save([pwd,'/',file_name,'_',num2str(N_Perm),'_',tail,'_tail'],'z_map','z_map_th','all_perms','p_value_b','z_map_th_neg','p_value_b_neg','z_rand','Keep_Pix','Z_th_N','Z_th_P','-v7.3')
end
