function B= Smooth_SMaps(A,H,dim,Keep_Pix)



%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%        Author: Mahlega Hassanpour
%        Date : Sat September 08 15:47:14 EDT 2018
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if ~isempty(Keep_Pix)
    
    A2=zeros(dim);A2(Keep_Pix)=A;
else
    A2=A;
end

if isempty(find(isnan(A2)))
    
    B=imfilter(reshape(A2,dim),H);
else
    B=A;
end

if ~isempty(Keep_Pix)
    
    B= B(Keep_Pix);
end
B= reshape(B,size(A));
'Done'
end


