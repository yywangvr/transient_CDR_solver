function nu_sc = shockcapturing(para,a,Projectedgrtem,grtem)
% This function is used to capture shock


OrthogonalGrTem=grtem-Projectedgrtem;

NormgrTem=norm(grtem);
NormOGrtem=norm(OrthogonalGrTem);
%S=eye(length(a))-a'*a/norm(a);

if NormgrTem >  1e-12
    nu_sc=0.5*para.alpha*(norm(a)*para.h+para.sigma*para.h^2)*NormOGrtem/NormgrTem;
else
    nu_sc=0;
end




end

