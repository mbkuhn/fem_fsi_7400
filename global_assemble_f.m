function [gradG,resG] = global_assemble_f(gradG,resG,gradGlf,gradGlm,...
    resGlf,resGlm,ne,econn,Jr,Jf,BC)
% Assembles the global tangent matrix and residual for fluid elements
% based on the element information supplied.

% Boundary conditions
nBC = size(BC,1);
mflag = 0;
for iBC=1:nBC
    %  Only care about the matching BC here
    if (BC{iBC,1}(1)=='m')
        % Check element list for BC
        if (fastintersect(ne,BC{iBC,2}))
            mflag = 1;
        end
    end
end

npe = length(econn);
for nna=1:npe
    A = econn(nna);
    
    Rindf = 5*A-4:5*A-2; rindf = 3*nna-2:3*nna;
    Rindm = 5*A-1:5*A  ; rindm = 2*nna-1:2*nna;
    
    % If matching BC exists for this element
    if (mflag)
        % Check node list for BC
        if (fastintersect(A,BC{iBC,3}))
            % Do not allow mesh eq's to contribute at matching BC
            Rindm = []; rindm = [];
        end
    end
    
    % Global residual vector
    resG(Rindf) = resG(Rindf) + resGlf(rindf)*Jf;
    resG(Rindm) = resG(Rindm) + resGlm(rindm)*Jr;
    for nnb=1:npe
        B = econn(nnb);
        
        Cindf = 5*B-4:5*B-2; cindf = 3*nnb-2:3*nnb;
        Cindm = 5*B-1:5*B  ; cindm = 2*nnb-1:2*nnb;
        % Global tangent matrix
        gradG(Rindf,Cindf) = gradG(Rindf,Cindf) ...
            + gradGlf(rindf,cindf)*Jf;
        gradG(Rindm,Cindm) = gradG(Rindm,Cindm) ...
            + gradGlm(rindm,cindm)*Jr;
    end
    
end