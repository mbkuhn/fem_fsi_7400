function [gradG,resG] = global_assemble_s(gradG,resG,gradGl,resGl,...
    ne,econn,J,BC)
% Assembles the global tangent matrix and residual for solid elements
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
    
    Rind = 5*A-4:5*A-3; rind = 2*nna-1:2*nna;
    Pind = 5*A-2; % pressure index
    
    % Set pressure to 0 in internal solid nodes
    pflag = 1; 
    % If matching BC exists for this element
    if (mflag)
        % Check node list for BC
        if (fastintersect(A,BC{iBC,3}))
            % Do not touch pressure eq for a matching node
            pflag = 0;
        end
    end
    
    % Modify pressure eq if needed
    if (pflag)
        resG(Pind) = 0;
        gradG(Pind,:) = 0;
        gradG(Pind,Pind) = 1;
    end
    
    % Global residual vector
    resG(Rind) = resG(Rind) + resGl(rind)*J;
    % Global tangent matrix
    for nnb=1:npe
        B = econn(nnb);
        Cind = 5*B-4:5*B-3; cind = 2*nnb-1:2*nnb;
        gradG(Rind,Cind) = gradG(Rind,Cind) ...
            + gradGl(rind,cind)*J;
    end
    
    % Mesh motion indices and residual vector
    MindA = 5*A-1:5*A;
    resG(MindA) = 0.0;
    % Set mesh motion equal to solid/fluid solution at every solid node
    for nn=1:2
        gradG(MindA(nn),:) = 0.0;
        gradG(MindA(nn),MindA(nn)) =  1.0;
        gradG(MindA(nn),MindA(nn)-3) = -1.0;
    end
    
end