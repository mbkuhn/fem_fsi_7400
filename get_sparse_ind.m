function [gradG,resG] = get_spars_ind(gradG,resG,gradGlf,gradGlm,...
    resGlf,resGlm,ne,econn,Jr,Jf,BC)
% Determines the indices needed for the sparse treatment of gradG

indflagvec = zeros((5*nnode+1)*(5*nnode),1);

nentry = sum(indflagvec);
trow = zeros(nentry,1); tcol = trow; gradGvec = trow;
for ie = 1:nentry
    trow(ie) = irow;
    tcol(ie) = icol;
end

for iBC=1:nBC
    %  Only care about the matching BC here
    if (BC{iBC,1}(1)=='m')
        mBC = [mBC iBC];
    end
end

nacc = 0;

for ne=1:nelem
    
    mflag = 0;
    for iBC=1:length(mBC)
        %  Only care about the matching BC here
        if (fastintersect(ne,BC{iBC,2}))
            mflag = 1;
        end
    end

    if (solver.mesh.tag(ne)) % solid
        
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
                nacc = nacc+1;
                trow(nacc) = Pind;
                tcol(nacc) = Pind;
            end
            
            % Global tangent matrix
            for nnb=1:npe
                for ii = 1:2
                    for jj = 1:2
                        nacc = nacc+1;
                        trow(nacc) = Rind(ii);
                        tcol(nacc) = Cind(jj);
                    end
                end
            end
            
            % Mesh motion indices and residual vector
            MindA = 5*A-1:5*A;
            % Set mesh motion equal to solid/fluid solution at every solid node
            for nn=1:2
                nacc = nacc+1;
                trow(nacc) = MindA(ii);
                tcol(nacc) = MindA(ii);
                trow(nacc) = MindA(ii);
                tcol(nacc) = MindA(ii)-3;
            end
            
        end
    else % fluid
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
                for nn=1:3
                    
            end
            
        end
end