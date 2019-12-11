function [gradG,resG] = global_BCs(gradG,resG,BC)
% Enforce dirichlet boundary conditions according to the nodes specified
nBC = size(BC,1);
% Boundary conditions
for iBC=1:nBC
    % If this BC is a dirichlet condition
    if (BC{iBC,1}(1)=='d')
        % Check the direction/type to get index modifier
        switch (BC{iBC,1}(2:3))
            case('xf')
                im = -4;
            case('yf')
                im = -3;
            case('pf')
                im = -2;
            case('xs')
                im = -4;
            case('ys')
                im = -3;
            case('xm')
                im = -1;
            case('ym')
                im = 0;
        end                
        % Loop through nodes in BC
        for nn = 1:length(BC{iBC,3})
            ind = 5*BC{iBC,3}(nn)+im;
            resG(ind) = 0.0;
            gradG(ind,:) = 0.0;
            gradG(ind,ind) = 1.0;
        end
    end
end