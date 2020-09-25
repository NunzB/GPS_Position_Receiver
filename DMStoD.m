
% DMS = [DEGREE MINUTE SECOND]

% D = [DEGREE]



function [D] = DMStoD (DMS)

    D = zeros(size(DMS,1),1);

    for i=1:size(DMS,1)
        
        for j=1:size(DMS,2)

            D(i) = D(i) + DMS(i,j) / 60^(j-1);
    
        end
    end

end