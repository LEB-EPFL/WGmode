%------------------------------------------------------------------%
%             findRootBisection - Main function                    %
%------------------------------------------------------------------%

%  [File name]	"findRootBisection.m"
%  [Author]		"Anna Archetti"
%  [Version]		"1.0"
%  [Modified by]	"Anna Archetti"
%  [Date]			"27 March 2017"

% findRootBisection uses the bisection method to find the roots of a 
%                   function of discrate points.
% findRootBisection is a recursive function until a root is found within the 
%                   tolerance set by the function call.
%
% INPUTs:
% - y: vector with the numeric values of the funciton of x of which we have
%      to find roots
% - x: vector with the numeric values of independent variable (same length as y)
% - rootNum: number of roots
% - NMAX: maximum numbero fo iterations
% - errTol: maximum error tolerance on the yc

% OUTPUTs:
% - xc: x values of the roots vectors
% - yc: y values of the roots vectors


function [xRoot, yRoot] = findRootBisection(x, y, rootNum, errTol, NMAX)

    % Check inputs
    if (nargin == 2)
       rootNum = 1;
       NMAX = 10^6;
       errTol = 0.3;
    elseif (nargin == 3)
       NMAX = 10^6;
       errTol = 0.3;
     elseif (nargin == 4)
       NMAX = 10^6; 
    end
    
    % Set initial interval extremes of the first search domain [x1, x2]
    val = 0;
    idxX1 = 0;
    idxX2 = 0;
    if y(1) <= 0
        [val, idxX1] = min(y);
        [val, idxX2] = max(y);
        
    elseif y(1) > 0
        [val, idxX1] = max(y);
        [val, idxX2] = min(y);
    end
    
    if idxX1 > idxX2
        idxTemp = idxX1;
        idxX1 = idxX2;
        idxX2 = idxTemp;
    elseif idxX1 == idxX2
        idxX2 = idxX1 + 1;
    end
    x1 = x(idxX1);
    x2 = x(idxX2);
    c = 0;
    itNum = 0;
    err = max(y);

    % For each root
    for rootIdx = 1: rootNum
        
        idxCPos = 1;
        % while the number of iterations is minor than NMAX
        while itNum <= NMAX;   
            % interval midpoint: root guess
            c = (x1 + x2)/2;
            [val, idxCPos] = min(abs(x(idxX1:idxX2) - c)); 
            idxCPos = idxX1 + idxCPos - 1;
            err = abs(y(idxCPos));
            if err <= errTol 
                itNum = 0;
                break;
            elseif  y(idxCPos) > 0
                x1 = x(idxCPos);
                idxX1 = idxCPos;
            else  
                x2 = x(idxCPos);
                idxX2 = idxCPos;
            end    
            itNum = itNum + 1;
            
        end
        
        if itNum >= NMAX
            error('The function reached the maximum number of iterations')
        end
        
        % Save root
        xRoot(rootIdx) = x(idxCPos);
        yRoot(rootIdx) = y(idxCPos);
        
        % Update initial interval extremes of the new serach domain
        if y(idxCPos + 1) <= 0
            [val, idxX1] = min(y(idxCPos + 1 : end));
            [val, idxX2] = max(y(idxCPos + 2 : end));
        elseif y(idxCPos + 1) > 0
            [val, idxX1] = max(y(idxCPos + 1 : end));
            [val, idxX2] = min(y(idxCPos + 2 : end));
        end
        err = max(y);
        
        idxX1 = idxX1 + idxCPos;
        idxX2 = idxX2 + idxCPos + 1;
        if idxX1 > idxX2
            idxTemp = idxX1;
            idxX1 = idxX2;
            idxX2 = idxTemp;
        elseif idxX1 == idxX2
            idxX2 = idxX1 + 1;
        end
        
        x1 = x(idxX1);
        x2 = x(idxX2);      
    end
    
    figure,
    plot(x, y, 'k*'),
    for rootIdx = 1: rootNum
        hold on
        plot(repmat(xRoot(rootIdx), 1, length(y)), y, 'r--'),
        plot(x, repmat(yRoot(rootIdx), 1, length(y)), 'r--'),
        disp(['The root ' num2str(rootIdx) ' of the function is: ' num2str(xRoot(rootIdx)) ])
    end
    xlabel('x')
    ylabel('y')
    legend('function values', 'roots')
end

