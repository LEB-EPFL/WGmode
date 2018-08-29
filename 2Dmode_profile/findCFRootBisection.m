%------------------------------------------------------------------%
%             findCFRootBisect - Main function                    %
%------------------------------------------------------------------%

%  [File name]	"findCFRootBisection.m"
%  [Author]		"Anna Archetti"
%  [Version]		"1.0"
%  [Modified by]	"Anna Archetti"
%  [Date]			"27 March 2017"

% findRootBisection uses the bisection method to find the roots of a 
%                   CONTINUOUS function.
% findRootBisection is a recursive function until a root is found within the 
%                   tolerance set by the function call.
%
% INPUTs:
% - myFunc: the anonymus function 
% - domInt: a vecotr of two number, upper and lower limit [L, U]
% - rootNum: number of roots
% - NMAX: maximum numbero fo iterations
% - errTol: maximum error tolerance on the yc

% OUTPUTs:
% - xc: x values of the roots vectors
% - yc: y values of the roots vectors


function [xRoot, yRoot] = findCFRootBisection(myFunc, domInt, errTol, NMAX)

    % Check inputs
    if (nargin == 2)
       NMAX = 10^4;
       errTol = 0.1;
    elseif (nargin == 4)
       NMAX = 10^4;
    end
    
    % Check initial interval extremes of the first search domain [x1, x2]    
    if domInt(1) > domInt(2)
        error('The lower bound limit is higher than the upper one')
    elseif domInt(1) == domInt(2)
        error('The lower bound limit is equal to the upper one')        
    end
    
    x1 = domInt(1);
    x2 = domInt(2);
    c = 0;
    itNum = 0;
    err = 1000;

        
    % while the number of iterations is minor than NMAX
    while itNum <= NMAX;   
        % interval midpoint: root guess
        c = (x1 + x2)/2; 
        err = abs(myFunc(c));
        if err <= errTol 
            itNum = 0;
            break;
        elseif  myFunc(x1)*myFunc(c) >= 0 
            x1 = c;
        else  
            x2 = c;
        end    
        itNum = itNum + 1;

    end

    if itNum >= NMAX
        error('The function reached the maximum number of iterations')
    end
    
    clear itNum;
    % Save root
    xRoot = c;
    yRoot = myFunc(c);
    
    x = (domInt(1):0.001:domInt(2));
    y = myFunc(x);
    plot(x, myFunc(x), 'k'),

        hold on
        plot(repmat(xRoot, 1, length(y)), y, 'r--'),
        plot(x, repmat(yRoot, 1, length(y)), 'r--'),
     %   disp(['The root of the function is: ' num2str(xRoot) ])
        
    xlabel('x')
    ylabel('y')
    legend('function values', 'roots')
end

