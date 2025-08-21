function [X,solinfo,sol] = ...
    solveAndContinue(x0,fun_residual,lam_s,lam_e,sOpt,solOpt)

%[text] ## dir helps in fixing stepping direction, towards higher frequency or towards lower frequency
% Determine continuation direction
if lam_s>lam_e
    dir             = -1;
else
    dir             = 1;
end
%[text] Initialize result vectors
X                   = zeros(length(x0)+1,sOpt.stepmax); % Matrix whose columns are the solution points to R(x)
sol                 = struct;
sol(1:sOpt.stepmax) = struct;
solinfo.NIT         = zeros(sOpt.stepmax,1);            % Number of iterations in each step
solinfo.IEx         = zeros(sOpt.stepmax,1);            % Exit flags of the nonlinear solver in each step
solinfo.FC          = zeros(sOpt.stepmax,1);            % Number of function evaluations in each step
solinfo.ctime       = 0;                                % Total computation time
solinfo.FCtotal     = 0;                                % Total number of function evaluations
%[text] Initial guess,  reference vector, tangent vector and additional row for extended jacobian
X0                  = [x0; lam_s];                      % Combine x0 and lam_s into single augmented solution vector
                                                        % R(y,omega) =0, is what this function is solving. 
zref                = [zeros(length(x0),1); 1];         % Initializes reference tnagent vector (assumes the solution path 
                                                        % primarily changes along the \lambda-direction.
Xref                = X0;                               % Sets the reference solution for pseudo arclength constraint to initial solution X0,
Xold                = X0;                               % saves the current solution X0 as Xold before X0 is updated with the newly computed solution
%[text] Relax solver option for first iteration
solOptTemp        = solOpt;                             % Save courrent solver options in a temperory variable for reseting the options after first iteration
solOpt            = optimset(solOptTemp,'MaxIter',500,'Display','iter');    % Change only MaxIter to 500, and display to iter. Other values remain same                                               
%[text] Solve for first point with relaxed options using fsolve ( solve system of nonlinear equations)
%[text] \[Xp, fval, iEx, output, J\] = fsolve(@(X) extended\_residual(X, Xref, zref, fun\_residual,sOpt), X0, solOpt); 
%[text] Xp     :  The solution vector where fun(x) =0.
%[text] fval    :  The value of function at the solution.
%[text] iEx    :  Exit flag, 1 =\> solution found.  0 =\> Too many iterations or function evaluations. -2  =\> No solution found.
[Xp, ~, ~, ~, J] = fsolve(@(X) extended_residual(X, Xref, zref, fun_residual,sOpt), X0, solOpt); 
X0               = Xp;        
%[text] Save first point of solution. 
X(:,1)                   = X0;  % X is the output vector, first point obtained from fsolve is added to first column of X.
%[text] Dispaly continuation on command window ( For inspection)
disp(['Continuation at ' num2str(X0(end),'%.4f') ...
    ', step size ' num2str(sOpt.ds) '.']);  
%[text] Reset solver option back.  ( Relaxation in number of iteration is given only for first solution point)
solOpt                   = solOptTemp;
%[text] Starting continuation
istep                    = 2;
while istep<=sOpt.stepmax
    % Determine the predictor direction
    % max(abs(Xref), 1e-4) compares each element of abs(Xref) with 1e-4 and returns the element-wise maximum.
    % If an element of abs(Xref) is greater than 1e-4, it stays unchanged , If an element is smaller than 1e-4, it is replaced with 1e-4.
    % This avoids division by values that are too small (helps with numerical stability)
    [~,kk] = sort(abs(zref./max(abs(Xref),1e-4)),...
        1,'descend');
    % Temporarily switch off warning
    warning('off','MATLAB:nearlySingularMatrix');
    warning('off','MATLAB:singularMatrix');
    for ik=1:length(kk)
        k = kk(ik);
        c = zeros(length(X0),1); 
        c(k) = 1;
        % Determine unit tangent to the solution
        % path (Eq. 4.8)
        ztmp = [J(1:end-1,:);c']\...
            [zeros(size(J,1)-1,1);1];
        if ~any(isnan(ztmp))
            % Successful!
            break;
        end
    end
    % Switch warning on again
    warning('on','MATLAB:nearlySingularMatrix');
    warning('on','MATLAB:singularMatrix');

    if any(isnan(ztmp))
        error('Could not determine predictor direction.');
    end
    % Apply linear scaling to tangent
    z = ztmp/norm(ztmp);
    % Take step
    XP = X0+dir*sOpt.ds*z;

    % Ensure forward stepping along solution path
    if (istep > 2) && ...
            (transpose(X0-Xold)*(dir*sOpt.ds*z) < 0)
        XP = X0-dir*sOpt.ds*z;
    end

    % corrector
    % For orthogonal or normal parametrization, the predicted 
    % solution is relevant
    Xref = XP;
    zref = z;
    % Solve extended nonlinear system of equations
    [Xtmp,~,iEx,output,Jtmp] = ...
            fsolve(@(X) extended_residual( ...
            X,Xref,zref,fun_residual,sOpt),XP,solOpt);
    % POSTPROCESS SOLUTION POINT
    % Save previous solution
    Xold = X0;
    % Update solution vector
    X0 = Xtmp;
    J = Jtmp;

    % Save solution point
    X(:,istep)         = X0;
    solinfo.IEx(istep) = iEx;
    solinfo.NIT(istep) = output.iterations;
    solinfo.FC(istep)  = output.funcCount;
    
    if istep==sOpt.stepmax
        disp(['Terminating continuation since maximum number of ' ...
            'solution points ' num2str(sOpt.stepmax) ' is reached.']);
    end
     % Increment loop count
    istep = istep+1;
end


%[appendix]{"version":"1.0"}
%---
