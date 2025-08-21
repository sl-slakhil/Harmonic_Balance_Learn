%[text] ## Two degree of freedom system Harmonic balance base code
%[text] The code is heavily inspired from NL-Vib. Functions are edited and copied from NL-Vib, the modification removed many options. This code provide only single options. 
%[text] There is only on predictor method and one corrector method.
%[text] The example used is also from the Book : Harmonic balance mehod for nonlinear vibrations.
%[text] ## Aim of this code is to assist in learning coding Harmonic balance method.
%[text] 
%   The example is from:
%   M. Krack, J. Gross: Harmonic Balance for Nonlinear Vibration
%   Problems. Springer, 2019. https://doi.org/10.1007/978-3-030-14023-6.

close all
clear variables
clc
%[text] M is the mass matrix, K is the stiffness matrix, D is the damping and Fex1 is the external force amplitude.
system.M          = [1   0; 0  0.05];
system.K          = [1.0453 -0.0453; -0.0453 0.0453 ];
system.D          = [0.0150 -0.0130; -0.0130 0.0130];
system.Fex1       = [0.11;  0];
%[text] Nonlinear spring in the two degrees of freedom is given as:
% Unilateral spring applied to 1st DOF 
nonlinear_elements{1} = struct('type','cubicSpring',...
    'stiffness',1,'force_direction',[1;0]);
% Unilateral spring applied to 2nd DOF
nonlinear_elements{2} = struct('type','cubicSpring',... 
    'stiffness',.0042,'force_direction',[1;-1]);
% Adding nonlinear element to system
system.nonlinear_elements = nonlinear_elements;
%[text] Inital conditions and solutions
analysis  = 'FRF';

% Analysis Parameters
H              = 7;
N              = 4*H+1;
Om_s           = 0.7;
Om_e           = 1.4;

% Initial guess ( Solution of Underlying linear system)
Q1             = (-Om_s^2*system.M+ 1i*Om_s*system.D + system.K)\system.Fex1;
y0 = zeros((2*H+1)*length(Q1),1);
%[text] Why ? In the below line of code, for 2 degree of freedom and for 7 harmonics. The length of y0 is 30. 2(2\*7+1) =30. The vector y0 is the stacked form of all these coefficient. First two are dc components. Third and fourth is real coefficient of harmonics one. ( for both degrees of freedom ) fifth and sixth are the imaginary coefficient of harmonic 1. (( for both degrees of freedom ) ..... Similarly 29 and 30 th are the imaginary to see maping run the below line in command window: \[1:30; zeros(1,2) ones(1,4) 2\*ones(1,4) 3\*ones(1,4) 4\*ones(1,4) 5\*ones(1,4) 6\*ones(1,4) 7\*ones(1,4)\]
y0(length(Q1)+(1:2*length(Q1))) = [real(Q1);-imag(Q1)];
ds          = 0.1;
%[text] 
% sOpt.eps : Tolerence of the norm of the residual in order ro identify
% failure in convergence.
sOpt.eps            = 1e-4;
% Maximum number of steps before termination
sOpt.stepmax        = 1e2;
% ds : Normalized step size
sOpt.ds             = ds;
%  sOpt is the solution options dscale is the scaling value, initially set
%  for one(no scaling)
sOpt.Dscale         = ones(length(y0)+1,1);
sOpt.flag           = 1;
%[text] - Each keyword modifies how the solvers run
%[text] - - Display - 'off '  no output shown (silent mode)
%[text] - -' iter '   display results at each iteration
%[text] - - 'final ' display only final result
%[text] - - ' off '  supresses solver messages. \
solOpt              = optimset('Display','off',...%'iter',...
                    'Jacobian','on','MaxIter',50,...
	                'TolFun',sOpt.eps,'TolX',sOpt.eps);%,'DerivativeCheck','on');
%[text] # Solve and continue using predictor corrector
[X_HB, Solinfo_HB]        = solveAndContinue(y0, @(X) HB_residual(X,system,H,N),Om_s,Om_e,sOpt,solOpt); %[output:1aa40d04]
%[text] Interpret solver output
Om_HB                     = X_HB(end,:);
Q_HB                      = X_HB(1:end-1,:);
%[text] Define amplitude as magnitude of the fundamental harmonic of the second 
%[text] coordinate's displacement
n=2;
a_HB                    = sqrt(Q_HB(n+2,:).^2 + Q_HB(2*n+2,:).^2);
a_rms_HB                = sqrt(sum(Q_HB(1:2:end,:).^2))/sqrt(2);
%[text] Plot results
figure; hold on; %[output:0928a474]
plot(Om_HB,a_rms_HB,'g-'); %[output:0928a474]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":26.7}
%---
%[output:1aa40d04]
%   data: {"dataType":"text","outputData":{"text":"\n                                             Norm of      First-order   Trust-region\n Iteration  Func-count     ||f(x)||^2           step       optimality         radius\n     0          1         0.000113672                          0.0106              1\n     1          2         2.49651e-08       0.040723         0.000134              1\n     2          3         1.33102e-15    0.000638394         3.23e-08              1\n\n<a href = \"matlab: helpview('optim','eqn_solved','CSHelpWindow');\">Equation solved<\/a>.\n\nfsolve completed because the vector of function values is near zero\nas measured by the value of the <a href = \"matlab: helpview('optim','fcn_tolerance_fsolve','CSHelpWindow');\">function tolerance<\/a>, and\nthe <a href = \"matlab: helpview('optim','appears_regular','CSHelpWindow');\">problem appears regular<\/a> as measured by the gradient.\n\n<<a href = \"matlab: createExitMsg({'optim:fsolve:Exit1basic','fsolve'},{'optim:fsolve:Exit1detailed','3.227515e-08','1.000000e-04','1.331023e-15','1.000000e-02'},true,true);;\">stopping criteria details<\/a>>\nContinuation at 0.7000, step size 0.1.\nTerminating continuation since maximum number of solution points 100 is reached.\n","truncated":false}}
%---
%[output:0928a474]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAPIAAACSCAYAAAB\/nJf+AAAAAXNSR0IArs4c6QAAD3dJREFUeF7tXU1IHckWPpPEJDKJDJmFKHlDBMH1SxZZBIJrFw66kIkb3bgIjjuJIHIDcglkcGcEQYW4cUAC4d2Fa5GXhQtv1hcEYQauJOENQ8ZBzGTyHqfzSvre3J+u7vrvryGLeLur6nznfH1OnTpd9dWNGzf+S7iAABDwGoGvQGSv9YfBA4EIARAZhgAEAkAARA5AiRABCIDIsAEgEAACIHIASoQIQABEhg0AgQAQAJEDUCJEAAIgMmwACASAAIgcgBIhAhAAkWED2hH49N0n+uveX\/Tx3ke68MuFqL\/Onzq195unDqwRuaenh4aGhmhnZ4eOj4\/zhHluZGUCn\/1wRqdzpxGBL\/z6mcRM6EuvLlHXcFdusNAtqDUi3759m1ZWVmh0dBRE1q1lC+2fPjo9J\/CVn6\/UeGAm+O+vfweZFeolFZHv3r1LxWKRrl27RtVqlWZnZ+no6KjpsJ4\/fx79Njk5eX4PiKxQiw41JcLoP5\/9SZ1PO5uG0Bxq\/1H6g64PX6eOVx0OSSA\/FJZZTBnkn1bzRCoiMzHfvHlDq6urtLS0RLu7u7S8vNxwRDMzMzQ+Pk6VSgVEVqMzZ1thg37\/r\/fR+Oq9cKNBvy+9p45\/d3g9X2YZPv3jUyJ5dSpOmsjCG5dKpYi8T58+pe7u7hqSigHzvY8ePaJ3797R1atXQWSdmrTcdpzEXd93JfJQZw\/OiD33jW9vWB69fPdiesCe+Osfv7YeXaQicqFQoM3NTdre3o6I3N\/f3zC8Zs99cHBAN2\/e\/ILsCK3ljcfVJ+JGnZTEQpbf\/vNbRAT24L5ccXm\/+ec30bB5zp8kCtElozYic0h9586dyAs38tqCyOvr67SxsaFLPrSrGYFGRi3TJROA58hMZh+uZvJyiH3xl4vW5EhFZE50tQut2RsPDAzU6CY+TxZEnp6epnK57IMOMcY6BIRRsydKS0QOrflK+7xJpcSnD8ITi\/69IzIPXCbZxfe38shYfjJpiur6UkFiHo0v82SW9+TZSZTYajR98C60ZvDjy09xL8uE5Wtubq7GYkBkdQRyoaWs4XRcBkFk9nC2l3CaYcvy8ro4L5k1ywF4SWQVxoRklwoUzbfRKrxMOxrXE15MYn7hcPjfbM3btgzSc+S0yqp\/DkRWhaS5dtIsMSUZnW0StBqjiBhaFbeICMVm9h1ETmJpuIfazRGzQORq5jppHkDcZ7NKDUTOYoE5eVYniRlCFzPXMlMIUW5qc54PIueEjFnETDJHVNF+\/ZJOljazPCv74hJEtlmhBiJn0XgOnhVfMbWaI2aFwbUlKNkXl7jf5osIRM5qhQE\/nyTRo0J8F0JTIYd4cckkrkBkfI+sggda2kia6FHVuQuZ67QvLp7j\/\/3d31Y3SoBHVmWJAbUjk+hRJTZnrnn7H1sfT2R5cdkuz2QdgMiqLDGQdmQTParEtrkElfXFZbuqC0RWZYUBtSOb6FEluq3wVMWLi6cFOpOBSTCGR06CUk7uESS2EeLayFwnqaFup3oXikHgkdtpKUe\/p030qILIxscTKqIPEBm7aKriQOZ2siR6Mncea8Bk5lrV+rgrS2cIrVVaoodtZU30qBTZVMJLFYlZdhequhBaq7RCT9viJFOr72xNisVj4RBbZ6mj6imEC8UgILJJK3WwLxVzRJVi6Z4nqyYxyw4iY46skgPSbekwaulBNHhA1zxZl7y2ls3qocMcWYX1edaGK8mtRrCxh+PxqdyMT+WcuH7MXNXFl+1zrEBkz0iYdbjxAgibX+s0k0N1eK2TxCyDC+WZmCNnZYWHz7uU3GoGn4rsNb+wxEmQOquuQGTMkY2\/BtJ8omd8kLEtctPuuBGPOnSf\/uBCnTU8sg0rtdSnWO\/U6Z1UisYE4fOUZeeeLKfYOkj2+Jo04weR4ZHT2E2qZ0TRRxpipOpQwUOCkEk9ajyUNnmIuq4suyyESHbJIubZ\/Sq+7rElcpIlI0FgvpevpMRXIZML2+AKOUBkFRp1uA3Xij5koWKCsgyCpBxViIt\/4yNcTBNY9O\/KBxOYI8talWf3J\/FoPogUD5vFePl4GSY1e2Dbu4rY3M8aHtkHC84wRpeLPjKI5dSj8MhEhCNj9NmkS1806ZPSfsveEzl+GmO1WqXZ2Vk6OjqqQbavr4+Wlpaot7c3+vvW1hYtLy+f3wMi6zNEH4o+9ElvrmXviZzkfOT4UaozMzM0PDxMCwsLtL+\/HyENIusxODEvltmXWc9Iwm\/VayILb1wqlSIP2+js43oVMpEHBwdrPDeIrN7QhWH5UvShHgGzLbqyO0iqrDUTuVAo0ObmJm1vb0dE7u\/vbxteI7TWa2Q+Fn3oRUR\/6+LFmbaUVOUIpdeRZYgsBlrvxRFaq1QhaT3yVO1Iw2rNe49cLBZJJrQWia\/Dw0Oam5urmSOvr6\/TxsZGWBo2LI3vRR+G4VLWnddzZEYhabKL72Xijo2N0dTUFK2trUXheNwjT09PU7lcVgZu3hoKpejDR715T+T48lOlUqHJyclIDzxfFuTF8pN+00TRh36MW\/XgPZFVwIesdTYUXd\/pI5t0fjwNImMdObOlYl6cGUIlDeAzRnyPnNqQUPSRGjrlD4LIIHIqo8K8OBVs2h4CkUFkaePCxxDSkGl\/AFv9gMhSRubzTh9Sgnp2s4odP1WILF3ZpaLT+Dry6OgoHR8fq2o22HaQ3HJTtWKjP5Ub6qeRFEROg5rhZ1D0YRhwie5wZAxC60TmguRWIpis3YRD3EDktsaH5FZbiKzfIKIlnUfBJhESoXUSlCzdg50+LAEv0a3qs6okuq65FUROi5zm55Dc0gywouZd+ZQRRFakUJXNILmlEk39bfESVOdPnda25WUJQWT9epbqAcktKbicuNmFohAQ2QlT+DwIJLccUobEUFxYSwaRJRSm81ZUbulEV2\/b4rham5lrEFmvjhO3juRWYqicu9GFzDWI7IBZiDc6trF1QBkphuDCBgMgcgrFqXwEGWqVaNpry\/bnjCCyPd1HyS2R8bRddG8RhiC6tv0VFIhsyYyQobYEvKZubX88ASJrUmyrZuMk7vq+i\/isX1x+I2C75hpENmw\/WGYyDLjB7mzOk0Fkg4rmrrDMZBhwg93ZnCeDyAYVLUhsuy7XoMi56spmYQiIbMjUsFZsCGiL3Ygvoa4PX6eOVx1GRwIiG4AbJDYAsiNd2AqvQWTNBoCCD80AO9Y8L0Oxzk3XXYPIGg0BJNYIrqNN26q7BpE1GQQrlENqniuhaksTyI42a2MZKhWR48eqVqtVmp2dpaOjoxpY4\/fwD3t7e+eHnPP\/Qz6NUbyVr\/x8BSR2lGw6h2Xj++RURE5y0Lm4Rxx0zgeav3jxgpaXlyMMQyUywmmdFPGjbWEDJrPX0kQWnrZUKkWk5MPNu7u7zw87bwS1OPR8d3c3aCIjO+0H0UyM0nT2OhWRC4UCbW5u0vb2dkTk\/v7+huG1AGxsbIwmJiZocXGR9vf3g\/TIILEJevjTh+mkl3YiswePE1+oQoTW6+vrtLGx4Y+GGoxUVGzxnJirtnABAUbApFdOReRisUhJQutGnrieyDx3LpfLXmqeP4A4++GMTudOCbt7eKlCrYM26ZWlicySJ0l2MYlHRkZofn7+i4x2CMkuQWJWFjyxVj543TgvRZl4yacicnxpqVKpnCe6eL7M1+rqKi0tLVFvb2+NEra2toJIdsW\/J8YHEF7zTPvgTVV6pSKyCul9XX7CpgAqtJ+vNkwUiIDIEjYlvm7hHT2ws4cEcDm\/1YRXBpETGFk8qXXp1SXqGu5K8BRuAQKfERCbLHKpLudTdFwgchtUWQm8vMTeGEktHSaYjzbZfi7+elHb\/mwgcgs7wnw4HyQLQUoQuYkWRaUWh9LXfrym7U0aghFBBvsIgMh1OsD6sH2jxAjkEQCRY5gJL8xZaawPyxsTnrCHAIj8\/6ziybMT+njvo5EqHHvqRs+hIpBrIseXldgL8\/KA6d0PQzUsyGUWgdwSOZ6RxrKSWaNDb+oRyB2R670wKrTUGxVaNI9A7ojM34jyBS9s3tjQoz4Eckdk3RU2+lSFloFAcwRyR2QYAxAIEQEQOUStQqbcIQAi507lEDhEBEDkELUKmXKHAIicO5VD4BARAJFD1Cpkyh0CIHLuVA6BQ0QARA5Rq5ApdwiAyLlTOQQOEQEQOUStQqbcIQAi507lEDhEBEDkELUKmXKHAIicO5VD4BARAJFD1Cpkyh0CIHLuVA6BQ0QgFZHjpzFWq1WanZ1teHQqA9buoPPR0VE6Pj4OEVvIBASMIZCKyEnOR2YJ+IxkPsicr5WVFdre3j4XzNfTGI1pBh0BAQkEpIksvHGpVIrOOuYzkbu7u8\/PSBZ9830PHz6kvb09evDgAa2trYHIEorBrUBABoFURC4UCrS5uRkRk4nc39\/fNLxmrzw1NdWUyOyxQw+th4aG6PXr18HLyYaXF1ll5dRt49aI3NPTQwsLC8QhNi4gEDoC6+vrtLGxoU3MVEQuFovULrQWI27mkfl3JjP\/wwUEQkeAPbJOryxNZAY8abJLJLwahdahKw7yAQGTCKQicnz5qVKpnCe6eL7M19zc3LkMrTyySUHRFxAIGYFURA4ZEMgGBHxEQDuRZ2ZmaHx8PMKGl6Li3joOGIfrAwMD0Z+2traipS3friSyxqOZDx8+fLG+7pvM8fE2K\/7xWSYee19fHy0tLdHu7m5Du4zr3ZZOtRI5rlgGZGJighYXF2l\/f79Gt\/G1aA7FR0ZGaH5+vmm1mIuGkVRWfmHxNTk52XbpzkU5m42pVfGPT3LUjzX+4m3kYOpfXmzLvBLDKzL1dq4TB61E5jfV8PBwJNTbt28bvtX4bffkyRN6+fJlTcGITqF1tJ1EVu43\/tJqVkyjY3w622xX\/KOzb51ts20+fvyYdnZ2oloIsVLTqk9+oTVzWDrHqp3Ig4ODUbEIXxyeHB4e1oTXgsj8+61bt7wNrZnI7WQVimQC379\/v+VUQ6fSdbUdamKzvpqxFX62Xs7WiSxAKpfLEcHjns1kaJLVuJMQWcy1xMusXVVc1jGZfj7vRI7bwNHRkVH4tRM5SWgdTyT4agxJQuv6+ZSvsraaJ4dYM5DEI9vyxEIXWomcNAEUB8FXj5xE1nqP7KusIHItAo3qJ4y6YyLSSmQWJp6aj2f94sILA+\/t7SVb6XsVwCeRNeTlp9AiDGETjTyysN+Dg4PoU93Lly+fm9DJyUlYWWsV5EAbQAAItEdAu0duPwTcAQSAQFYEQOSsCOJ5IOAAAiCyA0rAEIBAVgT+BwIOFZDyH\/5qAAAAAElFTkSuQmCC","height":146,"width":242}}
%---
