%% Laplace Transform Lab: Solving ODEs using Laplace Transform in MATLAB
%
% This lab will teach you to solve ODEs using a built in MATLAB Laplace 
% transform function |laplace|.
%
% There are five (5) exercises in this lab that are to be handed in.  
% Write your solutions in a separate file, including appropriate descriptions 
% in each step.
%
% Include your name and student number in the submitted file.
%
%% Student Information
%
%  Student Name: Yoonho Kim
%
%  Student Number: 1008035635
%

%% Using symbolic variables to define functions
% 
% In this exercise we will use symbolic variables and functions.

syms t s x y

f = cos(t)
h = exp(2*x)


%% Laplace transform and its inverse

% The routine |laplace| computes the Laplace transform of a function

F=laplace(f)

%%
% By default it uses the variable |s| for the Laplace transform
% But we can specify which variable we want:

H=laplace(h)
laplace(h,y)

% Observe that the results are identical: one in the variable |s| and the
% other in the variable |y|

%% 
% We can also specify which variable to use to compute the Laplace
% transform:

j = exp(x*t)
laplace(j)
laplace(j,x,s)

% By default, MATLAB assumes that the Laplace transform is to be computed
% using the variable |t|, unless we specify that we should use the variable
% |x|

%% 
% We can also use inline functions with |laplace|. When using inline
% functions, we always have to specify the variable of the function.

l = @(t) t^2+t+1
laplace(l(t))

%% 
% MATLAB also has the routine |ilaplace| to compute the inverse Laplace
% transform

ilaplace(F)
ilaplace(H)
ilaplace(laplace(f))

%% 
% If |laplace| cannot compute the Laplace transform, it returns an
% unevaluated call.

g = 1/sqrt(t^2+1)
G = laplace(g)

%% 
% But MATLAB "knows" that it is supposed to be a Laplace transform of a
% function. So if we compute the inverse Laplace transform, we obtain the
% original function

ilaplace(G)

%%
% The Laplace transform of a function is related to the Laplace transform 
% of its derivative:

syms g(t)
laplace(diff(g,t),t,s)


%% Exercise 1
%
% Objective: Compute the Laplace transform and use it to show that MATLAB
% 'knows' some of its properties.
%
% Details:  
%
% (a) Define the function |f(t)=exp(2t)*t^3|, and compute its Laplace
%   transform |F(s)|.
% (b) Find a function |f(t)| such that its Laplace transform is
%   |(s - 1)*(s - 2))/(s*(s + 2)*(s - 3)|
% (c) Show that MATLAB 'knows' that if |F(s)| is the Laplace transform of
%   |f(t)|, then the Laplace transform of |exp(at)f(t)| is |F(s-a)| 
% 
% (in your answer, explain part (c) using comments).      
%
% Observe that MATLAB splits the rational function automatically when
% solving the inverse Laplace transform.

% Part a)
syms f F L s t
f = exp(2*t) * t^3;
F = laplace(f);

% Part b)
L = ((s-1)*(s-2))/(s*(s+2)*(s-3));
f = ilaplace(L);

% Part c)
syms a f(t)
F_new = laplace(exp(a*t)*f(t))

% Since we stated that F is the laplace transform of f, we can see that
% MATLAB understood that taking the laplace of f multiplied by exp(a*t)
% would translate F by a, since F_new's result was laplace(f(t), t, s-a).

%% Heaviside and Dirac functions
%
% These two functions are builtin to MATLAB: |heaviside| is the Heaviside
% function |u_0(t)| at |0|
%
% To define |u_2(t)|, we need to write

f=heaviside(t-2)
ezplot(f,[-1,5])

% The Dirac delta function (at |0|) is also defined with the routine |dirac|

g = dirac(t-3)

% MATLAB "knows" how to compute the Laplace transform of these functions

laplace(f)
laplace(g)


%% Exercise 2
%
% Objective: Find a formula comparing the Laplace transform of a 
%   translation of |f(t)| by |t-a| with the Laplace transform of |f(t)|
%
% Details:  
%
% * Give a value to |a|
% * Let |G(s)| be the Laplace transform of |g(t)=u_a(t)f(t-a)| 
%   and |F(s)| is the Laplace transform of |f(t)|, then find a 
%   formula relating |G(s)| and |F(s)|
%
% In your answer, explain the 'proof' using comments.

syms f(t) g(t) G(s) F(s)

num = 10;
for j = 1:num
    a = j;
    g(t) = heaviside(t-a) * f(t-a);
    F = laplace(f(t));
    G = laplace(g(t));
    fprintf("a = %g\n", a);
    disp(F);
    disp(G);
end

% Through this step, we can see that G(s) = exp(-a*s) * F(s), which is what
% we expected since in lecture we discussed that a translation by a in
% cartesian coordinates would result in the multiplication by a factor of
% exp(-a*s) in the Laplace coordinates.

%% Solving IVPs using Laplace transforms
%
% Consider the following IVP, |y''-3y = 5t| with the initial
% conditions |y(0)=1| and |y'(0)=2|.
% We can use MATLAB to solve this problem using Laplace transforms:

% First we define the unknown function and its variable and the Laplace
% tranform of the unknown

syms y(t) t Y s

% Then we define the ODE

ODE=diff(y(t),t,2)-3*y(t)-5*t == 0

% Now we compute the Laplace transform of the ODE.

L_ODE = laplace(ODE)

% Use the initial conditions

L_ODE=subs(L_ODE,y(0),1)
L_ODE=subs(L_ODE,subs(diff(y(t), t), t, 0),2)

% We then need to factor out the Laplace transform of |y(t)|

L_ODE = subs(L_ODE,laplace(y(t), t, s), Y)
Y=solve(L_ODE,Y)

% We now need to use the inverse Laplace transform to obtain the solution
% to the original IVP

y = ilaplace(Y)

% We can plot the solution

ezplot(y,[0,20])

% We can check that this is indeed the solution

diff(y,t,2)-3*y


%% Exercise 3
%
% Objective: Solve an IVP using the Laplace transform
%
% Details: Explain your steps using comments
%
%
% * Solve the IVP
% *   |y'''+2y''+y'+2*y=-cos(t)|
% *   |y(0)=0|, |y'(0)=0|, and |y''(0)=0|
% * for |t| in |[0,10*pi]|
% * Is there an initial condition for which |y| remains bounded as |t| goes to infinity? If so, find it.

syms y(t) Y t s

func1 = diff(y(t), t, 1);
func2 = diff(y(t), t, 2);
func3 = diff(y(t), t, 3);

ODE = func3 + 2*func2 + func1 + 2*y(t) == -1*cos(t);

% Initial Conditions
y0_zero = 0; 
y1_zero = 0;
y2_zero = 0;

L = laplace(ODE);
L = subs(L, y(0), y0_zero);
L = subs(L, subs(func1, t, 0), y1_zero);
L = subs(L, subs(func2, t, 0), y2_zero);

L = subs(L, laplace(y(t), t, s), F);
Y = solve(ODE, Y);

y = ilaplace(Y);

% For the above ODE, we can find that the general solution is:
% y = (3*sin(t))/50 - (2*cos(t))/25 + (t*cos(t))/10 + (4*y(0)*cos(t))/5 - (t*sin(t))/5 + (2*y(0)*sin(t))/5 + exp(-2*t)*(y(0)/5 + 2/25)
% where all the terms either remain constant or are decaying. This is true
% except for the two terms (t*cos(t))/10 and (t*sin(t))/5. Since both of
% these terms are unbounded as t goes to infinity and do not depend on the
% initial conditions, this would imply that the entire solution would be
% unbounded, regardless of what initial condition is put on the solution.
% Therefore, there is no initial condition which can keep y bounded as t
% goes to infinity.
%% Exercise 4
%
% Objective: Solve an IVP using the Laplace transform
%
% Details:  
% 
% * Define 
% *   |g(t) = 3 if 0 < t < 2|
% *   |g(t) = t+1 if 2 < t < 5|
% *   |g(t) = 5 if t > 5|
%
% * Solve the IVP
% *   |y''+2y'+5y=g(t)|
% *   |y(0)=2 and y'(0)=1|
%
% * Plot the solution for |t| in |[0,12]| and |y| in |[0,2.25]|.
%
% In your answer, explain your steps using comments.

syms y(t) Y t s

g = @(t) (-t+4)*heaviside(t-5) + (t-2)*heaviside(t-2) + 3*heaviside(t);
diff1 = diff(y(t), t, 1);
diff2 = diff(y(t), t, 2);

ODE = diff2 + 2*diff1 + 5*y(t) == g(t);
y0_zero = 2;
y1_zero = 1;

Lapl = laplace(ODE);
Lapl = subs(Lapl, y(0), y0_zero);
Lapl = subs(Lapl, subs(diff1, t, 0), y1_zero);

Lapl = subs(Lapl, laplace(y(t), t, s), Y);
Y = solve(Lapl, Y);

y = ilaplace(Y);

% plot the solution
t = linspace(0, 12, 300);
y_graph = subs(y);
plot(t, y_graph);
ylim([0, 2.25]);
title("Solution for Laplace Piecewise ODE");
xlabel('t');
ylabel('y');

%% Exercise 5
%
% Objective: Use the Laplace transform to solve an integral equation
% 
% Verify that MATLAB knowns about the convolution theorem by explaining why the following transform is computed correctly.
syms t tau y(tau) s
I=int(exp(-2*(t-tau))*y(tau),tau,0,t)
laplace(I,t,s)

% The convolution theorem states that the laplace transform of the
% convolution of f(t) and g(t) is simply the product F(s)G(s). We know that
% the laplace of exp(-2*t) to be 1/(s+2). MATLAB computes the convolution
% of exp(-2*t) and some function y(t) as the product of the respective
% laplace transforms, seen through the result which is laplace(y(t), t,
% s)/(s+2), which can be written as Y(s)/(s+2). This confirms that MATLAB
% knows about the convolution theorem.
