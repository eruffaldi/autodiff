
%% basic test
clear all
a = matexp(eye(3));
b = matexp('b',ones(3));
c = a*(2+b);
update(c)
value(c)
set(b,ones(3)*2);
update(c)
value(c)
%
vars = collectvars(c);
vars

autodiff(c)

%% example from paper
X = matexp('X',ones(3));
F = trace((inv(I+X)*X')*X);