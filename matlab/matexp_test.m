
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
X = matexp('X',magic(3));
F = trace((inv(eye(3)+X)*X')*X);
update(F)
value(F)
autodiff(F)
adjoint(X)
vX = value(X);
T1 = eye(3)+vX;
T2 = inv(T1);
T3 = vX';
T4 = T2*T3;
T5 = T4*vX;
% by example is: (R0 T4)' + ((XR0T2)')' - (-T2T3XR0T2)'
R0 = 1;
R1 = R0*vX;
R2 = R0*T4;
R3 = T3*R1;
R4 = R1*T2;
R5 = -T2*R3*T2;
R6 = R4';
R8 = R5;
Q = R2' + R6' + R8'; %R2 ok

value(F)
trace(T5)
aX = adjoint(X)
Q
%% example with sym
%clear all
Xs = sym(sym('x',[3,3]),'real');
X = matexp('X',Xs);
F = trace((inv(eye(length(Xs))+X)*X')*X);
update(F)
value(F)
autodiff(F)
adjoint(X)

%% then we can compute the jacobian directly from F
J = jacobian(value(F),Xs(:));
Xv = magic(3);
t = { Xs(:),Xv(:)};
Xc = cell(numel(Xv),1);
ts = Xc;
for I=1:length(ts)
    Xc{I} = t{2}(I);
    ts{I} = t{1}(I);
end
Jv = double(subs(J,ts,Xc));

assert(all(abs(Jv-Q(:)')<1e-3),'Jv corresponds to manual');
assert(all(abs(Jv-aX(:)')<1e-3),'Jv corresponds to auto');
