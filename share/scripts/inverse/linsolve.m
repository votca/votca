A = load('$name.gmc');
b = load('$name.imc');

I=eye(size(A));

x(:,1)=b(:,1);
x(:,2)=-(A'*inv((A'*A)+$reg*I))*b(:,2);

save '$name_out' x '-ascii'
quit
