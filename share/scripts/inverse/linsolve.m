A = load('%NAME.gmc');
b = load('%NAME.imc');

b(:,2)=b(:,2)
x = b;
x(:,2)=0;

% x(25:end,2) = linsolve(A(25:end,25:end),b(25:end,2));
x(25:end-1,2) = 0.1*linsolve(A(25:end-1,25:end-1),b(25:end-1,2));

save '%NAME.dpot.new' x '-ascii'
quit