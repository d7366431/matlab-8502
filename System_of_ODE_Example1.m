t = 0:0.1:10;
n=3
y_init = [ones(n,1); zeros(n,1)]
[t,y] = ode45(@rhs,t,y_init);
plot(t,y(:,end)); 

function dydt = rhs(t,y)
n=3
A=2
f=1
        y1     = zeros(n,1)
        y2     = zeros(n,1)
        dydt_1 = zeros(n,1)
        dydt_2 = zeros(n,1)
        dydt   = zeros(2*n,1)
        y1 = y(1:n)
        y2 = y(n+1:2*n)
        dydt_1 = y2
        dydt_2 = A*y1 + 2*y2 + f
        dydt = [dydt_1; dydt_2]
  end