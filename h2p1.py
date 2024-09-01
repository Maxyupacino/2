from numpy import ones,copy,cos,tan,pi,linspace,exp
import matplotlib.pyplot as plt

def gaussxw(N):
    # Initial approximation to roots of the Legendre polynomial
    a = linspace(3,4*N-1,N)/(4*N+2)
    x = cos(pi*a+1/(8*N*N*tan(a)))

    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta > epsilon:
        p0 = ones(N,float)
        p1 = copy(x)
        for k in range(1,N):
            p0,p1 = p1,((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))

    #Calculate the weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)
    return x,w

V = 0.001
kb = 1.38e-23
p = 6.022e28
xitaD = 428

y = []
for T in range(5,500):
    def Cv(x):
        return 9*V*p*kb*((T/xitaD)**3)*(x**4)*exp(x)/((exp(x)-1)**2)
    N = 50
    a = 0.0
    b = xitaD/T

#Calculate the sample points and weights,then map them
#to the required integration domain
    x,w = gaussxw(N)
    xp = 0.5*(b-a)*x + 0.5*(b+a)
    wp = 0.5*(b-a)*w

#Perform the integration
    s = 0
    for k in range(N):
        s += wp[k]*Cv(xp[k])
    y.append(s)

#以1K为最小单位画出所对应的各个温度的Cv并且连线
plt.plot(range(5,500),y,color='blue')
plt.xlabel('T/K')
plt.ylabel('Cv/(J/mol*K)')
plt.grid(True)
plt.legend()
plt.show()