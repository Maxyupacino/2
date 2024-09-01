from math import tanh
import matplotlib.pyplot as plt
import numpy as np
#分别从math库和matplotlib.pylot中import函数
def f(x):#引入函数f(x)
    return 1+0.5*tanh(2*x)

h = 0.25#设出间隔数
def df0(x):#定义用中心差值计算
    return (f(x+h/2)-f(x-h/2))/h

def df1(x):#定义用解析法计算
    return 1-tanh(2*x)*tanh(2*x)

cd = []#中心差值df0的列表
tr = []#解析法df1得到的列表
for i in range(-200,200):
    j = i/100
    tr1 = df1(j)
    tr.append(tr1)
    cd1 = df0(j)
    cd.append(cd1)

x = np.linspace(-2.00,2.00,400)
#标刻出x的刻度
plt.plot(x,cd,color='red',marker='.')
plt.plot(x,tr,color='black')
#画出df0和df1关于x的图像，其中df0打点方便辨认
plt.xlabel('x')
plt.ylabel('df')
plt.grid(True)
plt.legend()
plt.show()