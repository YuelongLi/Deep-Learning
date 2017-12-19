import numpy as np

def initializeNetwork(layerSizes = [1,2,3,4,1]):
    l = len(layerSizes)
    parameters = {}
    for i in range(1,l):
        parameters['W'+str(i)] = np.random.randn(layerSizes[i],layerSizes[i-1])
        parameters['b'+str(i)] = np.empty((layerSizes[i],1))
    return parameters

def forwardProp(X, parameters):
    As = {}
    A = X
    l = len(parameters)//2+1
    for i in range(1, l):
        A = relu(np.dot(parameters['W'+str(i)],A)+parameters['b'+str(i)])
        As['A'+str(i)] = A
    return As

def getGradients(As, parameters, Y):
    gradients = {}

def sigmoid(Z):
    return 1/(1+np.exp(-Z))

def dSigmoid(A):
    return A-A**2

def relu(Z):
    return np.maximum(0,Z)

def dRelu(A):
    return 0 if A<0 else 1

parameters = initializeNetwork()
forwardProp(np.array([[1,2,3,1]]),parameters)

x = np.arange(-10,10,0.05)
X = x.reshape(1,-1)

result = forwardProp(X, parameters)
Y = result['A4'][0]

import matplotlib.pyplot as plt

plt.plot(x, Y)
plt.show()
