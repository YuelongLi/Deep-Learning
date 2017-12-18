import numpy as np

def initializeNetwork(layerSizes = [1,2,3,4,1]):
    l = len(layerSizes)
    parameters = {}
    for i in range(1,l):
        parameters['W'+str(i)] = np.random.randn(layerSizes[i],layerSizes[i-1])*0.1
        parameters['b'+str(i)] = np.empty((i,1))
    return parameters

def forwardProp(X, parameters):
    As = {}
    A = X
    l = len(parameters)//2
    for i in range(1, l):
        A = np.dot(parameters['W'+str(i)],A)
        As['A'+str(i)] = A
    return As

parameters = initializeNetwork()
forwardProp(np.array([[1,2,3,1]]),parameters)
