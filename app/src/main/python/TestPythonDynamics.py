#Purpose:This is for testing. It should output the xHat, and yHat arrays after taking in the subject data
#Author:Noah Lape
#3/2/2021

import mat4py as loadmat
import numpy as np
from numpy import random

import scipy as sci
from scipy import signal
from scipy.fft import fft, ifft
from scipy.special import comb

from os.path import dirname, join

import math as math
import scipy.io

def simulateDynamics(stateLength, t, y, Ad, Bd, Cd, L):
    xHat = np.zeros([stateLength,len(t)], dtype = float)
    xHat[:,0] = 1
    yHat = np.zeros([1,len(t)])
    L = np.matrix([[.00028762], [.00011015], [.00051347]])

    # first loop through columns,start one column in
    for j in range(1,len(t)):#columns

        #this array temporarily holds a single xHat column and allows for input into the matrix multiplication
        xHatTempColm = np.zeros((len(xHat),1)) #create a temp array for holding variables
        for i in range(0,len(xHat)):#rows
            xHatTempColm[i,0] = xHat[i,j-1]

        #Calculates the current state
        curXColm = np.matmul((Ad - L*Cd), xHatTempColm) + L*y[j-1]
        #sets the xHat colm
        xHat[:,[j]] = curXColm
        #obtains the yHat scalar
        yHat[0,[j]] = np.matmul(Cd,curXColm)

    #this is a matrix
    return [xHat,yHat]

def loadData(subject):
    t = []
    y = []
    tPlusY = []
    subjectData = []


    filename = join(dirname(__file__), "A3.mat")
    data = scipy.io.loadmat(filename)
    # data = scipy.io.loadmat('A' + str(subject)+ '.mat')
    #takes the data from dictionary "data"
    subjectData = data['Data']

    i = 0
    #loop goes until the last row of the matrix
    while i < len(subjectData):
        t.append(subjectData[i][0])
        y.append(subjectData[i][len(subjectData[0])-1])
        i = i + 1
    t = np.transpose(t)

    #converting to numpy arrays for performance
    t = np.array(t)
    y = np.array(y)

    tPlusY = [t,y]
    return tPlusY

def createStateSpace (order, stateLength, t, omg, zeta, gamma_d,gamma_omg):
    Ac = np.zeros([stateLength,stateLength], dtype = float)
    Bc = np.zeros([stateLength,1],dtype = float)
    Cc = np.zeros([1,stateLength],dtype = float)
    Dc = 0

    #This is to populate State matrices
    #   -Python takes the LHS which forms a 2x2 matrix, matches it to the RHS which, to python looks like a
    #    2x2 matrix and assigns each element accordingly. Process follows suite for Bc matrix
    for k in range(order):
        i = (1 + k) * 2
        k = k + 1 #handles the one off error
        Ac[i - 2:i, i - 2:i] = [[0, 1], [float(-(k * omg) ** 2), 0]]
        Bc[i-2:i] = [[0],[(k*omg)**2]]
        Cc[0][i-1] = (2 * zeta) / (k * omg)

    #Assigns edge values
    Bc[len(Bc)-1][0] = gamma_d
    Cc[0][len(Cc[0])-1] = 1

    #uses the "scipy" library and converts the continuous matrixes into a SS
    #   -Then takes the continuous system and decritizes it (NOTE: method for c2d may not match model)
    dt = float(t[1]-t[0])
    contSystem = signal.StateSpace(Ac, Bc, Cc, Dc)
    discSystem = contSystem.to_discrete(dt)

    #Posible solution if .to_discrete system does not work.
    #   -function for computing matrix exponenetal (Future implimentation)
    #   -might need to maintain stability
    #discSystem = signal.cont2discrete(contSystem,dt,'impulse')

    return discSystem

#initlial conditions
def TestPythonDynamics():
    pi = math.pi
    omg = (2*pi)/24
    zeta = 1
    gamma_d = 1
    gamma_omg = 0
    order = 1
    stateLength = (2*order + 1)
    subject = 3

    t = loadData(subject)[0] #loads time into t
    y = loadData(subject)[1] #loads data into y
    discSystem = createStateSpace(order, stateLength,t, omg, zeta,gamma_d,gamma_omg )

    Ad = discSystem.A
    Bd = discSystem.B
    Cd = discSystem.C
    Dc = discSystem.D

    #FOR TESTING. FORCING L ARRAY
    L = np.matrix([[.00028762],[.00011015],[.00051347]])

    xyHat = simulateDynamics(stateLength,t,y,Ad,Bd,Cd,L)
    xHat = xyHat[0]
    yHat = xyHat[1]
    result = "xHat: {} yHat: {}".format(xHat,yHat);
    return(result)
    # print("Finish")

