import pandas as pd
import numpy as np

import matplotlib

import matplotlib.pyplot as plt

import Vector as vctr
import Quaternion as qtrn

import cufflinks as cf

import plotly
import plotly.offline as py
import plotly.graph_objs as go

class BiVector:
    def __init__(self, vector1, vector2):
        self.u = vector1.unit()
        self.v = vector2.unit()

    def __str__(self):
        return "u = {}  ;    v = {} ".format(self.u,self.v)


    def contents(self):
        print("u = ",self.u , ";    v = ",self.v)

    def array(self):
        array = np.array([[self.u.array()] ,[self.v.array()]])
        return array

    def Rot_by_Mat(self,Rot_Mat):
        u2 = vctr.array2vec(Rot_Mat@(self.u.array())).unit()
        v2 = vctr.array2vec(Rot_Mat@(self.v.array())).unit()
        Bvctr2 = BiVector(u2,v2)
        return Bvctr2

#   FUNCTION
def Det4(u1,u2,v1,v2):
    det = vctr.Dot(u1.unit(),v2.unit()) - vctr.Dot(u2.unit(),v1.unit())
    return -np.sign(det)

def CoPlanar_Index(u1,u2,v1,v2):
    a = (vctr.ScTrProduct(u1,v1,v2))**2 + (vctr.ScTrProduct(u2,v1,v2))**2
    return np.sign(a)

def Axis_Angle(Bvctr1,Bvctr2):
    u1 = Bvctr1.u.unit()
    u2 = Bvctr2.u.unit()
    v1 = Bvctr1.v.unit()
    v2 = Bvctr2.v.unit()

    det =  Det4(u1,u2,v1,v2)
    if det == 0:
        if CoPlanar_Index(u1,u2,v1,v2) == 0:
            if Det4(u1,v1,v2,u2) == 0:
                n = vctr.Cross(u1,u2)
                cosTheta = vctr.Dot(u1,u2)
                Theta = np.arccos(cosTheta)
            else:
                n = (u1.unit() + u2.unit() + v1.unit() + v2.unit()).unit()
                Theta = pi
        else:
            u1 = vctr.Cross(u1,v1).unit()
            u2 = vctr.Cross(u2,v2).unit()
            Bvctr1 = BiVector(u1,v1)
            Bvctr2 = BiVector(u2,v2)
            return Axis_Angle(Bvctr1,Bvctr2)
    else:
        n = vctr.Cross((u1-u2),(v1-v2))
        eta = vctr.Dot(n,n)
        u0 = vctr.Dot(u1,u2)
        un = vctr.Dot(n,u1)
        cosTheta = (u0*eta - un*un)/(eta - un*un)
        Theta = np.arccos(cosTheta)
        return det*n,Theta
