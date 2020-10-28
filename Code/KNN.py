import numpy as np
import pandas as pd

def CDML(x, y):
    M = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]]) #Mudar isso depois
    sub = np.array(x-y)
    result = np.sum(sub * M * sub.transpose())
    return result