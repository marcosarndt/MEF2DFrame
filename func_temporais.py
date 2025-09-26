# Funções temporais para descrição das forças variando com o tempo em Análises Transientes

import numpy as np

def heaviside(t,arg):
    f = 1
    return f

def senoidal(t,freq):
    f = np.sin(freq*t)
    return f

def cosenoidal(t,freq):
    f = np.cos(freq*t)
    return f

def forca_tempo(F,t,tipofunc,arg):
    opcoes = {'hav': heaviside, 'sen': senoidal, 'cos': cosenoidal}
    Ft = opcoes.get(tipofunc)(t,arg) * F
    return Ft
