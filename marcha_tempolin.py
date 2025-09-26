# Métodos de Marcha no Tempo para Problemas da Dinâmica Linear
import numpy as np
import func_temporais as ft
from scipy import linalg

def Newmark(gama,beta,K,M,C,F,u0,up0,tipo_ft,argft,dt,T,arquivo):
    # abertura e preparação dos arquivos de saída
    arq_desloc = open(arquivo+'.des', 'w', encoding='utf-8')
    arq_veloc = open(arquivo+'.vel', 'w', encoding='utf-8')
    arq_acel = open(arquivo+'.ace', 'w', encoding='utf-8')
    arq_desloc.write(f'Deslocamentos - Arquivo: {arquivo} ; Método: Newmark gama= {gama} beta= {beta} ; Delta t: {dt}\n')
    arq_veloc.write(f'Velocidades - Arquivo: {arquivo} ; Método: Newmark gama= {gama} beta= {beta} ; Delta t: {dt}\n')
    arq_acel.write(f'Acelerações - Arquivo: {arquivo} ; Método: Newmark gama= {gama} beta= {beta} ; Delta t: {dt}\n')

    # Cálculo da aceleração inicial
    t=0.
    F0 = ft.forca_tempo(F,t,tipo_ft,argft)
    upp0 = linalg.solve(M,F0-np.dot(C,up0)-np.dot(K,u0))

    # Inclusão das condições iniciais nos arquivos de saída
    arq_desloc.write(f'{t} {' '.join(str(numero) for numero in u0)}\n')
    arq_veloc.write(f'{t} {' '.join(str(numero) for numero in up0)}\n')
    arq_acel.write(f'{t} {' '.join(str(numero) for numero in upp0)}\n')

    # Cálculo dos parâmetros constantes
    a1 = (1/(beta*(dt**2)))*M + (gama/(beta*dt))*C
    a2 = (1/(beta*dt))*M + ((gama/beta)-1)*C
    a3 = ((1/(2*beta))-1)*M + dt*((gama/(2*beta))-1)*C
    Kc = K + a1
    
    # Processo iterativo
    while t <= T:
        t += dt
        F1 = ft.forca_tempo(F,t,tipo_ft,argft)
        F1c = F1 + np.dot(a1,u0) + np.dot(a2,up0) + np.dot(a3,upp0)
        u1 = linalg.solve(Kc,F1c)
        up1 = (gama/(beta*dt))*(u1-u0) + (1-(gama/beta))*up0 + dt*(1-(gama/(2*beta)))*upp0
        upp1 = (1/(beta*(dt**2)))*(u1-u0) - (1/(beta*dt))*up0 - ((1/(2*beta))-1)*upp0
        # atualização dos arquivos de saída
        arq_desloc.write(f'{t} {' '.join(str(numero) for numero in u1)}\n')
        arq_veloc.write(f'{t} {' '.join(str(numero) for numero in up1)}\n')
        arq_acel.write(f'{t} {' '.join(str(numero) for numero in upp1)}\n')
        # atualização das variáveis i
        u0 = u1
        up0 = up1
        upp0 = upp1
    
    # Fechamento dos arquivos de saída
    arq_desloc.close()
    arq_veloc.close()
    arq_acel.close()
    return

def difcentral(K,M,C,F,u0,up0,tipo_ft,argft,dt,T,arquivo):
    # abertura e preparação dos arquivos de saída
    arq_desloc = open(arquivo+'.des', 'w', encoding='utf-8')
    arq_veloc = open(arquivo+'.vel', 'w', encoding='utf-8')
    arq_acel = open(arquivo+'.ace', 'w', encoding='utf-8')
    arq_desloc.write(f'Deslocamentos - Arquivo: {arquivo} ; Método: Diferença Central; Delta t: {dt}\n')
    arq_veloc.write(f'Velocidades - Arquivo: {arquivo} ; Método: Diferença Central ; Delta t: {dt}\n')
    arq_acel.write(f'Acelerações - Arquivo: {arquivo} ; Método: Diferença Central ; Delta t: {dt}\n')

    # Cálculo da aceleração inicial
    t=0.
    F0 = ft.forca_tempo(F,t,tipo_ft,argft)
    upp0 = linalg.solve(M,F0-np.dot(C,up0)-np.dot(K,u0))

    # Inclusão das condições iniciais nos arquivos de saída
    arq_desloc.write(f'{t} {' '.join(str(numero) for numero in u0)}\n')
    arq_veloc.write(f'{t} {' '.join(str(numero) for numero in up0)}\n')
    arq_acel.write(f'{t} {' '.join(str(numero) for numero in upp0)}\n')

    # Cálculo dos parâmetros constantes
    um1 = u0 - dt*up0 + ((dt**2)/2)*upp0
    Kc = (1/(dt**2))*M + (1/(2*dt))*C
    a = (1/(dt**2))*M - (1/(2*dt))*C
    b = K - (2/(dt**2))*M
    
    # Processo iterativo
    while t <= T:
        Fc = F0 - np.dot(a,um1) - np.dot(b,u0)
        u1 = linalg.solve(Kc,Fc)
        up0 = (1/(2*dt))*(u1 - um1)
        upp0 = (1/(dt**2))*(u1 - 2*u0 + um1)
        t += dt
        # atualização dos arquivos de saída
        arq_desloc.write(f'{t} {' '.join(str(numero) for numero in u1)}\n')
        if t >= 2*dt:
            arq_veloc.write(f'{t-dt} {' '.join(str(numero) for numero in up0)}\n')
            arq_acel.write(f'{t-dt} {' '.join(str(numero) for numero in upp0)}\n')
        # atualização das variáveis
        um1 = u0
        u0 = u1
        F0 = ft.forca_tempo(F,t,tipo_ft,argft)
    
    # Fechamento dos arquivos de saída
    arq_desloc.close()
    arq_veloc.close()
    arq_acel.close()
    return

def marcha_tempo(metodo,K,M,C,F,u0,up0,tipo_ft,argft,dt,T,arquivo):
    opcoes = {'DC': difcentral, 'NAC': Newmark, 'NAL': Newmark}
    if metodo == 'NAC' or metodo == 'NAL':
        if metodo == 'NAC':
            gm = 1./2.
            bt = 1./4.
        else:
            gm = 1./2.
            bt = 1./6.
        opcoes.get(metodo)(gm,bt,K,M,C,F,u0,up0,tipo_ft,argft,dt,T,arquivo+metodo)
    else:
        opcoes.get(metodo)(K,M,C,F,u0,up0,tipo_ft,argft,dt,T,arquivo+metodo)
    return

# K = np.array([[6.,-2.],[-2.,4.]])
# M = np.array([[2.,0.],[0.,1.]])
# C = np.zeros((2,2))
# F = np.array([0.,10.])
# u0 = np.array([0.,0.])
# up0 = np.array([0.,0.])
# Newmark(1./2.,1./4.,K,M,C,F,u0,up0,'hav',0.,0.05,20,'testeN')
# difcentral(K,M,C,F,u0,up0,'hav',0.,0.05,20,'testeDC')