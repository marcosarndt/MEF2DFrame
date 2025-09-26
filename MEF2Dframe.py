import numpy as np
from scipy import linalg
import ler_dados
import marcha_tempolin as mtemp
import plotmef2d as pltmef

# Função para encontrar o índice da sublista cujo primeiro elemento é o valor procurado
def pos_lista(lista_de_listas, valor):
    for i, sublista in enumerate(lista_de_listas):
        if sublista[0] == valor:
            return i
        return -1 # Retorna -1 se o valor não for encontrado como primeiro elemento de nenhuma sublista
    
# Contrução das Matrizes de Transformação, Comprimentos elementares e Posições
def transf(tipoe, coord, inci):
    L = []
    T = []
    pos = []
    for i in range(inci.shape[0]):
        noi = inci[i,1]
        nof = inci[i,2]
        xi = coord[noi,0]
        xf = coord[nof,0]
        yi = coord[noi,1]
        yf = coord[nof,1]
        Le = ((xf - xi)**2 + (yf - yi)**2)**(1/2)
        L.append(Le)
        c = (xf - xi)/Le
        s = (yf - yi)/Le
        if tipoe == 'truss':
            Tr = np.array([[c,s,0,0],[0,0,c,s]])
            posicao = [2*noi, 2*noi+1, 2*nof, 2*nof+1]
        if tipoe == 'frame':
            Tr = np.array([[c,s,0,0,0,0],[-s,c,0,0,0,0],[0,0,1,0,0,0],[0,0,0,c,s,0],[0,0,0,-s,c,0],[0,0,0,0,0,1]])
            posicao = [3*noi, 3*noi+1, 3*noi+2, 3*nof, 3*nof+1, 3*nof+2]
        T.append(Tr)
        pos.append(posicao)
    return(L,T,pos)

# Entrada de dados
# Tipo de entrada de dados
arquivo = input('Qual o nome do arquivo de modelo que você deseja abrir? ')
tpin = int(input('Tipo de entrada de dados de malha (0 para arquivo texto e 1 para leitura de malha Gmsh): '))

# Entrada de dados da malha
if tpin == 1: 
    Coord,Inci = ler_dados.read_gmsh(arquivo + '.msh',1) # Entrada via Gmsh
else:
    Coord,Inci = ler_dados.read_geo(arquivo + '.mef') # Entrada via arquivo de texto

# Entrada de dados de materiais, contorno e carregamento
tp,Mat,Contor,Fnod,Felem,ApEl,Dinic,Vinic = ler_dados.read_load(arquivo + '.loa')

tpelem = tp[0]
tpana = tp[1]
# Número de graus de liberdade por nó
if tpelem == 'truss':
    gln = 2
else:
    gln = 3

# Vetor de forças nodais no sistema global
# Fnodal [gdl, força ou momento]
Fnodal = []
for i in range(len(Fnod)):
    gl = gln*Fnod[i][0]
    for j in range(gln):
        Fnodal.append([gl+j,Fnod[i][j+1]])
print(Fnodal)

# Vetor de forças elementares no sistema local (carga uniformemente distribuída em x e trapezoidal em y)
# Felem [elemento, valor da carga na direção x, valor da carga na direção y no nó1, valor da carga na direção y no nó2]

# Condições de contorno [gdl, deslocamento prescrito]
CC = []
for i in range(len(Contor)):
    CC.append([Contor[i][0] * gln + Contor[i][1], Contor[i][2]])

if tpana == 'transiente':
    # Vetor de deslocamentos iniciais no sistema global para análise transiente
    # Dinic_gl [gdl, deslocamento ou rotação]
    Dinic_gl = []
    for i in range(len(Dinic)):
        gl = gln*Dinic[i][0]
        for j in range(gln):
            Dinic_gl.append([gl+j,Dinic[i][j+1]])
    print(Dinic_gl)
    # Vetor de velocidades iniciais no sistema global para análise transiente
    # Vinic_gl [gdl, velocidade ou veloc. angular]
    Vinic_gl = []
    for i in range(len(Vinic)):
        gl = gln*Vinic[i][0]
        for j in range(gln):
            Vinic_gl.append([gl+j,Vinic[i][j+1]])
    print(Vinic_gl)

# Plotagem da malha
pltmef.plotframe(Coord,Inci)

# Dimensão do problema e ordem de reordenação do sistema
nnos = Coord.shape[0] # número de nós
nelem = Inci.shape[0] # número de elementos
numgdl = gln*nnos
ncon = len(CC) # número de condições de contorno
dim = numgdl - ncon
ordem = [i for i in range(numgdl)]
for i in range(ncon):
    ordem.remove(CC[i][0])
    ordem.append(CC[i][0])

# Matrizes de Rigidez e Massa
Kg = np.zeros((numgdl,numgdl))
Mg = np.zeros((numgdl,numgdl))
list_K = []
list_M = []
Le,T,pos = transf(tpelem,Coord,Inci)
for i in range(nelem):
    nmat = Inci[i,0]
    I = Mat[nmat,3] # momento de inércia
    E = Mat[nmat,2] # módulo de elasticidade
    A = Mat[nmat,1] # area da seção
    ro = Mat[nmat,0] # densidade
    # Matriz de rigidez do elemento no sistema local
    if tpelem == 'truss':
        Ke = E*A/Le[i] * np.array([[1.,-1.],[-1.,1.]])
    else:
        if tpelem == 'frame':
            Ke = np.array([[E*A/Le[i],0,0,-E*A/Le[i],0,0],[0,12*E*I/(Le[i]**3),6*E*I/(Le[i]**2),0,-12*E*I/(Le[i]**3),6*E*I/(Le[i]**2)],
                           [0,6*E*I/(Le[i]**2),4*E*I/Le[i],0,-6*E*I/(Le[i]**2),2*E*I/Le[i]],[-E*A/Le[i],0,0,E*A/Le[i],0,0],
                           [0,-12*E*I/(Le[i]**3),-6*E*I/(Le[i]**2),0,12*E*I/(Le[i]**3),-6*E*I/(Le[i]**2)],
                           [0,6*E*I/(Le[i]**2),2*E*I/Le[i],0,-6*E*I/(Le[i]**2),4*E*I/Le[i]]])
    list_K.append(Ke)
    # Matriz de massa do elemento no sistema local
    if tpelem == 'truss':
        Me = ro*A*Le[i] * np.array([[1/3,1/6],[1/6,1/3]])
    else:
        if tpelem == 'frame':
            Me = ro*A*Le[i]/420. * np.array([[420/3,0,0,420/6,0,0],[0,156,22*Le[i],0,54,-13*Le[i]],[0,22*Le[i],4*Le[i]**2,0,13*Le[i],-3*Le[i]**2],
                                             [420/6,0,0,420/3,0,0],[0,54,13*Le[i],0,156,-22*Le[i]],[0,-13*Le[i],-3*Le[i]**2,0,-22*Le[i],4*Le[i]**2]])
    list_M.append(Me)
    # Matrizes elementares no sistema global
    Keg = np.dot(np.transpose(T[i]),np.dot(Ke,T[i]))
    Meg = np.dot(np.transpose(T[i]),np.dot(Me,T[i]))
    # Inclusão nas matrizes globais
    dimm = 2*gln
    for j in range(dimm):
        for k in range(dimm):
            Kg[pos[i][j],pos[i][k]] += Keg[j,k]
            Mg[pos[i][j],pos[i][k]] += Meg[j,k]

# Inclusão da rigidez dos apoios elásticos
for i in range(len(ApEl)):
    gdl = ApEl[i][0] * gln
    for j in range(gln):
        Kg[gdl+j,gdl+j] += ApEl[i][j+1]
    
# Vetor de forças
F = np.zeros((numgdl))
# Inclusão das forças elementares
temfelem = [0] * nelem
for i in range(len(Felem)):
    L = Le[Felem[i][0]]
    temfelem[Felem[i][0]] = 1
    if tpelem == 'truss':
        Fe = np.array([Felem[i][1]*L/2,Felem[i][1]*L/2])
    else:
        if tpelem == 'frame':
            Fe = np.array([Felem[i][1]*L/2, 7*Felem[i][2]*L/20 + 3*Felem[i][3]*L/20, Felem[i][2]*L**2/20 + Felem[i][3]*L**2/30,
                   Felem[i][1]*L/2, 3*Felem[i][2]*L/20 + 7*Felem[i][3]*L/20, -Felem[i][2]*L**2/30 - Felem[i][3]*L**2/20])
    # Vetor elementar no sistema global
    Feg = np.dot(np.transpose(T[Felem[i][0]]),Fe)
    for j in range(dimm):
        F[pos[Felem[i][0]][j]] += Feg[j]
# Inclusão das forças aplicadas globais
for i in range(len(Fnodal)):
    F[Fnodal[i][0]] += Fnodal[i][1]

# Vetor de deslocamentos
d = np.zeros((numgdl))
# Graus de liberdade prescritos
for i in range(ncon):
    d[CC[i][0]] = CC[i][1]

if tpana == 'transiente':
    # Vetores de deslocamentos e velocidades iniciais para análise transiente
    u0g = np.zeros((numgdl))
    up0g = np.zeros((numgdl))
    for i in range(len(Dinic_gl)):
        u0g[Dinic_gl[i][0]] += Dinic_gl[i][1]
    for i in range(len(Vinic_gl)):
        up0g[Vinic_gl[i][0]] += Vinic_gl[i][1]

# Reordenação das matrizes
Kgord = np.zeros((numgdl,numgdl))
Mgord = np.zeros((numgdl,numgdl))
Ford = np.zeros((numgdl))
dord = np.zeros((numgdl))
if tpana == 'transiente':
    u0gord = np.zeros((numgdl))
    up0gord = np.zeros((numgdl))
for i in range(numgdl):
   Ford[i] = F[ordem[i]]
   dord[i] = d[ordem[i]]
   if tpana == 'transiente':
       u0gord[i] = u0g[ordem[i]]
       up0gord[i] = up0g[ordem[i]]
   for j in range(numgdl):
       Kgord[i,j] = Kg[ordem[i],ordem[j]]
       Mgord[i,j] = Mg[ordem[i],ordem[j]]

# Sub Matrizes de rigidez e massa para solução do problema estático (K u = F), modal ou transiente (M a + C v + K u = F)
# com inclusão das condições de contorno
K = Kgord[:dim,:dim]
M = Mgord[:dim,:dim]
Fl = Ford[:dim]

if tpana == 'estatica':
    dl = linalg.solve(K,Fl - np.dot(Kgord[:dim,dim:],dord[dim:]))
    for i in range(dim):
        d[ordem[i]] = dl[i]
    print (f'Deslocamentos: {d}')
    # Cálculo das reações de apoio
    Fp = np.dot(Kgord[dim:,:dim],dl) + np.dot(Kgord[dim:,dim:],dord[dim:]) - Ford[dim:]
    print (f'Reacoes de apoio: {Fp}')
    # Esforços Internos
    if tpelem == 'truss':
        print('Esforços internos nodais: [Ni,  Nj]')
    else:
        print('Esforços internos nodais: [Ni, Qi, Mi, Nj, Qj, Mj]')
    for i in range(nelem):
        Ue = np.zeros((dimm))
        for j in range(dimm):
            Ue[j] = d[pos[i][j]]
        ue = np.dot(T[i],Ue)
        if temfelem[i] == 1:
            k = pos_lista(Felem,i)
            L = Le[Felem[k][0]]
            if tpelem == 'truss':
                Fe = np.array([Felem[k][1]*L/2,Felem[k][1]*L/2])
            else:
                if tpelem == 'frame':
                    Fe = np.array([Felem[k][1]*L/2, 7*Felem[k][2]*L/20 + 3*Felem[k][3]*L/20, Felem[k][2]*L**2/20 + Felem[k][3]*L**2/30,
                                Felem[k][1]*L/2, 3*Felem[k][2]*L/20 + 7*Felem[k][3]*L/20, -Felem[k][2]*L**2/30 - Felem[k][3]*L**2/20])
        else:
            if tpelem == 'truss':
                Fe = np.zeros((2))
            else:
                if tpelem == 'frame':
                    Fe = np.zeros((6))
        EsfInt = np.dot(list_K[i], ue) - Fe
        if tpelem == 'truss':
            EsfInt[0] = EsfInt[0] * (-1)
        else:
            if tpelem == 'frame':
                c = T[i][0,0]
                s = T[i][0,1]
                if (c>=0 and s>=0) or ((c>0 and s<0) and abs(c)>abs(s)) or ((c<0 and s>0) and abs(c)<abs(s)):
                    signei = np.array([-1.,1.,-1.,1.,-1.,1.])
                else:
                    signei = np.array([-1.,1.,1.,1.,-1.,-1.])
                for m in range(6):
                    EsfInt[m] = signei[m] * EsfInt[m]
        print(f'Elemento {i}: {EsfInt}')
    resp = input('Deseja plotar a estrutura deformada (S ou N)? ')
    if resp == 's' or resp == 'S':
        pltmef.plotdeform(tpelem, Coord, Inci, Le, T, pos, d)

if tpana == 'modal' or (tpana == 'transiente' and float(tp[3]) > 0):
    autovalores, autovetores = linalg.eigh(K, M)
    freq = autovalores ** (1/2)
    print(f'Frequências: {freq}')
    # Normalização dos modos em relação à massa
    Mmodal = np.array(autovetores)
    modo = []
    for i in range(dim):
        Mn = np.dot(np.transpose(Mmodal[:,i]),np.dot(M,Mmodal[:,i]))
        modo.append(1/(Mn**(1/2)) * Mmodal[:,i])
        print(f'Modo {i}: {modo[i]}')
    resp = 'S'
    while resp == 's' or resp == 'S':
        resp = input('Deseja plotar algum modo de vibração da estrutura (S ou N)? ')
        if resp == 's' or resp == 'S':
            dmodo = np.zeros((numgdl))
            nmodo = int(input(f'Modo que deseja plotar (0 a {dim-1}):'))
            for i in range(dim):
                dmodo[ordem[i]] = modo[nmodo][i]
            pltmef.plotdeform(tpelem, Coord, Inci, Le, T, pos, dmodo)

if tpana == 'transiente':
    txam = float(tp[3])
    met_int = tp[2]
    tp_force = tp[4]
    arg_force = float(tp[5])
    delta_t = float(tp[6])
    tempo_final = float(tp[7])
    if txam > 0:
        omega1 = freq[0]
        if dim == 1:
            omega2 = freq[0]
        else:
            omega2 = freq[1]
        a0 = txam*2*omega1*omega2/(omega1+omega2)
        a1 = txam*2/(omega1+omega2)
        C = a0 * M + a1 * K
    else:
        C = np.zeros((dim,dim))
    u0 = u0gord[:dim]
    up0 = up0gord[:dim]
    mtemp.marcha_tempo(met_int,K,M,C,Fl,u0,up0,tp_force,arg_force,delta_t,tempo_final,arquivo)
    resp = 'S'
    while resp == 's' or resp == 'S':
        resp = input('Deseja plotar a resposta transiente da estrutura (S ou N)? ')
        if resp == 's' or resp == 'S':
            var = input(f'Qual variável (deslocamento, velocidade ou aceleracao)? ')
            ntran = int(input(f'Qual o nó (0 a {nnos-1}): '))
            dir = int(input(f'Qual a direção (0 para x, 1 para y ou 2 para rotação): '))
            gld = gln * ntran + dir
            try:
                # Encontrando o índice do valor 
                indice = ordem.index(gld)
                pltmef.plot_transi(arquivo,indice,var,met_int)
            except ValueError: 
                print(f'O grau de liberdade é restrito.')
            
