# Leitura de arquivo .msh criado pelo Gmsh ou texto .mef e .loa. Montagem das matrizes de coordenadas nodais, incidência de elementos e materiais.
# Listas de condições de contorno e carregamentos
# Versão 4 - 21/11/2024

import numpy as np

# Converte linha do arquivo texto em lista de números inteiros (tipo = 1) ou float (tipo diferente de 1)
def conv_listnum (tipo, linha):
    listastr = linha.split()
    linhac = []
    for i in range(len(listastr)):
        if tipo == 1:
            linhac.append(int(listastr[i]))
        else:
            linhac.append(float(listastr[i]))
    return(linhac)

# Converte linha do arquivo texto em lista de números inteiros (tipo = 1) ou float (tipo diferente de 1) removendo o primeiro item
def conv_listnum2 (tipo, linha):
    listastr = linha.split()
    listastr.pop(0)
    linhac = []
    for i in listastr:
        if tipo == 1:
            linhac.append(int(i))
        else:
            linhac.append(float(i))
    return(linhac)

# Gera uma lista de nós não repetidos
def lista_nos(dados):
    nos = []
    for i in range(dados.shape[0]):
        for j in range(dados.shape[1]):
            nos.append(dados[i,j])
    nos_sem_repeticao = list(set(nos))
    return(nos_sem_repeticao)

# Leitura de malha para um tipo específico de elemento diretamente de arquivo .msh gerado pelo GMSH
def read_gmsh(arquivo,tpe):
    # Tipo de elemento tpe: (1) barra 2 nós, (2) triangular 3 nós ou (3) quadrangular 4 nós
    try:
        file = open(arquivo, 'r', encoding='utf-8')
        linha = (file.readline())
        while linha != '$Nodes\n':
            linha = (file.readline())
        info_nodes = conv_listnum(1,file.readline())
        print(f'Número de blocos de nós: {info_nodes[0]}')
        print(f'Número de nós: {info_nodes[1]}')
        coord = []
        for i in range(info_nodes[0]):
            print(f'Bloco {i}')
            info_bloco = conv_listnum(1,file.readline())
            for j in range(info_bloco[3]):
                file.readline()
            for j in range(info_bloco[3]):
                coord.append(conv_listnum(2,file.readline()))
        # cada linha de coordenadas corresponde às coordenadas [x, y, z] de um nó
        coordenadas = np.array(coord)
        print(coordenadas)
        print(f'A matriz de coordenadas possui {coordenadas.shape[0]} linhas e {coordenadas.shape[1]} colunas.')
        while linha != '$Elements\n':
            linha = (file.readline())
        info_elem= conv_listnum(1,file.readline())
        print(f'Número de blocos de elementos: {info_elem[0]}')
        print(f'Número de elementos: {info_elem[1]}')
        elem = []
        for i in range(info_elem[0]):
            print(f'Bloco {i}')
            info_bloco = conv_listnum(1,file.readline())
            print(f'Tipo de elemento: (1) barra 2 nós, (2) triangular 3 nós ou (3) quadrangular 4 nós: {info_bloco[2]}')
            for j in range(info_bloco[3]):
                linha = conv_listnum(1,file.readline())
                if info_bloco[2] == tpe:
                    linha[0] = 0
                    for k in range(info_bloco[2]+1):
                        linha[k+1] = linha[k+1] - 1
                    # print(linha)
                    elem.append(linha)
        # print(elem)
        # A matriz incidência em cada linha as informações de um elemento na forma [material, nó 1, nó 2, ...]
        incidencia = np.array(elem)
        print(incidencia)
        print(f'A matriz incidencia possui {incidencia.shape[0]} linhas e {incidencia.shape[1]} colunas.')
    except IOError:
        print(f'Não foi possível abrir o arquivo!')
    else:
        file.close()
    return(coordenadas,incidencia)

# Leitura de malha de treliças e pórticos 2D a partir de arquivo texto .mef
def read_geo(arquivo):
    try:
        file = open(arquivo, 'r', encoding='utf-8')
        linha = (file.readline())
        while linha != 'Nodes\n':
            linha = (file.readline())
        coord = []
        linha = (file.readline())
        while linha.split()[0] != 'end':
            coord.append(conv_listnum(2,linha))
            linha = (file.readline())
        # cada linha de coordenadas corresponde às coordenadas [x, y, z] de um nó
        coordenadas = np.array(coord)
        print(coordenadas)
        print(f'A matriz de coordenadas possui {coordenadas.shape[0]} linhas e {coordenadas.shape[1]} colunas.')
        while linha != 'Elements\n':
            linha = (file.readline())
        elem = []
        linha = (file.readline())
        while linha.split()[0] != 'end':
            elem.append(conv_listnum(1,linha))
            linha = (file.readline())
        # A matriz incidência em cada linha as informações de um elemento na forma [material, nó 1, nó 2, ...]
        incidencia = np.array(elem)
        print(incidencia)
        print(f'A matriz incidencia possui {incidencia.shape[0]} linhas e {incidencia.shape[1]} colunas.')
    except IOError:
        print(f'Não foi possível abrir o arquivo!')
    else:
        file.close()
    return(coordenadas,incidencia)

# Leitura de material, condições de contorno (apoios e carregamentos) e informações sobre a análise
# de treliças e pórticos 2D de arquivo texto com extensão .loa
def read_load(arquivo):
    try:
        file = open(arquivo, 'r', encoding='utf-8')
        linha = (file.readline())
        while linha != 'Dados Analise\n':
            linha = (file.readline())
        linha = (file.readline())
        tp = linha.split()
        while linha != 'Material\n':
            linha = (file.readline())
        mat = []
        linha = (file.readline())
        while linha.split()[0] != 'end':
            mat.append(conv_listnum(2,linha))
            linha = (file.readline())
        # cada linha de materiais corresponde às propriedades [densidade, área, mód. de elasticidade, inércia]
        materiais = np.array(mat)
        print(materiais)
        print(f'A matriz de materiais possui {materiais.shape[0]} linhas e {materiais.shape[1]} colunas.')
        while linha != 'C Contorno\n':
            linha = (file.readline())
        contorno = []
        linha = (file.readline())
        while linha.split()[0] != 'end':
            listastr = linha.split()
            linhac = []
            for i in range(len(listastr)):
                if i <= 1:
                    linhac.append(int(listastr[i]))
                else:
                    linhac.append(float(listastr[i]))
            contorno.append(linhac)
            linha = file.readline()
        # A lista contorno contém em cada posição [nó, direção, valor prescrito]
        print(contorno)
        while linha != 'F nodais\n':
            linha = (file.readline())
        fnodais = []
        linha = (file.readline())
        while linha.split()[0] != 'end':
            listastr = linha.split()
            linhac = []
            for i in range(len(listastr)):
                if i == 0:
                    linhac.append(int(listastr[i]))
                else:
                    linhac.append(float(listastr[i]))
            fnodais.append(linhac)
            linha = file.readline()
        # A lista fnodais contém em cada posição [nó, Fx, Fy, Mz]
        print(fnodais)
        while linha != 'F elementares\n':
            linha = (file.readline())
        felem = []
        linha = (file.readline())
        while linha.split()[0] != 'end':
            listastr = linha.split()
            linhac = []
            for i in range(len(listastr)):
                if i == 0:
                    linhac.append(int(listastr[i]))
                else:
                    linhac.append(float(listastr[i]))
            felem.append(linhac)
            linha = file.readline()
        # A lista felem contém em cada posição [elemento, px, qyi, qyj]
        print(felem)
        while linha != 'Apoios elasticos\n':
            linha = (file.readline())
        apoioel = []
        linha = (file.readline())
        while linha.split()[0] != 'end':
            listastr = linha.split()
            linhac = []
            for i in range(len(listastr)):
                if i == 0:
                    linhac.append(int(listastr[i]))
                else:
                    linhac.append(float(listastr[i]))
            apoioel.append(linhac)
            linha = file.readline()
        # A lista apoioel contém em cada posição [nó, kx, ky, ktheta]
        print(apoioel)
        dini = []
        vini = []
        if tp[1] == 'transiente':
            while linha != 'C inicial desloc\n':
                linha = (file.readline())
            linha = (file.readline())
            while linha.split()[0] != 'end':
                listastr = linha.split()
                linhac = []
                for i in range(len(listastr)):
                    if i == 0:
                        linhac.append(int(listastr[i]))
                    else:
                        linhac.append(float(listastr[i]))
                dini.append(linhac)
                linha = file.readline()
            # A lista dini contém em cada posição [nó, desloc inicial em x, desloc inicial em y, rot inicial em torno de z]
            print(dini)
            while linha != 'C inicial vel\n':
                linha = (file.readline())
            linha = (file.readline())
            while linha.split()[0] != 'end':
                listastr = linha.split()
                linhac = []
                for i in range(len(listastr)):
                    if i == 0:
                        linhac.append(int(listastr[i]))
                    else:
                        linhac.append(float(listastr[i]))
                vini.append(linhac)
                linha = file.readline()
            # A lista vini contém em cada posição [nó, veloc inicial em x, veloc inicial em y, veloc rot inicial em torno de z]
            print(vini)
    except IOError:
        print(f'Não foi possível abrir o arquivo!')
    else:
        file.close()
    return(tp,materiais,contorno,fnodais,felem,apoioel,dini,vini)

# Leitura de malha 2D e regiões com condições de contorno geradas no Gmsh utilizando a biblioteca meshio
def read_mesh2D(path):
    # Dados da malha gmsh
    import meshio
    mesh = meshio.read(path + '.msh')
    # Obtenção da matriz de coordenadas nodais
    coord = mesh.points
    nn = coord.shape[0]
    print(f'Coordenadas dos {nn} nós: {coord}')
    # Obtenção das chaves da condições de contorno e seus números
    field_data = mesh.field_data
    list_cc = []
    for name in field_data.keys():
        numf = field_data[name]
        list_cc.append([name,numf[0],numf[1]])
    # Obtenção da matriz de incidência e condições de contorno
    cells = mesh.cells
    gfisicos = mesh.cell_data['gmsh:physical']
    elemdata = None
    mat = 0
    i = 0
    condc = {}
    material = []
    for cell_block in cells:
        if cell_block.type == "triangle" or cell_block.type == "quad":
            if mat == 0:
                elemtype = cell_block.type
                elemdata = cell_block.data
                aux = np.full((elemdata.shape[0], 1), mat)
                inci = np.hstack((aux, elemdata))
                mat += 1
            else:
                elemdata = cell_block.data
                aux = np.full((elemdata.shape[0], 1), mat)
                aux = np.hstack((aux, elemdata))
                inci = np.vstack((inci,aux))
                mat += 1
            for j in list_cc:
                if j[1] == gfisicos[i][0]:
                    material.append(j[0])
        else:
            for j in list_cc:
                if j[1] == gfisicos[i][0]:
                    condc[j[0]] = [cell_block.type,cell_block.data]   
        i += 1
    ne = inci.shape[0]
    print(f'Incidência dos {ne} elementos do tipo {elemtype}: {inci}')
    print(f'Materiais: {material}')
    print(f'Condições de contorno: {condc}')
    return(coord,inci,elemtype,material,condc)

# Leitura de materiais, condições de contorno (apoios e carregamentos) e informações sobre a análise
# de problemas de estado plano de arquivo texto com extensão .mcc
def read_mccEP(path,lista_mat,condc,ngl,coord):
    try:
        file = open(path + '.mcc', 'r', encoding='utf-8')
        linha = (file.readline())
        while linha != 'Dados Analise\n':
            linha = (file.readline())
        linha = (file.readline())
        tp = linha.split()
        while linha != 'Material\n':
            linha = (file.readline())
        mat = {}
        linha = (file.readline())
        while linha.split()[0] != 'end':
            mat[linha.split()[0]] = conv_listnum2(2,linha)
            linha = (file.readline())
        # cada linha de materiais corresponde às propriedades [espessura, mód. de elasticidade, coef. de Poisson e densidade]
        list_mat = []
        for chave in lista_mat:
            list_mat.append(mat[chave])
        materiais = np.array(list_mat)
        print(materiais)
        print(f'A matriz de materiais possui {materiais.shape[0]} linhas e {materiais.shape[1]} colunas.')
        while linha != 'C Contorno\n':
            linha = (file.readline())
        contorno = []
        desl = np.zeros((ngl))
        linha = (file.readline())
        while linha.split()[0] != 'end':
            listastr = linha.split()
            dir = int(listastr[1])
            data = condc[listastr[0]][1]
            if condc[listastr[0]][0] == 'vertex':
                if dir != 2:
                    contorno.append(2*data[0,0]+dir)
                    desl[2*data[0,0]+dir] = float(listastr[2])
                else:
                    contorno.append(2*data[0,0])
                    contorno.append(2*data[0,0]+1)
                    desl[2*data[0,0]] = float(listastr[2])
                    desl[2*data[0,0]+1] = float(listastr[3])
            else:
                if condc[listastr[0]][0] == 'line':
                    nos = lista_nos(data)
                    if dir != 2:
                        for i in range(len(nos)):
                            contorno.append(2*nos[i]+dir)
                            desl[2*nos[i]+dir] = float(listastr[2])
                    else:
                        for i in range(len(nos)):
                            contorno.append(2*nos[i])
                            contorno.append(2*nos[i]+1)
                            desl[2*nos[i]] = float(listastr[2])
                            desl[2*nos[i]+1] = float(listastr[3])
            linha = file.readline()
        contorno = list(set(contorno))
        # A lista contorno contém os graus de liberdade prescritos
        # O vetor desl contém os deslocamentos prescritos e os demais nulos
        force = np.zeros((ngl))
        while linha != 'F nodais\n':
            linha = (file.readline())
        linha = (file.readline())
        while linha.split()[0] != 'end':
            listastr = linha.split()
            data = condc[listastr[0]][1]
            if condc[listastr[0]][0] == 'vertex':
                force[data[0,0]*2] = float(listastr[1])
                force[data[0,0]*2+1] = float(listastr[2])
            linha = file.readline()
        # O vetor force até aqui contém as forças nodais relacionadas a cada grau de liberdade
        while linha != 'F superficie\n':
            linha = (file.readline())
        linha = (file.readline())
        while linha.split()[0] != 'end':
            listastr = linha.split()
            qx = float(listastr[1])
            qy = float(listastr[2])
            data = condc[listastr[0]][1]
            if condc[listastr[0]][0] == 'line':
                for i in range(data.shape[0]):
                    noi = data[i,0]
                    noj = data[i,1]
                    aresta = np.sqrt((coord[noj,0]-coord[noi,0])**2 + (coord[noj,1]-coord[noi,1])**2)
                    force[noi*2] += qx * aresta / 2
                    force[noj*2] += qx * aresta / 2
                    force[noi*2+1] += qy * aresta / 2
                    force[noj*2+1] += qy * aresta / 2
            linha = file.readline()
        # O vetor force até aqui contém as forças nodais e de superfícies relacionadas a cada grau de liberdade
    except IOError:
        print(f'Não foi possível abrir o arquivo!')
    else:
        file.close()
    return(tp,materiais,contorno,desl,force)