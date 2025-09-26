
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# Plotagem de malha de estruturas reticuladas 2D (barras, vigas, treliças e pórticos)
def plotframe(coord, inci):
    plt.plot(coord[:,0],coord[:,1],'ro', markersize=5)
    for i in range(coord.shape[0]):
        plt.annotate(f'{i}', xy=(coord[i,0], coord[i,1]), xytext=(3,3), textcoords='offset points', fontsize=8)
    for i in range(inci.shape[0]):
        x = np.array([coord[inci[i,1],0],coord[inci[i,2],0]])
        y = np.array([coord[inci[i,1],1],coord[inci[i,2],1]])
        plt.plot(x,y)
    plt.show()
    return

# Plotagem da estrutura deformada e indeformada (treliça e pórtico)
def plotdeform(tp, coord, inci, Le, T, pos, desloc):
    n = inci.shape[0]
    npts = 20
    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.1, bottom=0.25)
    # Configurando a posição do slider
    ax_slider = plt.axes([0.1, 0.1, 0.8, 0.05])
    slider = Slider(ax_slider, 'Escala', 0.1, 100.0, valinit=1.0)
    # Plotagem da estrutura indeformada
    undefor = []
    defor = []
    for i in range(n):
        x = np.array([coord[inci[i,1],0],coord[inci[i,2],0]])
        y = np.array([coord[inci[i,1],1],coord[inci[i,2],1]])
        undefor.append(ax.plot(x,y,color='blue',linestyle='--')[0])
    # PLotagem da estrutura deformada
    if tp == 'truss':
        for i in range(n):
            x = np.array([coord[inci[i,1],0] + desloc[2*inci[i,1]],coord[inci[i,2],0] + desloc[2*inci[i,2]]])
            y = np.array([coord[inci[i,1],1] + desloc[2*inci[i,1]+1],coord[inci[i,2],1] + desloc[2*inci[i,2]+1]])
            defor.append(ax.plot(x,y,color='red')[0])
    else:
        for i in range(n):
            c = T[i][0,0]
            s = T[i][0,1]
            if c != 0:
                m = s/c
                b = coord[inci[i,1],1] - m * coord[inci[i,1],0]
                x = np.linspace(coord[inci[i,1],0],coord[inci[i,2],0],npts)
                y = m * x + b
            else:
                y = np.linspace(coord[inci[i,1],1],coord[inci[i,2],1],npts)
                x = np.linspace(coord[inci[i,1],0],coord[inci[i,2],0],npts)
            Ue = np.zeros((6))
            for j in range(6):
                Ue[j] = desloc[pos[i][j]]
            ue = np.dot(T[i],Ue)
            for j in range(npts):
                xi = np.linspace(-1.,1.,npts)
                v = ue[1] * (1/2-3/4*xi[j]+1/4*xi[j]**3) + ue[2] * Le[i]/8*(1-xi[j]-xi[j]**2+xi[j]**3)
                v += ue[4] * (1/2+3/4*xi[j]-1/4*xi[j]**3) + ue[5] * Le[i]/8*(-1-xi[j]+xi[j]**2+xi[j]**3)
                u = ue[0] * (1-xi[j])/2 + ue[3] * (1+xi[j])/2
                dx = c * u - s * v
                dy = s * u + c * v
                x[j] += dx
                y[j] += dy
            defor.append(ax.plot(x,y,color='red')[0])                  
    # Função de atualização do gráfico
    def update(val):
        scale = slider.val
        for i, line in enumerate(defor):
            if tp == 'truss':
                x = np.array([coord[inci[i,1],0] + scale * desloc[2*inci[i,1]],coord[inci[i,2],0] + scale * desloc[2*inci[i,2]]])
                y = np.array([coord[inci[i,1],1] + scale * desloc[2*inci[i,1]+1],coord[inci[i,2],1] + scale * desloc[2*inci[i,2]+1]])
            else:
                c = T[i][0,0]
                s = T[i][0,1]
                if c != 0:
                    m = s/c
                    b = coord[inci[i,1],1] - m * coord[inci[i,1],0]
                    x = np.linspace(coord[inci[i,1],0],coord[inci[i,2],0],npts)
                    y = m * x + b
                else:
                    y = np.linspace(coord[inci[i,1],1],coord[inci[i,2],1],npts)
                    x = np.linspace(coord[inci[i,1],0],coord[inci[i,2],0],npts)
                Ue = np.zeros((6))
                for j in range(6):
                    Ue[j] = desloc[pos[i][j]]
                ue = np.dot(T[i],Ue)
                for k in range(npts):
                    xi = np.linspace(-1.,1.,npts)
                    v = ue[1] * (1/2-3/4*xi[k]+1/4*xi[k]**3) + ue[2] * Le[i]/8*(1-xi[k]-xi[k]**2+xi[k]**3)
                    v += ue[4] * (1/2+3/4*xi[k]-1/4*xi[k]**3) + ue[5] * Le[i]/8*(-1-xi[k]+xi[k]**2+xi[k]**3)
                    u = ue[0] * (1-xi[k])/2 + ue[3] * (1+xi[k])/2
                    dx = c * u - s * v
                    dy = s * u + c * v
                    x[k] += dx * scale 
                    y[k] += dy * scale 
            line.set_xdata(x)
            line.set_ydata(y)
        ax.relim() 
        ax.autoscale_view()
        fig.canvas.draw_idle()
    # Ligando o slider à função de atualização
    slider.on_changed(update)
    # Apresentado o resultado
    plt.show()
    return

# Plotagem de resposta transiente
def plot_transi(arquivo,gdl,variavel,metodo):
    tplinha = ['-'] # '-' solid, '--' tracejado, '-.' traço ponto , ':' pontilhado
    cores = ['blue']
    ext = '.des'
    if variavel == 'velocidade':
        ext = '.vel'
    else:
        if variavel == 'aceleracao':
            ext = '.ace'
    tempo = np.genfromtxt(arquivo + metodo + ext, delimiter=' ', skip_header=1, usecols = 0)
    var = np.genfromtxt(arquivo + metodo + ext, delimiter=' ', skip_header=1, usecols = gdl)
    plt.plot(tempo, var, linestyle=tplinha[0], linewidth=1.2, color=cores[0])
    plt.grid()
    plt.xlabel('tempo')
    plt.ylabel(variavel)
    plt.show()
    return

# PLotagem de malha 2D formada por elementos triangulares ou quadrilaterais
def plot_malha2D(coord, inci):
    plt.scatter(coord[:, 0], coord[:, 1], color='red')
    for i in range(coord.shape[0]):
        plt.annotate(f'{i}', xy=(coord[i,0], coord[i,1]), xytext=(3,3), textcoords='offset points', fontsize=8)
    for q in inci[:,1:]:
        plt.fill(coord[q,0], coord[q,1], edgecolor='blue', fill=False)
    plt.show()
    return

# Plotagem da estrutura deformada e indeformada (Estado Plano)
def plotdef2D(coord, inci, desloc):
    scale = 1.
    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.1, bottom=0.25)
    # Configurando a posição do slider
    ax_slider = plt.axes([0.1, 0.1, 0.8, 0.05])
    slider = Slider(ax_slider, 'Escala', 0.1, 100.0, valinit=scale)
    defor = []
    undefor = []
    # Plotagem da estrutura indeformada
    for q in inci[:,1:]:
        undefor.append(ax.fill(coord[q,0], coord[q,1], edgecolor='blue', linestyle='--', fill=False))
    # PLotagem da estrutura deformada
    for q in inci[:,1:]:
        defor.append(ax.fill(coord[q,0] + desloc[2*q]*scale, coord[q,1] + desloc[2*q+1]*scale, edgecolor='red', fill=False))               
    # Função para atualizar o gráfico
    def update(val):
        scale = slider.val
        ax.clear()
        undefor.clear()
        defor.clear()
        for q in inci[:,1:]:
            undefor.append(ax.fill(coord[q,0], coord[q,1], edgecolor='blue', linestyle='--', fill=False))
        for q in inci[:,1:]:
            defor.append(ax.fill(coord[q,0] + desloc[2*q]*scale, coord[q,1] + desloc[2*q+1]*scale, edgecolor='red', fill=False))
        ax.relim() 
        ax.autoscale_view()
        fig.canvas.draw_idle()
    # Conectando os sliders à função de atualização
    slider.on_changed(update)
    # Apresentado o resultado
    plt.show()
    return

# Gráfico de contorno 2D
def plot_contor2D(coord,result,text):
    # Criando o gráfico de contorno
    plt.figure()
    contour = plt.tricontourf(coord[:,0], coord[:,1], result, cmap='Reds')
    plt.title(text)
    plt.xlabel('X')
    plt.ylabel('Y')
    # Adicionando a barra de cores
    plt.colorbar(contour)
    # Encontrando os valores máximos e mínimos
    z_max = np.max(result)
    z_min = np.min(result)
    # Encontrando os índices dos valores máximos e mínimos
    max_index = np.argmax(result)
    min_index = np.argmin(result)
    # Coordenadas dos pontos máximos e mínimos
    x_max, y_max = coord[max_index,0], coord[max_index,1]
    x_min, y_min = coord[min_index,0], coord[min_index,1]
    # Adicionando os pontos máximos e mínimos ao gráfico
    plt.plot(x_max, y_max, 'go', label=f'Máx: {z_max:.6f}')
    plt.plot(x_min, y_min, 'bo', label=f'Mín: {z_min:.6f}')
    # Exibindo a legenda
    plt.legend(framealpha=0.3)
    # Exibindo o gráfico
    plt.show()