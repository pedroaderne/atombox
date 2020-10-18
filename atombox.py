
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation


epsilon = eval(input('ε em Kcal/mol: '))

sigma = eval(input('σ em Angstrom: '))

massa = eval(input('Massa do átomo em Kg: '))

delta_t = 0.001

num_atomos = 20

tamanho_da_caixa = 50

elev = 30
azim = 45

xlim = (0, tamanho_da_caixa)
ylim = (0, tamanho_da_caixa)
zlim = (0, tamanho_da_caixa)

fig = plt.figure()
ax = Axes3D(fig)

ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_zlim(zlim)

ax.view_init(elev, azim)


class Ball:
    def __init__(self, xyz, v, fmt, epsilon, sigma, massa):
        self.xyz = np.array(xyz)
        self.v = np.array(v)
        self.epsilon = epsilon
        self.sigma = sigma
        self.massa = massa
        self.aceleracao = [0,0,0]

        self.scatter, = ax.plot([], [], [], fmt, animated=True)



    def interage(self, atomo):
        coord2 = atomo.getcoord()
        direcao = atomo.getv()
        distancia = np.sqrt((self.xyz[0] - coord2[0])**2 + (self.xyz[1]- coord2[1])**2 + (self.xyz[2]- coord2[2])**2)


        def reduzPra1(n):
            if(n > 0):
                n = 0
            elif(n < 0):
                n = -1
            else:
                n = 0

            return (4 * epsilon * (((sigma/10**10) ** 12 / (distancia/10**10) ** 13) * 12 - ((sigma/10**10) ** 6 / (distancia/10**10) ** 7) * 6))*n





        if distancia < 2**(1/6)*sigma:
            Pot_LJ = 4 * epsilon * ((sigma / distancia) ** 12 - (sigma / distancia) ** 6)  # Fórmula para o potencial
            For_LJ = map(reduzPra1,direcao)  # Fórmula para a força
            for a, n  in enumerate(For_LJ):
                self.aceleracao[a] = n / self.massa
                print(self.aceleracao)


    def getcoord(self):
        return self.xyz

    def getv(self):
        return self.v

    def update(self):

        if self.xyz[0] <= xlim[0]:
            self.v[0] = np.abs(self.v[0])
            self.aceleracao[0] = np.abs(self.aceleracao[0])

        elif self.xyz[0] >= xlim[1]:
            self.v[0] =  -np.abs(self.v[0])
            self.aceleracao[0] = -np.abs(self.aceleracao[0])

        if self.xyz[1] <= ylim[0]:
            self.v[1] =  np.abs(self.v[1])
            self.aceleracao[1] = np.abs(self.aceleracao[1])

        elif self.xyz[1] >= ylim[1]:
            self.v[1] = -np.abs(self.v[1])
            self.aceleracao[1] = -np.abs(self.aceleracao[1])

        if self.xyz[2] <= zlim[0]:
            self.v[2] =  np.abs(self.v[2])
            self.aceleracao[2] = np.abs(self.aceleracao[2])

        elif self.xyz[2] >= zlim[1]:
            self.v[2] =  -np.abs(self.v[2])
            self.aceleracao[2] = -np.abs(self.aceleracao[2])

        delta_v = delta_t
        self.v += delta_v

        self.xyz += self.v


        self.xyz[0] = np.clip(self.xyz[0], xlim[0], xlim[1])
        self.xyz[1] = np.clip(self.xyz[1], ylim[0], ylim[1])
        self.xyz[2] = np.clip(self.xyz[2], zlim[0], zlim[1])


        self.scatter.set_xdata(self.xyz[0])
        self.scatter.set_ydata(self.xyz[1])
        self.scatter.set_3d_properties(self.xyz[2])

        return self.xyz



atomos = []

for i in np.arange(0, num_atomos):
    xyz = np.random.rand(1, 3)[0] * tamanho_da_caixa
    v = np.random.rand(1, 3)[0] * 0.1
    cor_do_atomo = str('b' + 'o')
    atomos.append(Ball(xyz, v, cor_do_atomo, epsilon, sigma, massa))


def init():
    return []


def update(t):
    global elev, azim

    for atomo in atomos:
        for atomo2 in atomos:
            if atomo2 == atomo:
                continue
            atomo.interage(atomo2)


        atomo.update()

    artists = [atomo.scatter for atomo in atomos]

    return artists


ani = FuncAnimation(fig, update, frames=np.arange(0, 10, delta_t), init_func=init, interval=10, blit=True,
                    repeat=False)

plt.show()