from numpy import loadtxt
import matplotlib.pyplot as plt


net_size = 3
subplots = (net_size, 2)
substances = ['A', 'B', 'C']
fig, ax = plt.subplots(subplots[0], subplots[1], sharex=True)

ax[0, 0].set_title('Abundance')
ax[0, 1].set_title('Absolute difference')
ax[0, 0].set_xscale('log')
ax[0, 1].set_xscale('log')

ax[-1, 0].set_xlabel('Time [s]')
ax[-1, 1].set_xlabel('Time [s]')
ax[-1, 0].set_xticks([10**e for e in range(-5, 6, 2)])


euler = loadtxt('data/sol0.txt')
xtfc = loadtxt('data/sol1.txt')

for i in range(subplots[0]):
    ax[i, 0].plot(euler[:, 0], euler[:, i + 1], color='red', label='Backward Euler')
    ax[i, 0].plot(xtfc[:, 0], xtfc[:, i + 1], linestyle='dashed', color='blue', label='X-TFC')
    ax[i, 0].set_ylabel(substances[i])
    ax[i, 1].plot(euler[:, 0], abs(euler[:, i + 1] - xtfc[:, i + 1]))
    ax[i, 0].grid()
    ax[i, 1].grid()
    ax[i, 1].set_yscale('log')


ax[0, 0].legend(['Backward Euler', 'X-TFC'])
plt.savefig('rober.pdf', bbox_inches='tight')
plt.show()


fig, ax = plt.subplots(net_size, 1, sharex=True)
data = loadtxt('data/sol2.txt')
for i in range(net_size):
    ax[i].plot(data[:, 0], data[:, i + 1], marker='.')
    ax[i].set_ylabel(f'$y_{substances[i]}$')
    ax[i].grid()
ax[0].set_title('X-TFC generalized')
plt.xscale('log')
plt.xlabel('Time [s]')
plt.show()


