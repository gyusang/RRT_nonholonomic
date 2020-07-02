import matplotlib.pyplot as plt
import numpy as np

from result import nodes_list, points

fig, ax = plt.subplots()

goal = plt.Circle((0.9, 0.9), 0.05, color='b')
obstacle = plt.Rectangle((0.3, 0.3), 0.4, 0.4, color='#00FFFF', fill=True)
background = plt.Rectangle((0, 0), 1, 1, color='#000000', lw=1, fill=False)
ax.add_artist(background)
ax.add_artist(goal)
ax.add_artist(obstacle)

for node in nodes_list:
    if node[1] is not None:
        ax.plot(node[3], node[4], 'k-')
    ax.plot([node[2][0]], [node[2][1]], 'k.')
ax.plot([nodes_list[0][2][0]], [nodes_list[0][2][1]], 'ro')

for node in points:
    if node[1] is not None:
        ax.plot(node[3], node[4], 'r-')
    ax.plot([node[2][0]], [node[2][1]], 'r.')

plt.axis('equal')
# plt.xlim((-0.01, 1.01))
# plt.ylim((-0.01, 1.01))
plt.show()