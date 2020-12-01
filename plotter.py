from matplotlib import pyplot as plt
import numpy as np


def transform_array(a, Y, X):
    a_mod = []
    for i in range(X):
        a_mod.append([])
    for i in range(X):
        for j in range(Y):
            a_mod[i].append(a[j][i])
    return a_mod


def prepare_array(a, Y, X):
    a_pr = []
    for i in range(Y):
        a_pr.append([])
    average = 0.0
    N = 0
    for i in range(Y):
        for j in range(X):
            if a[i][j] != 0.0:
                average += a[i][j]
                N += 1
    average /= N
    for i in range(Y):
        for j in range(X):
            if a[i][j] <= average * 5:
                a_pr[i].append(a[i][j])
            else:
                a_pr[i].append(average)
    return a_pr


data = open("data.txt", 'r')
Y, X = map(int, data.readline().split())
a = []
for i in range(Y):
    a.append(list(map(float, data.readline().split())))


a_pr = prepare_array(a, Y, X)
a_mod = transform_array(a_pr, Y, X)

fig, ax = plt.subplots()
ax.imshow(a_mod, interpolation='none', cmap='OrRd', origin='lower')

fig.set_figwidth(6)    #  ширина и
fig.set_figheight(6)    #  высота "Figure"

plt.show()
data.close()
