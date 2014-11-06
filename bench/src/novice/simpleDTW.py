import numpy as np
import matplotlib.pyplot as plt
#matplotlib inline

a = np.array([1, 1, 2, 3, 2, 0])
b = np.array([0, 1, 1, 2, 3, 2, 1])

plt.plot(a,'r', label='a')
plt.plot(b, 'g', label='b')
plt.legend()

#plt.close()

distances = np.zeros((len(b), len(a)))

for i in range(len(b)):
    for j in range(len(a)):
        distances[i,j] = (a[j]-b[i])**2 

def distance_cost_plot(distances):
    im = plt.imshow(distances, interpolation='nearest', cmap='Reds') 
    plt.gca().invert_yaxis()
    plt.xlabel("A")
    plt.ylabel("B")
    plt.grid()
    plt.colorbar()

distance_cost_plot(distances)
#plt.close()
    

accumulated_cost = np.zeros((len(b), len(a)))

accumulated_cost[0,0] = distances[0,0]

distance_cost_plot(accumulated_cost)

for i in range(1, len(x)):
    accumulated_cost[0,i] = distances[0,i] + accumulated_cost[0, i-1] 

distance_cost_plot(accumulated_cost)

for i in range(1, len(y)):
    accumulated_cost[i,0] = distances[i, 0] + accumulated_cost[i-1, 0]    

distance_cost_plot(accumulated_cost)

for i in range(1, len(y)):
    for j in range(1, len(x)):
        accumulated_cost[i, j] = 
        min(accumulated_cost[i-1, j-1], accumulated_cost[i-1, j], 
            accumulated_cost[i  , j-1]) + distances[i, j]


path = [[len(x)-1, len(y)-1]]
i = len(y)-1
j = len(x)-1

while i>0 and j>0:
    print i,j
    if accumulated_cost[i-1,j] == min(accumulated_cost[i-1, j-1], 
                                      accumulated_cost[i-1, j], 
                                      accumulated_cost[i, j-1]):
        i = i-1
    elif accumulated_cost[i,j-1] == min(accumulated_cost[i-1, j-1], 
                                        accumulated_cost[i-1, j], 
                                        accumulated_cost[i, j-1]):
        j = j-1
    else:
        i = i-1
        j= j-1
    path.append([j,i])
    for [y, x] in path:
        cost = cost +distances[x, y]
    return path, cost    

path_x = [point[0] for point in path]
path_y = [point[1] for point in path]

distance_cost_plot(accumulated_cost)
plt.plot(path_x, path_y);
