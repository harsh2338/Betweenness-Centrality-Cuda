import matplotlib.pyplot as plt
# line 1 points

x=[500,1000,2000,3000,4000,5000,10000,50000]
serial= [40.354,174.183,696.3,1598.29,2740.26,4462.15,18020.1,536997]
vertex= [95.316,233.816,529.076,1019.63,1527.83,2467.21,8782.2,207194]
edge= [55.596,178.061,403.688,670.56,1069.92,1719.09,6734.76,174811]
we= [114.307,266.456,439.739,684.583,908.601,1353.45,3717.3,65036.8]
wd= [24.413,25.735,34.552,49.149,66.664,95.955,269.268,15660.6]
serial=[x/1000.0 for x in serial]
vertex=[x/1000.0 for x in vertex]   
edge=[x/1000.0 for x in edge]
we=[x/1000.0 for x in we]
wd=[x/1000.0 for x in wd]
x=[i/1000 for i in x]

xt=[128,256,512,1024]
st=[4462.15,4462.15,4462.15,4462.15]
vt=[6987.55,4104.33,2828.54,2467.21]
et=[6385.65,3548.33,2241.53,1719.09]
wet=[2328.56,1690.65,1418.06,1353.45]
wdt=[134.329,80.017,74.527,95.955]



pgv=[]
pge=[]
pgwe=[]
pgwd=[]

suv=[]
sue=[]
suwe=[]
suwd=[]

for i in range(len(x)):
    suv=1/((1-vertex[i])+vertex[i]/xt[n])

for i in range(len(x)):
    pgv.append((serial[i]-vertex[i])*100/serial[i])

for i in range(len(x)):
    pge.append((serial[i]-edge[i])*100/serial[i])

for i in range(len(x)):
    pgwe.append((serial[i]-we[i])*100/serial[i])
for i in range(len(x)):
    pgwd.append((serial[i]-wd[i])*100/serial[i])

print(pgv)
print(pge)
print(pgwe)
print(pgwd)

[-136.1996332457749, -34.23583242911191, 24.01608502082435, 36.204944033936265, 44.24507163553823, 44.7080443284067, 51.26442139610767, 61.41617178494481]
[-37.77072904792585, -2.226394079789642, 42.02384029872181, 58.045160765568205, 60.955529767248365, 61.47395313918177, 62.626400519419974, 67.44655929176513]
[-183.2606433067354, -52.974744952148036, 36.846330604624434, 57.167785570828826, 66.84252589170372, 69.66820927131539, 79.37136863835384, 87.888796399235]
[39.50289934083362, 85.22530901408288, 95.03777107568577, 96.92490098793087, 97.5672381452855, 97.8495792387078, 98.50573526229043, 97.0836708584964]

# plt.plot(x, serial, label = "Serial")
# plt.plot(x, vertex, label = "Vertex Parallel")
# plt.plot(x, edge, label = "Edge Parallel")
# plt.plot(x, we, label = "Work efficient")
# plt.plot(x, wd, label = "Work Distributed")

# plt.xlabel('No of vertices and edges (x1000)')
# # Set the y axis label of the current axis.
# plt.ylabel('Run Time (s)')
# # Set a title of the current axes.
# plt.title('Run Time vs No of vertices and edges')
# show a legend on the plot

# plt.plot(xt, st, label = "Serial")
# plt.plot(xt, vt, label = "Vertex Parallel")
# plt.plot(xt, et, label = "Edge Parallel")
# plt.plot(xt, wet, label = "Work efficient")
# plt.plot(xt, wdt, label = "Work Distributed")

# plt.xlabel('No of threads per block')
# # Set the y axis label of the current axis.
# plt.ylabel('Run Time (ms)')
# # Set a title of the current axes.
# plt.title('Run Time for 5000 nodes and threads vs No of threads')



plt.plot(x, pgv, label = "Vertex Parallel")
plt.plot(x, pge, label = "Edge Parallel")
plt.plot(x, pgwe, label = "Work efficient")
plt.plot(x, pgwd, label = "Work Distributed")

plt.xlabel('No of vertices and edges (x1000)')
# Set the y axis label of the current axis.
plt.ylabel('Performance Gain (%)')
# Set a title of the current axes.
plt.title('Performance Gain wrt Serial')


plt.legend()
# Display a figure.
plt.show()