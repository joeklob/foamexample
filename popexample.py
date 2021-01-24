
#For making a figure showing popping process for Voronoi and lattice initial conditions.  
#Outputs are initial Voronoi tesselation and foam after given number of ruptures


#Updated for 01/20/2021





import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi
from shapely.geometry import Polygon
import random



np.random.seed([193840])




#How many cells and how many pops
N = 250
numpops = 100









def voronoi_finite_polygons_2d(vor, radius=None):


    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()*2

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]
        


        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)







# # For generating Voronoi diagram
points = np.random.uniform(0, 1, size=(N,2))      


# compute Voronoi tesselation
vor = Voronoi(points)

# plot
regions, vertices = voronoi_finite_polygons_2d(vor)


min_x = 0
max_x = 1
min_y = 0
max_y = 1

mins = np.tile((min_x, min_y), (vertices.shape[0], 1))
bounded_vertices = np.max((vertices, mins), axis=0)
maxs = np.tile((max_x, max_y), (vertices.shape[0], 1))
bounded_vertices = np.min((bounded_vertices, maxs), axis=0)


box = Polygon([[min_x, min_y], [min_x, max_y], [max_x, max_y], [max_x, min_y]])

# colorize
q = 0
polys = []
goodvertices = []
polyind = [[] for i in range(len(regions))]
p = 0
edgeind = []

#If you want a pretty picture of initial conditions, use this
for region in regions:    
    polygon = vertices[region]
    # Clipping polygon
    poly = Polygon(polygon)
    poly = poly.intersection(box)
    polygon = [p for p in poly.exterior.coords]
    for v in polygon[0:-1]:
        if v in goodvertices:
            polyind[p].append(goodvertices.index(v))
        else:
            goodvertices.append(v)
            polyind[p].append(q)
            q+= 1
        # plt.fill(*zip(*polygon), alpha=0.4)
    p+= 1

for pol in polyind:
    for i in range(len(pol)):
        e = {pol[i], pol[np.mod(i+1, len(pol))]} 
        if e not in edgeind:
            edgeind.append(e)

ledgeind = [list(i) for i in edgeind]
for e in ledgeind:
    plt.plot([goodvertices[e[0]][0], goodvertices[e[1]][0]], [goodvertices[e[0]][1], goodvertices[e[1]][1]], 'k-')
plt.xlim(min_x-.01, max_x+.01)
plt.ylim(min_y-.01, max_y+.01)
plt.axis('off')
plt.axes().set_aspect('equal', 'datalim')
plt.show()
plt.close()
            





# #Let's put things back in language for normal Voronoi diagrams

polys = [i for i in polyind]
edges = [i for i in edgeind]



  







#Gotta have my pops...

pops = 1
sampno = 1
while pops < numpops:
    outfaceind = []
    inface = []
#    print(pops)
    #first we have to pick a random edge
    noside = 0
    edgewalls = 0
    while noside == 0:
        edge = list(np.random.choice(edges))
        #Check to see if vertex isn't on boundary (should be pretty cheap move)
        if ((0 not in vertices[edge[0]]) and (1 not in vertices[edge[0]])) or ((0 not in vertices[edge[1]]) and (1 not in vertices[edge[1]])):
            noside = 1
            if (0 in vertices[edge[0]]) or (0 in vertices[edge[1]]):
                edgewalls +=1 
            if (1 in vertices[edge[0]]) or (1 in vertices[edge[1]]):
                edgewalls += 1
#    print(edge)
    Q = []
    #tags mean 0 for vertex neighbor, and 1 for edge neighbor
    Qtag = []
    esum = 0
    vmin = 10
    QEind = []
    QVind = []
    #What faces contain edge vertices?
    for i in range(len(polys)):
        if (edge[0] in polys[i]) and (edge[1] in polys[i]):
            outfaceind.append(i)
            Q.append(polys[i])
            Qtag.append(1)
            esum += len(polys[i])
            QEind.append(i)
        if (edge[0] in polys[i]) ^ (edge[1] in polys[i]):
            outfaceind.append(i)
            Q.append(polys[i])
            Qtag.append(0)
            vmin = min(vmin, len(polys[i]))
            QVind.append(i)

    QEind.reverse()        



#We can check for whether this edge can pop

    if ((len(Q) == 4) or ((len(Q) == 3) and edgewalls == 1) or ((len(Q) == 3) and edgewalls == 2)) and vmin>3 and esum>7:

        
        for j in QVind:
            if edge[0] in polys[j]:
                polys[j].remove(edge[0])
            else:
                polys[j].remove(edge[1])
        
        #Vertex shedding from polys
        


        #Remove edge neighbors from poly list        
        for ind in QEind:
            del polys[ind]

        
        #Time to add the big grain if there are two edge neighbors            
        qq = np.where(np.array(Qtag) == 1)[0]
        A = np.array(Q[qq[0]])
        B = np.array(Q[qq[1]])
        ll = int(np.where(A == edge[0])[0])
        if A[(ll+1)%len(A)] == edge[1]:
            arc1 = [A[(k+2+ll)%len(A)] for k in range(len(A)-2)]
        else:
            arc1 = [A[(ll-2-k)%len(A)] for k in range(len(A)-2)]
            
            
        ll = int(np.where(B == edge[1])[0])
        if B[(ll+1)%len(B)] == edge[0]:
            arc2 = [B[(k+2+ll)%len(B)] for k in range(len(B)-2)]
        else:
            arc2 = [B[(ll-2-k)%len(B)] for k in range(len(B)-2)]
    
        
        polys.append(arc1+arc2)

        

    
        uneigh = []
        vneigh = []
        #This loop finds u and v neighbors, and also deletes all edges containing u or v
        for k in reversed(range(len(edges))):
            if edge[0] in edges[k]:
                if edges[k] != set(edge):
#                    print(edges[k])
                    uneigh.append(list(edges[k].difference({edge[0]}))[0])
                del edges[k]
                continue
        
            if edge[1] in edges[k]:
                if edges[k] != set(edge):
#                    print(edges[k])
                    vneigh.append(list(edges[k].difference({edge[1]}))[0])
                del edges[k]
        edges.append(set(uneigh))
        edges.append(set(vneigh))
        
        


            
        pops += 1

            


#If you want a figure
ledgeind = [list(i) for i in edges]
for e in ledgeind:
    plt.plot([goodvertices[e[0]][0], goodvertices[e[1]][0]], [goodvertices[e[0]][1], goodvertices[e[1]][1]], 'k-')
plt.xlim(min_x-.01, max_x+.01)
plt.ylim(min_y-.01, max_y+.01) 
plt.axis('off')
plt.axes().set_aspect('equal', 'datalim')
plt.show()



