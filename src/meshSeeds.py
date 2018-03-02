import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def bilinearPlaneMap( v0, v1, v2, v3, r1, r2 ):
    result = np.zeros(2)
    result[0] = 0.25 * ( (1.0-r1)*(1.0+r2)*v0[0] + (1.0-r1)*(1.0-r2)*v1[0] + 
        (1.0+r1)*(1.0-r2)*v2[0] + (1.0+r1)*(1.0+r2)*v3[0])
    result[1] = 0.25 * ( (1.0-r1)*(1.0+r2)*v0[1] + (1.0-r1)*(1.0-r2)*v1[1] + 
        (1.0+r1)*(1.0-r2)*v2[1] + (1.0+r1)*(1.0+r2)*v3[1])
    return result

def cubedSphereSeed():
    oor3 = 1.0 / np.sqrt(3.0)
    xyz = np.zeros([14,3])
    xyz[0] = (oor3, -oor3, oor3)
    xyz[1] = (oor3, -oor3, -oor3)
    xyz[2] = (oor3, oor3, -oor3)
    xyz[3] = (oor3, oor3, oor3)
    xyz[4] = (-oor3, oor3, -oor3)
    xyz[5] = (-oor3, oor3, oor3)
    xyz[6] = (-oor3, -oor3, -oor3)
    xyz[7] = (-oor3, -oor3, oor3)
    xyz[8] = (1.0, 0.0, 0.0)
    xyz[9] = (0.0, 1.0, 0.0)
    xyz[10] = (-1.0, 0.0, 0.0)
    xyz[11] = (0.0, -1.0, 0.0)
    xyz[12] = (0.0, 0.0, 1.0)
    xyz[13] = (0.0, 0.0, -1.0)
    edgeOrigs = np.array([0,1,2,3,2,4,5,4,6,7,6,0],dtype=int)
    edgeDests = np.array([1,2,3,0,4,5,3,6,7,5,1,7],dtype=int)
    edgeLefts = np.array([0,0,0,0,1,1,1,2,2,2,3,3],dtype=int)
    edgeRights= np.array([3,6,1,4,5,2,4,5,3,4,5,4],dtype=int)
    edgeInteriors = None
    faceVerts = np.array([[0,1,2,3], [3,2,4,5], [5,4,6,7], [7,6,1,0], [7,0,3,5],[1,6,4,2]],dtype=int)
    faceEdges = np.array([[0,1,2,3], [2,4,5,6], [5,7,8,9], [8,10,0,11], [11,3,6,9],[10,7,4,1]],dtype=int)
    faceCenters = np.array([8,9,10,11,12,13],dtype=int)
    return xyz, edgeOrigs, edgeDests, edgeLefts, edgeRights, edgeInteriors, \
        faceVerts, faceCenters, faceEdges

def sphereTriCenter(v0, v1, v2):
    result = (v0 + v1 + v2)/3.0
    norm = np.sqrt(np.sum(np.square(result)))
    result /= norm
    return result
        
def icosTriSeed():
    xyz = np.zeros([32,3])
    xyz[0] = (0.0,0.0,1.0)
    xyz[1] = (0.723606797749978969640917366873,0.525731112119133606025669084848,0.447213595499957939281834733746)
    xyz[2] = (-0.276393202250021030359082633126,0.850650808352039932181540497063,0.447213595499957939281834733746)
    xyz[3] = (-0.894427190999915878563669467492,0.0,0.447213595499957939281834733746,)
    xyz[4] = (-0.276393202250021030359082633127,-0.850650808352039932181540497063,0.447213595499957939281834733746,)
    xyz[5] = (0.723606797749978969640917366873,-0.525731112119133606025669084848,0.447213595499957939281834733746,)
    xyz[6] = (0.894427190999915878563669467492,0.0,-0.447213595499957939281834733746,)
    xyz[7] = (0.276393202250021030359082633127,0.850650808352039932181540497063,-0.447213595499957939281834733746,)
    xyz[8] = (-0.723606797749978969640917366873,0.525731112119133606025669084848,-0.447213595499957939281834733746,)
    xyz[9] = (-0.723606797749978969640917366873,-0.525731112119133606025669084848,-0.447213595499957939281834733746,)
    xyz[10] = (0.276393202250021030359082633127,-0.850650808352039932181540497063,-0.447213595499957939281834733746,)
    xyz[11] = (0.0,0.0,-1.0,)
    edgeOrigs = np.array([0,1,2,2,0,3,4,4,5,5,1,6,7,7,7,8,8,8,3,9,9,10,10,10,6,11,11,8,11,10],dtype=int)
    edgeDests = np.array([1,2,0,3,3,4,0,5,0,1,6,7,1,2,8,2,3,9,9,4,10,4,5,6,5,6,7,11,9,11],dtype=int)
    edgeLefts = np.array([0,0,0,1,2,2,2,3,3,4,5,5,5,6,7,7,8,9,10,10,11,11,12,13,13,19,15,17,17,19],dtype=int)
    edgeRights= np.array([4,6,1,8,1,10,3,12,4,14,14,15,6,7,16,8,9,17,9,11,18,12,13,19,14,15,16,16,18,18],dtype=int)
    edgeInteriors = None
    faceVerts = np.zeros([20,3], dtype=int)
    faceVerts[0] = (0,1,2)
    faceVerts[1] = (0,2,3)
    faceVerts[2] = (0,3,4)
    faceVerts[3] = (0,4,5)
    faceVerts[4] = (0,5,1)
    faceVerts[5] = (1,6,7)
    faceVerts[6] = (7,2,1)
    faceVerts[7] = (2,7,8)
    faceVerts[8] = (8,3,2)
    faceVerts[9] = (3,8,9)
    faceVerts[10] = (9,4,3)
    faceVerts[11] = (4,9,10)
    faceVerts[12] = (10,5,4)
    faceVerts[13] = (5,10,6)
    faceVerts[14] = (6,1,5)
    faceVerts[15] = (11,7,6)
    faceVerts[16] = (11,8,7)
    faceVerts[17] = (11,9,8)
    faceVerts[18] = (11,10,9)
    faceVerts[19] = (11,6,10)
    faceEdges = np.zeros([20,3],dtype=int)
    faceEdges[0] = (0,1,2)
    faceEdges[1] = (2,3,4)
    faceEdges[2] = (4,5,6)
    faceEdges[3] = (6,7,8)
    faceEdges[4] = (8,9,0)
    faceEdges[5] = (10,11,12)
    faceEdges[6] = (13,1,12)
    faceEdges[7] = (13,14,15)
    faceEdges[8] = (16,3,15)
    faceEdges[9] = (16,17,18)
    faceEdges[10] = (19,5,18)
    faceEdges[11] = (19,20,21)
    faceEdges[12] = (22,7,21)
    faceEdges[13] = (22,23,24)
    faceEdges[14] = (10,9,24)
    faceEdges[15] = (26,11,25)
    faceEdges[16] = (27,14,26)
    faceEdges[17] = (28,17,27)
    faceEdges[18] = (29,20,28)
    faceEdges[19] = (25,23,29)
    faceCenters = range(12,32)
    for i in range(20):
        xyz[i+12] = sphereTriCenter(xyz[faceVerts[i][0]], xyz[faceVerts[i][1]], xyz[faceVerts[i][2]])
    return xyz, edgeOrigs, edgeDests, edgeLefts, edgeRights, edgeInteriors, \
        faceVerts, faceCenters, faceEdges
 
def triHexSeed():
    pio3 = np.pi/3.0
    xyz = np.zeros([13,2])
    xyz[0] = (0.0, 0.0)
    xyz[1] = (np.cos(pio3), np.sin(pio3))
    xyz[2] = (np.cos(2*pio3), np.sin(2*pio3))
    xyz[3] = (-1.0, 0.0)
    xyz[4] = (np.cos(4*pio3), np.sin(4*pio3))
    xyz[5] = (np.cos(5*pio3), np.sin(5*pio3))
    xyz[6] = (1.0, 0.0)
    xyz[7] = (xyz[0] + xyz[1] + xyz[2]) / 3.0
    xyz[8] = (xyz[0] + xyz[2] + xyz[3]) / 3.0
    xyz[9] = (xyz[0] + xyz[3] + xyz[4]) / 3.0
    xyz[10]= (xyz[0] + xyz[4] + xyz[5]) / 3.0
    xyz[11] = (xyz[0] + xyz[5] + xyz[6])/ 3.0
    xyz[12] = (xyz[0] + xyz[6] + xyz[1])/ 3.0
    edgeOrigs = np.array([0,1,2,3,4,5,6,2,3,4,0,0],dtype=int)
    edgeDests = np.array([1,2,3,4,5,6,1,0,0,0,5,6],dtype=int)
    edgeLefts = np.array([0,0,1,2,3,4,5,0,1,2,4,5],dtype=int)
    edgeRights= np.array([5,-1,-1,-1,-1,-1,-1,1,2,3,3,4],dtype=int)
    edgeInteriors = None
    faceVerts = np.array([[0,1,2], [2,3,0],[4,0,3], [0,4,5], [5,6,0], [1,0,6]],dtype=int)
    faceCenters = np.array([7,8,9,10,11,12],dtype=int)
    faceEdges = np.array([[0,1,7],[2,8,7],[9,8,3],[9,4,10],[5,11,10],[0,11,6]],dtype=int)
    return xyz, edgeOrigs, edgeDests, edgeLefts, edgeRights, edgeInteriors, \
        faceVerts, faceCenters, faceEdges

def quadCubicSeed():
    sqrt5 = np.sqrt(5.0)
    qps = np.array([-1.0, -1.0/sqrt5, 1.0/sqrt5, 1.0])
    qws = np.array([1.0/6.0, 5.0/6.0, 5.0/6.0, 1.0/6.0])
    f0verts = np.array(range(12),dtype=int)
    f1verts = np.array([3,16,17,18,19,20,21,22,23,6,5,4],dtype=int)
    f2verts = np.array([6,23,22,21,28,29,30,31,32,33,34,35],dtype=int)
    f3verts = np.array([9,8,7,6,35,34,33,40,41,42,43,44],dtype=int)
    
    faceVerts = np.array([f0verts, f1verts, f2verts, f3verts], dtype=int)
    
    f0centers = np.array([12,13,14,15],dtype=int)
    f1centers = np.array([24,25,26,27],dtype=int)
    f2centers = np.array([36,37,38,39],dtype=int)
    f3centers = np.array([45,46,47,48],dtype=int)
    
    faceCenters = np.array([f0centers, f1centers, f2centers, f3centers], dtype=int)
    
    faceEdges = np.array([[0,8,11,7], [1,2,10,8], [10,3,4,9], [11,9,5,6]], dtype=int)
    
    edgeOrigs = np.array([0,3,18,21,30,33,42,9,3,6,21,6],dtype=int)
    edgeDests = np.array([3,18,21,30,33,42,9,0,6,33,6,9],dtype=int)
    edgeInteriors = np.array([[1,2],[16,17], [19,20], [28,29], [31,32], [40,41], [43,44], [10,11],
        [4,5],[35,34],[22,23],[7,8]])
    edgeLefts = np.array([0,1,1,2,2,3,3,0,0,3,1,2],dtype=int)
    edgeRights = np.array([-1,-1,-1,-1,-1,-1,-1,-1,1,2,2,3],dtype=int)
    xyz = np.zeros([49,2])
    xyz[0] = (-1.0, 1.0)
    xyz[3] = (-1.0, 0.0)
    xyz[18] = (-1.0,-1.0)
    xyz[21] = (0.0, -1.0)
    xyz[30] = (1.0, -1.0)
    xyz[33] = (1.0, 0.0)
    xyz[42] = (1.0, 1.0)
    xyz[9] = (0.0, 1.0)
    xyz[6] = (0.0,0.0)
    
    xyz[1] = bilinearPlaneMap(xyz[0], xyz[3], xyz[6], xyz[9], -1.0, 1.0/sqrt5)
    xyz[2] = bilinearPlaneMap(xyz[0], xyz[3], xyz[6], xyz[9], -1.0, -1.0/sqrt5)
    xyz[4] = bilinearPlaneMap(xyz[0], xyz[3], xyz[6], xyz[9], -1.0/sqrt5, -1.0)
    xyz[5] = bilinearPlaneMap(xyz[0], xyz[3], xyz[6], xyz[9], 1.0/sqrt5, -1.0)
    xyz[7] = bilinearPlaneMap(xyz[0], xyz[3], xyz[6], xyz[9], 1.0, -1.0/sqrt5)
    xyz[8] = bilinearPlaneMap(xyz[0], xyz[3], xyz[6], xyz[9], 1.0, 1.0/sqrt5)
    xyz[10] = bilinearPlaneMap(xyz[0], xyz[3], xyz[6], xyz[9], 1.0/sqrt5, 1.0)
    xyz[11] = bilinearPlaneMap(xyz[0], xyz[3], xyz[6], xyz[9], -1.0/sqrt5, 1.0)
    xyz[12] = bilinearPlaneMap(xyz[0], xyz[3], xyz[6], xyz[9], -1.0/sqrt5, 1.0/sqrt5)
    xyz[13] = bilinearPlaneMap(xyz[0], xyz[3], xyz[6], xyz[9], -1.0/sqrt5, -1.0/sqrt5)
    xyz[14] = bilinearPlaneMap(xyz[0], xyz[3], xyz[6], xyz[9], 1.0/sqrt5, -1.0/sqrt5)
    xyz[15] = bilinearPlaneMap(xyz[0], xyz[3], xyz[6], xyz[9], 1.0/sqrt5, 1.0/sqrt5)
    
    xyz[16] = bilinearPlaneMap(xyz[3], xyz[18], xyz[21], xyz[6], -1.0, 1.0/sqrt5)
    xyz[17] = bilinearPlaneMap(xyz[3], xyz[18], xyz[21], xyz[6], -1.0, -1.0/sqrt5)
    xyz[19] = bilinearPlaneMap(xyz[3], xyz[18], xyz[21], xyz[6], -1.0/sqrt5, -1.0)
    xyz[20] = bilinearPlaneMap(xyz[3], xyz[18], xyz[21], xyz[6], 1.0/sqrt5, -1.0)
    xyz[22] = bilinearPlaneMap(xyz[3], xyz[18], xyz[21], xyz[6], 1.0, -1.0/sqrt5)
    xyz[23] = bilinearPlaneMap(xyz[3], xyz[18], xyz[21], xyz[6], 1.0, 1.0/sqrt5)
    xyz[24] = bilinearPlaneMap(xyz[3], xyz[18], xyz[21], xyz[6], -1.0/sqrt5, 1.0/sqrt5)
    xyz[25] = bilinearPlaneMap(xyz[3], xyz[18], xyz[21], xyz[6], -1.0/sqrt5, -1.0/sqrt5)
    xyz[26] = bilinearPlaneMap(xyz[3], xyz[18], xyz[21], xyz[6], 1.0/sqrt5, -1.0/sqrt5)
    xyz[27] = bilinearPlaneMap(xyz[3], xyz[18], xyz[21], xyz[6], 1.0/sqrt5, 1.0/sqrt5)
    
    xyz[28] = bilinearPlaneMap(xyz[6], xyz[21], xyz[30], xyz[33], -1.0/sqrt5, -1.0)
    xyz[29] = bilinearPlaneMap(xyz[6], xyz[21], xyz[30], xyz[33], 1.0/sqrt5, -1.0)
    xyz[31] = bilinearPlaneMap(xyz[6], xyz[21], xyz[30], xyz[33], 1.0, -1.0/sqrt5)
    xyz[32] = bilinearPlaneMap(xyz[6], xyz[21], xyz[30], xyz[33], 1.0, 1.0/sqrt5)
    xyz[34] = bilinearPlaneMap(xyz[6], xyz[21], xyz[30], xyz[33], 1.0/sqrt5, 1.0)
    xyz[35] = bilinearPlaneMap(xyz[6], xyz[21], xyz[30], xyz[33], -1.0/sqrt5, 1.0)
    xyz[36] = bilinearPlaneMap(xyz[6], xyz[21], xyz[30], xyz[33], -1.0/sqrt5, 1.0/sqrt5)
    xyz[37] = bilinearPlaneMap(xyz[6], xyz[21], xyz[30], xyz[33], -1.0/sqrt5, -1.0/sqrt5)
    xyz[38] = bilinearPlaneMap(xyz[6], xyz[21], xyz[30], xyz[33], 1.0/sqrt5, -1.0/sqrt5)
    xyz[39] = bilinearPlaneMap(xyz[6], xyz[21], xyz[30], xyz[33], 1.0/sqrt5, 1.0/sqrt5)
    
    xyz[40] = bilinearPlaneMap(xyz[9], xyz[6], xyz[33], xyz[42], 1.0, -1.0/sqrt5)
    xyz[41] = bilinearPlaneMap(xyz[9], xyz[6], xyz[33], xyz[42], 1.0, 1.0/sqrt5)
    xyz[43] = bilinearPlaneMap(xyz[9], xyz[6], xyz[33], xyz[42], 1.0/sqrt5, 1.0)
    xyz[44] = bilinearPlaneMap(xyz[9], xyz[6], xyz[33], xyz[42], -1.0/sqrt5, 1.0)
    xyz[45] = bilinearPlaneMap(xyz[9], xyz[6], xyz[33], xyz[42], -1.0/sqrt5, 1.0/sqrt5)
    xyz[46] = bilinearPlaneMap(xyz[9], xyz[6], xyz[33], xyz[42], -1.0/sqrt5, -1.0/sqrt5)
    xyz[47] = bilinearPlaneMap(xyz[9], xyz[6], xyz[33], xyz[42], 1.0/sqrt5, -1.0/sqrt5)
    xyz[48] = bilinearPlaneMap(xyz[9], xyz[6], xyz[33], xyz[42], 1.0/sqrt5, 1.0/sqrt5)
    
    return xyz, edgeOrigs, edgeDests, edgeLefts, edgeRights, edgeInteriors, \
        faceVerts, faceCenters, faceEdges

def edgeXyz(xyz, orig, dest, ints):
    ptsPerEdge = 2
    if ints is not None:
        ptsPerEdge += len(ints)
    result = np.zeros([ptsPerEdge,2])
    result[0] = xyz[orig]
    if ints is not None:
        for i in range(len(ints)):
            result[i+1] = xyz[ints[i]]
    result[-1] = xyz[dest]
    return result
    
def writeSeedFile(fname, xyz, origs, dests, lefts, rights, ints, faceVerts, faceCenters, faceEdges):
    nparticles = len(xyz)
    nedges = len(origs)
    nfaces = len(faceVerts)
    ncenters = len(faceCenters)
    nverts = len(faceVerts[0])
    nints = 0
    eformat = '%i   %i  %i  %i  '
    if ints is not None:
        nints = len(ints[0])
        for i in range(nints):
            eformat += '%i  '
    eformat += '\n'
    with open(fname, 'w') as f:
        if np.shape(xyz)[1]==2:
            f.write("x     y \n")
            for x, y in xyz:
                f.write('%.17f  %.17f\n'%(x,y))
        else:
            f.write("x   y   z\n")
            for x, y, z in xyz:
                f.write('%.17f  %.17f  %.17f\n'%(x,y,z))
        f.write("edgeO      edgeD       edgeLeft        edgeRight    edgeInts \n")
#         for i in range(nedges):
#             f.write(str(origs[i]) + " " + str(dests[i]) + " " + str(lefts[i]) + " " + str(rights[i]) + " ")
#             if ints is not None:
#                 for cs in ints:
#                     f.write(str(cs)[1:-1] + "\n")
#             else:
#                 f.write("\n")
        if ints is not None:
            for i in range(nedges):
                f.write(eformat%(origs[i], dests[i], lefts[i], rights[i], ints[i,0], ints[i,1]))
        else:
            for i in range(nedges):
                f.write(eformat%(origs[i], dests[i], lefts[i], rights[i]))
        
        f.write("faceverts\n")
        for v in faceVerts:
            f.write(str(v)[1:-1])
            f.write('\n')
        f.write('faceedges\n')
        for e in faceEdges:
            f.write(str(e)[1:-1]+"\n")
        f.write('facecenters\n')
        if ints is None:
            for c in faceCenters:
                f.write(str(c) + '\n')
        else:
            for c in faceCenters:
                f.write(str(c)[1:-1] + "\n")
            
# def plotSphereSeed(oname, xyz, origs, dests, lefts, rights, ints, faceVerts, faceCenters, faceEdges):
#     nparticles = len(xyz)
#     nedges = len(origs)
#     nfaces = len(faceVerts)
#     ncenters = len(faceCenters)
#     nverts = len(faceVerts[0])
#     nints = 0
#     if ints is not None:
#         nints = len(ints[0])
#     m_size = 4.0
#     m_width = m_size / 4.0
#     l_width= 2.0
# #     fig0 = plt.figure()
# #     ax00 = fig0.add_subplot(2,2,1, projection='3d')
# #     ax00.plot(xyz[:,0], xyz[:,1], xyz[:,2], 'ko')
# #     for i in range(nedges):
# #         exyz = sphereEdgeXyz(xyz, origs[i], dests[i])
# #         ax00.plot(exyz[:,0], exyz[:,1], exyz[:,2], 'r-')
# #     ax00 = fig0.add_subplot(2,2,1)
# #     lam = np.linspace(0.0, 2*np.pi, 180, endpoint=False)
# #     the = np.linspace(-0.5*np.pi, 0.5*np.pi, 91)
# #     lon, lat = np.meshgrid(lam,the)
#     crossesPi = np.zeros(nedges, dtype=bool)
#     if 'cubed' in oname:
#         crossesPi[7] = True
#         crossesPi[9] = True
#     elif 'icos' in oname:
#         pass
#     fig0, (ax0, ax1) = plt.subplots(2,1)
#     plt.tight_layout(pad=0.1, w_pad=6)
#     plon = np.arctan2(xyz[:,1], xyz[:,0])
#     plat = np.arctan2(xyz[:,2], np.sqrt(np.sum(np.square(xyz[:,0]) + np.square(xyz[:,1]))))
#     ax0.plot(plon, plat,'ko', markersize=m_size)
#     
#     for i in range(nedges):
#         exyz = sphereEdgeXyz(xyz, origs[i], dests[i], crossesPi[i])
#         elon = np.arctan2(exyz[:,1], exyz[:,0])
#         elat = np.arctan(exyz[:,2] / np.sqrt(np.sum(np.square(exyz[:,0]) + np.square(exyz[:,1]))))
#         if crossesPi[i]:
#             pos = np.where(elon>0.0)
#             neg = np.where(elon<0.0)
#             ax0.plot(elon[pos], elat[pos], 'r-')
#             ax0.plot(elon[neg], elat[neg], 'r-')
#         else:
#             ax0.plot(elon, elat, 'r-')
#     ax0.set(title='edges & particles', xlabel='x', ylabel='y')  
#     #ax0.set_aspect('equal','box') 
# #    ax00.plot_surface(0.99*np.cos(lat)*np.cos(lon), 0.99*np.cos(lat)*np.sin(lon), 0.99*np.sin(lat),
# #        color='white',edgecolor='white',facecolor='white',shade=False)
#     
#     fig0.savefig(oname, bbox_inches='tight')    
#     plt.close(fig0)
# 
# def sphereEdgeXyz(xyz, orig, dest, crossPi=False):
#     ptsPerEdge = 5
#     result = np.zeros([ptsPerEdge,3])
#     tparam = np.linspace(0.0, 1.0, ptsPerEdge)
#     evec = xyz[dest] - xyz[orig]
#     for i,t in enumerate(tparam):
#         result[i,0] = xyz[orig][0] + t * evec[0]
#         result[i,1] = xyz[orig][1] + t * evec[1]
#         result[i,2] = xyz[orig][2] + t * evec[2]
#         norm = np.sqrt(np.sum(np.square(result[i])))
#         result[i] /= norm
#     return result
    
def plotPlaneSeed(oname, xyz, origs, dests, lefts, rights, ints, faceVerts, faceCenters, faceEdges):
    m_size = 4.0
    m_width = m_size / 4.0
    l_width= 2.0
    fig0, (ax0, ax1) = plt.subplots(1,2)
    plt.tight_layout(pad=0.1, w_pad=6)
    ax0.plot(xyz[:,0], xyz[:,1], 'ko', markersize=m_size)
    ax0.set(title='edges & particles', xlabel='x', ylabel='y')
    ax0.set_aspect('equal','box')
    
    nparticles = np.shape(xyz)[0]
    nedges = np.shape(origs)[0]
    nfaces = np.shape(faceVerts)[0]
    ncenters = np.shape(faceCenters)[0]
    nverts = np.shape(faceVerts)[1]
    for i in range(nparticles):
        ax0.text(xyz[i,0], xyz[i,1], str(i), color='k')
    
    for i in range(nedges):
        if ints is not None:
            exy = edgeXyz(xyz, origs[i], dests[i], ints[i])
        else:
            exy = edgeXyz(xyz, origs[i], dests[i], None)
        dx = exy[1:,0] - exy[0:-1,0]
        dy = exy[1:,1] - exy[0:-1,1]
        midpt = 0.5 * (exy[0] + exy[-1])
#         ax0.arrow(exy[0,0], exy[0,1], dx[0], dy[0],head_width=0.1, head_length=0.05, fc='r', ec='r', overhang=0.2)
#         ax0.arrow(exy[1,0], exy[1,1], dx[1], dy[1],head_width=0.1, head_length=0.05, fc='r', ec='r')
#         ax0.arrow(exy[2,0], exy[2,1], dx[2], dy[2],head_width=0.1, head_length=0.05, fc='r', ec='r')
        ax0.arrow(exy[0,0], exy[0,1], midpt[0]-exy[0,0], midpt[1]-exy[0,1], head_width=0.1, 
            head_length=0.05, fc='r', ec='r', length_includes_head=False)
        ax0.plot([midpt[0],exy[-1,0]],[midpt[1],exy[-1,1]],'r-')
        ax0.text(midpt[0], midpt[1]+0.05, str(i), color='r')
    
    ax1.set_aspect('equal','box')
    ax1.set(title='faces & edges',xlabel='x',ylabel='y')
    for i in range(nedges):
        if ints is not None:
            exy = edgeXyz(xyz, origs[i], dests[i], ints[i])
        else:
            exy = edgeXyz(xyz, origs[i], dests[i], None)
        dx = exy[1:,0] - exy[0:-1,0]
        dy = exy[1:,1] - exy[0:-1,1]
        midpt = 0.5 * (exy[0] + exy[-1])
        ax1.arrow(exy[0,0], exy[0,1], midpt[0]-exy[0,0], midpt[1]-exy[0,1], head_width=0.1, 
            head_length=0.05, fc='r', ec='r', length_includes_head=False)
        ax1.plot([midpt[0],exy[-1,0]],[midpt[1],exy[-1,1]],'r-')
        ax1.text(midpt[0], midpt[1]+0.1, str(i), color='r')
    
    for i in range(nfaces):
        cntd = np.zeros(2)
        for j in range(nverts):
            cntd += xyz[faceVerts[i,j]]
        for j in range(ncenters):
            if len(np.shape(faceCenters)) == 2:
                cntd += xyz[faceCenters[i,j]]
            else:
                cntd += xyz[faceCenters[i]]
        cntd /= (nverts+ncenters)
        ax1.text(cntd[0], cntd[1], str(i), color='b', bbox=dict(facecolor='b', alpha=0.25))
    
    fig0.savefig(oname, bbox_inches='tight')
    plt.close(fig0)    
    
if (__name__ == "__main__"):
    
    print "tri hex seed"
    xyz, origs, dests, lefts, rights, ints, faceVerts, faceCenters, faceEdges = triHexSeed()
    writeSeedFile("triHexSeed.dat", xyz, origs, dests, lefts, rights, ints, faceVerts, faceCenters, faceEdges)
    plotPlaneSeed("triHexSeed.pdf", xyz, origs, dests, lefts, rights, ints, faceVerts, faceCenters, faceEdges)
    
    print "quad cubic seed"
    xyz, origs, dests, lefts, rights, ints, faceVerts, faceCenters, faceEdges = quadCubicSeed()
    writeSeedFile("quadCubicSeed.dat", xyz, origs, dests, lefts, rights, ints, faceVerts, faceCenters, faceEdges)
    plotPlaneSeed("quadCubicSeed.pdf", xyz, origs, dests, lefts, rights, ints, faceVerts, faceCenters, faceEdges)
    
    print "cubed sphere seed"
    xyz, origs, dests, lefts, rights, ints, faceVerts, faceCenters, faceEdges = cubedSphereSeed()
    writeSeedFile('cubedSphereSeed.dat', xyz, origs, dests, lefts, rights, ints, faceVerts, faceCenters, faceEdges)
    #plotSphereSeed('cubedSphereSeed.pdf',xyz, origs, dests, lefts, rights, ints, faceVerts, faceCenters, faceEdges)
    
    print "icosahedral triangle sphere seed"
    xyz, origs, dests, lefts, rights, ints, faceVerts, faceCenters, faceEdges = icosTriSeed()
    writeSeedFile('icosTriSphereSeed.dat', xyz, origs, dests, lefts, rights, ints, faceVerts, faceCenters, faceEdges)
    