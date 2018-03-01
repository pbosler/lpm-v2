import numpy as np
import matplotlib.pyplot as plt

def bilinearPlaneMap( v0, v1, v2, v3, r1, r2 ):
    result = np.zeros(2)
    result[0] = 0.25 * ( (1.0-r1)*(1.0+r2)*v0[0] + (1.0-r1)*(1.0-r2)*v1[0] + 
        (1.0+r1)*(1.0-r2)*v2[0] + (1.0+r1)*(1.0+r2)*v3[0])
    result[1] = 0.25 * ( (1.0-r1)*(1.0+r2)*v0[1] + (1.0-r1)*(1.0-r2)*v1[1] + 
        (1.0+r1)*(1.0-r2)*v2[1] + (1.0+r1)*(1.0+r2)*v3[1])
    return result
    


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

def edgeXy(xyz, orig, dest, ints):
    result = np.zeros([4,2])
    result[0] = xyz[orig]
    result[1] = xyz[ints[0]]
    result[2] = xyz[ints[1]]
    result[3] = xyz[dest]
    return result
    
def writeSeedFile():
    xyz, origs, dests, lefts, rights, ints, faceVerts, faceCenters, faceEdges = quadCubicSeed()
    with open("quadCubicSeed.dat", 'w') as f:
        f.write("x     y \n")
        for x, y in xyz:
            f.write(str(x) + "    " + str(y) + "\n")
        f.write("edgeO      edgeD       edgeLeft        edgeRight    edgeInts \n")
        for i in range(12):
            f.write(str(origs[i]) + "   " + str(dests[i]) + "   " + str(lefts[i]) + "   " + str(rights[i]) + "   "
                + str(ints[i,0]) + "   " + str(ints[i,1]) + "\n")
        f.write("faceverts\n")
        for v in faceVerts:
            f.write(str(v)[1:-1])
            f.write('\n')
        f.write('faceedges\n')
        for e in faceEdges:
            f.write(str(e[0]) + " " + str(e[1]) + " " + str(e[2]) + " " + str(e[3])+"\n")
        f.write('facecenters\n')
        for c in faceCenters:
            f.write(str(c)[1:-1] + "\n")
    
if (__name__ == "__main__"):
    print "quad maps"
    
    m_size = 4.0
    m_width = m_size / 4.0
    l_width= 2.0
    
    xyz, origs, dests, lefts, rights, ints, faceVerts, faceCenters, faceEdges = quadCubicSeed()
    
    fig0, (ax0, ax1) = plt.subplots(1,2)
    plt.tight_layout(pad=0.1, w_pad=6)
    ax0.plot(xyz[:,0], xyz[:,1], 'ko', markersize=m_size)
    ax0.set(title='edges & particles', xlabel='x', ylabel='y')
    ax0.set_aspect('equal','box')
    
    for i in range(49):
        ax0.text(xyz[i,0], xyz[i,1], str(i+1), color='k')
    
    for i in range(12):
        exy = edgeXy(xyz, origs[i], dests[i], ints[i])
        dx = exy[1:,0] - exy[0:-1,0]
        dy = exy[1:,1] - exy[0:-1,1]
        midpt = 0.5 * (exy[0] + exy[3])
#         ax0.arrow(exy[0,0], exy[0,1], dx[0], dy[0],head_width=0.1, head_length=0.05, fc='r', ec='r', overhang=0.2)
#         ax0.arrow(exy[1,0], exy[1,1], dx[1], dy[1],head_width=0.1, head_length=0.05, fc='r', ec='r')
#         ax0.arrow(exy[2,0], exy[2,1], dx[2], dy[2],head_width=0.1, head_length=0.05, fc='r', ec='r')
        ax0.arrow(exy[0,0], exy[0,1], midpt[0]-exy[0,0], midpt[1]-exy[0,1], head_width=0.1, 
            head_length=0.05, fc='r', ec='r', length_includes_head=False)
        ax0.plot([midpt[0],exy[3,0]],[midpt[1],exy[3,1]],'r-')
        ax0.text(midpt[0], midpt[1]+0.01, str(i), color='r')
    
    ax1.set_aspect('equal','box')
    ax1.set(title='faces & edges',xlabel='x',ylabel='y')
    for i in range(12):
        exy = edgeXy(xyz, origs[i], dests[i], ints[i])
        dx = exy[1:,0] - exy[0:-1,0]
        dy = exy[1:,1] - exy[0:-1,1]
        midpt = 0.5 * (exy[0] + exy[3])
        ax1.arrow(exy[0,0], exy[0,1], midpt[0]-exy[0,0], midpt[1]-exy[0,1], head_width=0.1, 
            head_length=0.05, fc='r', ec='r', length_includes_head=False)
        ax1.plot([midpt[0],exy[3,0]],[midpt[1],exy[3,1]],'r-')
        ax1.text(midpt[0], midpt[1]+0.01, str(i), color='r')
    
    for i in range(4):
        cntd = np.zeros(2)
        for j in range(12):
            cntd += xyz[faceVerts[i,j]]
        for j in range(4):
            cntd += xyz[faceCenters[i,j]]
        cntd /= 16.0
        ax1.text(cntd[0], cntd[1], str(i), color='b', bbox=dict(facecolor='b', alpha=0.25))
    
    fig0.savefig('quadCubicPlaneSeed.pdf', bbox_inches='tight')
    plt.close(fig0)
    
    writeSeedFile()