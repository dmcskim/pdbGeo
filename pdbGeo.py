__author__ = 'dim'
#!/usr/bin/env python
from glob import glob
from itertools import combinations_with_replacement
from MDAnalysis import Universe
from MDAnalysis.analysis.align import alignto
from numpy import dot, cross, arange, array
from numpy.linalg import norm
from math import ceil
import csv
import pickle


def getnormal(pqr):
    try:
        u = Universe(pqr)
    except:
        return None
    return __getnormal(u)


def __getnormal(univ, atoms=(58, 74, 133)):
    if len(atoms) != 3:
        print 'error in __getNormal'
        return None
    A = univ.selectAtoms('name CA and resid %d' % (atoms[0],))
    for atom in atoms[1:]:
        A += univ.selectAtoms('name CA and resid %d' % (atom,))
    mat = array([item.pos for item in A])
    if len(mat) != 3:
        return None
    a = mat[0, :]-mat[1, :]
    b = mat[1, :]-mat[2, :]
    axb = cross(a, b)
    return axb/norm(axb)


def alignstructures(pdbs, ref=None, resids=(58, 74, 133, 184, 220)):
    sele = []
    if ref == None:
        ref = pdbs.pop()
    for item in resids:
        sele.append('resid '+str(item))
    ref0 = Universe(ref, format='pdb')
    temp = ' or '.join(sele)
    sele = 'protein and name CA and ('+temp+')'
    for item in pdbs:
        if item != ref:
            mob = Universe(item, format='pdb')
            alignto(mob, ref0, select=sele, mass_weighted=True)
            mob.atoms.write(item[:-4]+'.aln.pdb')
    return


def getregion(pqr, start, stop):
    try:
        u = Universe(pqr, format='pqr')
    except:
        try:
            u = Universe(pqr, format='pdb')
        except:
            return None
    return __getregion(u, start, stop)


def __getregion(u, start, stop):
    A = u.selectAtoms('name CA and resid %d' % (start,))
    start += 1
    while start <= stop:
        A += u.selectAtoms('name CA and resid %d' % (start,))
        start += 1
    return A


def getvector(atoms):
    return getangle(atoms)


def getangle(atoms):
    angles = []
    for i in range(len(atoms)-1):
        temp = atoms[i+1].position+atoms[i].position
        if temp[0] == 'nan':
            return 0
        angles.append(temp/norm(temp))
    if len(angles) > 0:
        ret = sum(angles)
        return ret/norm(ret)
    return None


def getrevpdbs():
    return pickle.load(open('revpdbs.p', 'r'))
    #return pickle.load(open('/share/Data/aligned_PDBs/revpdbs.p','r'))


def angles(univ, plane, start, end):
    try:
        u = Universe(univ, format='pdb')
    except:
        try:
            u = Universe(univ, format='pqr')
        except:
            return None, None
    norm = __getnormal(u, atoms=plane)
    temp = __getregion(u, float(start), float(end))
    if norm is not None and temp is not None and len(temp) > 0:
        angle = getangle(temp)
        name = univ.split('/')[-1]
        if angle is not None:
            return name, dot(norm, angle)
    return None, None


def traceprot(u, plane, start, end, window):
    start = int(start)
    norm = __getnormal(u, atoms=plane)
    trace = []
    while start < int(end)-int(window):
        temp = __getregion(u, float(start), float(start)+float(window))
        if norm is not None and temp is not None and len(temp) > 0:
                angle = getangle(temp)
                if angle is not None:
                    trace.append((int(start), float(dot(norm, angle))))
                start += 1
    return trace


def getcolors(num):
    #returns a list containing (r,g,b) combinations for use in matplotlib
    #where each r,g,b can take one of num values
    temp = range(1/(2*num), num)/norm(num)
    return combinations_with_replacement(temp, 3)


def process_args(args):
    revpdb = getrevpdbs()
    start, end = args.region.split('-')
    plane = [int(x) for x in args.plane.split(',')]
    if args.testing:
        files = ['pdb2gs2.cnt_A.okd', 'pdb2itu.cnt_A.pdb']
        for item in files:
            reg = getregion(item, 40, 50)
            print item, reg
            print getvector(reg)
            #files = glob('*.okd')
            #alignStructures(files)
        return

    if args.t:
        import matplotlib.pyplot as plt
        files = []
        for kin in args.kinase:
            for item in revpdb[kin]:
                files += glob('okdevs/*'+item.lower()+'*.okd')
                #files += glob('/share/Data/aligned_PDBs/*'+item.lower()+'*.okd')
        ncolors = len(files)
        ndiv = int(ceil(pow(ncolors, 0.33)))
        colors = getcolors(ndiv)
        for c in colors:
            if len(files) > 0:
                item = files.pop()
                print item
            else:
                break
        try:
            u = Universe(item, format='pqr')
        except:
            try:
                u = Universe(item, format='pdb')
            except:
                print '\t'+item+' is not a properly formatted pdb or pqr file.'
        value = traceprot(u, plane, start, end, 5)
        x, y = zip(*value)
        plt.plot(x, y, color=c)
        plt.xticks(range(1, 240, 5))
        plt.show()
        return
    if args.classes:
        classes = getclasses()
        for item in classes.keys():
            print 'hi'
            # check and graph each class
    elif args.check:
        for kin in args.kinase:
            print kin, '\t', sorted(revpdb[kin])
    else:
        results, files, x, y = [], [], [], []
        if 'all' in args.kinase:
            args.kinase = revpdb.keys()
        for kin in args.kinase:
            for item in revpdb[kin]:
                files += glob('okdevs/*'+item.lower()+'*.okd')
                #files += glob('/share/Data/aligned_PDBs/*'+item.lower()+'*.okd')
        for item in files:
            name, value = angles(item, plane, start, end)
            if value != None:
                print name, '\t', value
            if args.graph:
                x.append(name)
                y.append(value)
            if args.graph and len(y) > 0:
                import matplotlib.pyplot as plt
                ysort = sorted(zip(x, y), key=lambda kv: float(kv[1]))
                x, y = zip(*ysort)
                ind = arange(len(x))
                plt.plot(y, ind, color='r')
                plt.yticks(ind, x)
                plt.show()
    return


def getclasses():
    wanted = {}
    reader = csv.reader(open('/share/Data/aligned_PDBs/class.txt', 'r'), delimiter='\t')
    for row in reader:
        if row[1] not in wanted:
            wanted[row[1]] = []
        wanted[row[1]].append(row[0])
    pickle.dump(wanted, open('class.p', 'w'))
    return wanted

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, epilog="""Examples, how to use, etc
    Typical usage:  ./pdbGeo EGFR
            ./pdbGeo -g EGFR
            ./pdbGeo -r 1-100 EGFR
            ./pdbGeo -p 2,3,4 -r 20-50 -g EGFR
    """,formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('kinase', metavar='K', type=str, nargs='+', help='name of kinase to evaluate')
    ap.add_argument('-c', '--check', action='store_true', help='check for pdb availability')
    ap.add_argument('-g', '--graph', action='store_true', help='sort and graph results')
    #default region - C-helix
    ap.add_argument('-r', '--region', type=str, action='store', default='39-53', help='region of interest')
    ap.add_argument('-p', '--plane', type=str, action='store', default='58,74,133', help='region of interest')
    ap.add_argument('--testing', action='store_true', help='testing')
    ap.add_argument('--classes', action='store_true', help='testing')
    ap.add_argument('-t', action='store_true', help='trace protein')
    args = ap.parse_args()
    process_args(args)


