#!/usr/bin/env python
from MDAnalysis import Universe
from numpy import dot, cross, arange, array
from numpy.linalg import norm

def get_normal(univ, atoms=(58, 74, 133)):
    if len(atoms) != 3:
        print 'error in get_normal'
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

def get_universe(infile):
	try: 
		u = Universe(infile)
		return u
	except:
		print '%s is not a recognied file type.'% (infile,)
	return None

def get_region(u, start, stop, atom='CA'):
	rlist = ['name %s and resid %d' % (atom,s) for s in range(int(start),int(stop))]
	A = u.selectAtoms(*rlist)
	return A

def get_vector(atoms):
    result = []
    for i in range(len(atoms)-1):
        temp = atoms[i+1].position+atoms[i].position
        if temp[0] == 'nan':
            return 0
        result.append(temp/norm(temp))
    if len(result) > 0:
        ret = sum(result)
        return ret/norm(ret)
    return None

def angles(univ, plane, start, end):
	""" Returns angle of <start, end> to plane """
	u = get_universe(univ)
	norm = get_normal(u, atoms=plane)
	temp = get_region(u, start, end)
	if norm is not None and temp is not None and len(temp) > 0:
		angle = get_vector(temp)
		name = univ.split('/')[-1]
		if angle is not None:
			return name, dot(norm, angle)
	return None, None

def process_trajectory(args):
	start, end = args.region.split('-')
	plane = [int(x) for x in args.plane.split(',')]
	univ = Universe(args.topo,args.traj)
	results = []
	for ts in univ.trajectory:
		norm = get_normal(univ, atoms=plane)
		temp = get_region(univ, start, end)
		if norm is not None and temp is not None and len(temp) > 0:
			angle = get_vector(temp)
			name = univ.split('/')[-1]
			if angle is not None:
				results.append(dot(norm, angle))
	if args.graph and len(results) > 0:
		import matplotlib.pyplot as plt
		ysort = sorted(list(enumerate(results,start=1)), key=lambda kv: float(kv[1]))
		x, y = zip(*ysort)
		ind = arange(len(x))
		plt.plot(y, ind, color='r')
		plt.yticks(ind, x)
		print ysort
		plt.show()
	return 

def process_files(args):
	start, end = args.region.split('-')
	plane = [int(x) for x in args.plane.split(',')]
	results, files, x, y = [], [], [], []

	if len(args.infiles) > 1:
		for item in args.infiles:
			name, value = angles(item, plane, start, end)
			if value != None:
				print name, '\t', value
			if args.graph:
				x.append(name)
				y.append(value)
	else:
		name, value = angles(args.infiles[0], plane, start, end)
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

def process_args(args):
	if args.topo and args.traj:
		process_trajectory(args)
	else:
		process_files(args)
	return

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, epilog="""Examples, how to use, etc
    Typical usage:  
            ./pdbGeo <infiles>
            ./pdbGeo -r 1-100 <infiles>
            ./pdbGeo -p 2,3,4 -r 20-50 -g <infiles>
    """,formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-i', '--infiles', type=str, nargs='+', help='files to evaluate')
    ap.add_argument('-g', '--graph', action='store_true', help='sort and graph results')

    #default region - C-helix
    ap.add_argument('-r', '--region', type=str, action='store', default='39-53', help='region of interest')
    ap.add_argument('-p', '--plane', type=str, action='store', default='58,74,133', help='region of interest')
    ap.add_argument('--topo', type=str, action='store', help='MD topology file.')
    ap.add_argument('--traj', type=str, action='store', help='MD trajectory file.')
    args = ap.parse_args()
    process_args(args)


