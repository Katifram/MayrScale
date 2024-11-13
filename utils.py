import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import cm

from skimage import measure

def parse_cube(filename):
    #from: https://github.com/psi4/psi4numpy/blob/6ed03e715689ec82bf96fbb23c1855fbe7835b90/Tutorials/14_Visualization/vizualize.ipynb
    """ Parses a cube file, returning a dict of the information contained.
        The cubefile itself is stored in a numpy array. """
    with open(filename) as fp:
        results = {}

        # skip over the title
        fp.readline()
        fp.readline()

        origin = fp.readline().split()
        natoms = int(origin[0])
        results['minx'] = minx = float(origin[1])
        results['miny'] = miny = float(origin[2])
        results['minz'] = minz = float(origin[3])

        infox = fp.readline().split()
        numx = int(infox[0])
        incx = float(infox[1])
        results['incx'] = incx
        results['numx'] = numx
        results['maxx'] = minx + incx * numx

        infoy = fp.readline().split()
        numy = int(infoy[0])
        incy = float(infoy[2])
        results['incy'] = incy
        results['numy'] = numy
        results['maxy'] = miny + incy * numy

        infoz = fp.readline().split()
        numz = int(infoz[0])
        incz = float(infoz[3])
        results['incz'] = incz
        results['numz'] = numz
        results['maxz'] = minz + incz * numz

        atnums = []
        coords = []
        for atom in range(natoms):
            coordinfo = fp.readline().split()
            atnums.append(int(coordinfo[0]))
            coords.append(list(map(float, coordinfo[2:])))
        results['atom_numbers'] = np.array(atnums)
        results['atom_coords'] = np.array(coords)

        data = np.array([ float(entry) for line in fp for entry in line.split() ])
        if len(data) != numx*numy*numz:
            raise Exception("Amount of parsed data is inconsistent with header in Cube file!")
        results['data'] = data.reshape((numx,numy,numz))

        return results
    

def draw_isosurface(parsed_cube, iso_value):

    vert, faces, norm, values= measure.marching_cubes(parsed_cube['data'], iso_value, spacing=(parsed_cube['incx'],parsed_cube['incy'],parsed_cube['incz']))

    print(parsed_cube['data'])

        
    # Set up a 3D plot
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Create the 3D mesh from vertices and faces
    mesh = Poly3DCollection(vert[faces], alpha=0.3, edgecolor='k', linewidth=0.1)
    mesh.set_facecolor([0.5, 0.5, 1, 0.9])  # Set color and transparency

    ax.add_collection3d(mesh)

    # Plot atoms with correct positions relative to the origin
    atom_coords = parsed_cube['atom_coords'] - np.array([parsed_cube['minx'], parsed_cube['miny'], parsed_cube['minz']])
    atom_numbers = parsed_cube['atom_numbers']
    colors = cm.get_cmap("viridis", len(set(atom_numbers)))  # Color map for atom types

    

    # Plot atoms and labels
    for i, (atom_num, coord) in enumerate(zip(atom_numbers, atom_coords)):
        ax.scatter(*coord, color=colors(i), s=100, label=f'Atom {atom_num}')
        ax.text(*coord, f'{atom_num}', color='black', fontsize=12, ha='center')  # Label with atomic number

    # Set plot limits based on the shape of the cube data
    ax.set_xlim(0, parsed_cube['data'].shape[0] * parsed_cube['incx'])
    ax.set_ylim(0, parsed_cube['data'].shape[1] * parsed_cube['incy'])
    ax.set_zlim(0, parsed_cube['data'].shape[2] * parsed_cube['incz'])

    

    # Set axis labels
    ax.set_xlabel("X / bohr")
    ax.set_ylabel("Y / bohr")
    ax.set_zlabel("Z / bohr")

    # Show the plot
    plt.show()