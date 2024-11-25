import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import cm

from skimage import measure

from scipy.interpolate import RegularGridInterpolator

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


def ray_cube_distance(cube_min, cube_max, point, direction):
    """
    Calculate the distance from a point to the edge of a cube in a given direction.
    
    Args:
    - cube_min: ndarray of shape (3,), minimum coordinates of the cube.
    - cube_max: ndarray of shape (3,), maximum coordinates of the cube.
    - point: ndarray of shape (3,), the starting point.
    - direction: ndarray of shape (3,), the direction vector.
    
    Returns:
    - Distance to the cube edge (float).
    """
    # Normalize the direction vector
    direction = direction / np.linalg.norm(direction)
    
    # Calculate t for each slab (x, y, z planes)
    t_min = (cube_min - point) / direction
    t_max = (cube_max - point) / direction
    
    # Sort t_min and t_max to ensure correct ordering
    t1 = np.minimum(t_min, t_max)
    t2 = np.maximum(t_min, t_max)
    
    # Find the largest t1 and smallest t2
    t_near = np.max(t1)
    t_far = np.min(t2)
    
    # If t_near > t_far or t_far < 0, the ray does not intersect the cube
    if t_near > t_far or t_far < 0:
        return None  # No intersection
    
    # Return the nearest positive t
    return t_near if t_near > 0 else t_far


def values_along_direction(parsed_cube, start, direction, num_points):
    """
    Extracts interpolated values along a line from `start` in `direction` in cube data and computes distances.

    Parameters:
    - cube: dict, parsed cube data with keys 'data', 'minx', 'miny', 'minz', 'incx', 'incy', 'incz'
    - start: list or array of [x, y, z] coordinates for the starting point in Angstroms.
    - direction: list or array of [dx, dy, dz], the direction vector.
    - num_points: Number of points to sample along the line.
    - step_size: Distance between points along the direction vector in Angstroms.

    Returns:
    - line_values: numpy array of interpolated values along the direction line.
    - line_points: numpy array of sampled points along the direction line.
    - distances: numpy array of distances from the starting point.
    """

    # Normalize the direction vector
    direction = np.array(direction) / np.linalg.norm(direction)

    # find distance between points from parsed_cube by
    min_array = [parsed_cube['minx'], parsed_cube['miny'], parsed_cube['minz']]
    max_array = [parsed_cube['maxx'], parsed_cube['maxy'], parsed_cube['maxz']]
    step_size = ray_cube_distance(min_array, max_array, start, direction) / num_points

    # Generate points along the line at each step_size interval
    line_points = np.array([start + i * step_size * direction for i in range(num_points)])

    # Calculate distances from the start point
    distances = np.arange(num_points) * step_size  # Distance = index * step_size

    # Define the grid axes based on cube grid
    x = np.linspace(start=parsed_cube['minx'], stop=parsed_cube['maxx'], num=parsed_cube['numx'])
    y = np.linspace(start=parsed_cube['miny'], stop=parsed_cube['maxy'], num=parsed_cube['numy'])
    z = np.linspace(start=parsed_cube['miny'], stop=parsed_cube['maxy'], num=parsed_cube['numy'])

    # Create an interpolator function for the 3D grid
    interpolator = RegularGridInterpolator((x, y, z), parsed_cube['data'])

    # Interpolate data values along the line points
    line_values = interpolator(line_points)

    return line_values, line_points, distances

