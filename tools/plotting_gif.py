import numpy as np
import matplotlib.pyplot as plt
import imageio
from PIL import Image
import os, sys
from mpl_toolkits.mplot3d import Axes3D
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-d", "--dir",    action="store", type="string", dest="d", default="tsunami")
options, args = parser.parse_args()
folder_name = options.d

print(folder_name)
if not os.path.isdir(folder_name):
    raise ValueError("Folder not valid. Run as 'python plotting.py --dir folder_name'")

i_plot = 0


folder_name = folder_name + "/"
    
delimiter_in = ' '
headerlines_in = 1

files_mesh = folder_name + 'MESH_*.dat'
files_solution = folder_name + 'SOLUTION_*.dat'

# Importing the MESH DATA
import_mesh = [file for file in os.listdir(folder_name) if file.startswith('MESH_')]
import_mesh.sort()
nfiles_mesh = len(import_mesh)
for i in range(nfiles_mesh):
    if i == 0:
        filename = import_mesh[i]
        mydata_mesh = np.loadtxt(folder_name + filename, delimiter=delimiter_in, skiprows=headerlines_in)
        x = mydata_mesh
        x_ini = min(x)
        x_end = max(x)
    if i == 1:
        filename = import_mesh[i]
        mydata_mesh = np.loadtxt(folder_name + filename, delimiter=delimiter_in, skiprows=headerlines_in)
        y = mydata_mesh
        y_ini = min(y)
        y_end = max(y)

X, Y = np.meshgrid(x, y)

# Importing the SOLUTION DATA
import_solution = [file for file in os.listdir(folder_name) if file.startswith('SOLUTION_OCT')]
import_solution.sort()
nfiles_solution = len(import_solution)


def plot_one_step(i):
    perc_simul = i/nfiles_solution
    fig = plt.figure(1, figsize=(17,12))
    ax_ro = fig.add_subplot(111, projection='3d')
    filename = import_solution[i]
    fig.suptitle(filename)
    mydata_solution = np.loadtxt(folder_name + filename, skiprows=headerlines_in)
    ro  = mydata_solution[:, 0]
    u   = mydata_solution[:, 1]
    v   = mydata_solution[:, 2]
    p   = mydata_solution[:, 3]
    phi = mydata_solution[:, 4]

    RO  = np.reshape(ro,  (X.shape[0], Y.shape[1]))
    U   = np.reshape(u,   (X.shape[0], Y.shape[1]))
    V   = np.reshape(v,   (X.shape[0], Y.shape[1]))
    P   = np.reshape(p,   (X.shape[0], Y.shape[1]))
    PHI = np.reshape(phi, (X.shape[0], Y.shape[1]))


    # Surface plots
    ax_ro.clear()
    ax_ro.plot_surface(X, Y, RO, facecolor=plt.cm.get_cmap('Blues')(0.7), alpha=0.7)
    ax_ro.set_title(r"$\rho$")
    ax_ro.set_xlabel('X')
    ax_ro.set_ylabel('Y')
    # ax_eta.set_zlabel(r'$\eta$')
    ax_ro.set_box_aspect([14, 4, 4])  # Set aspect ratio
    ax_ro.view_init(elev=25*(1-perc_simul)+5, azim=-90*perc_simul-45*(1-perc_simul)) 
    fig.tight_layout()
    return fig



def fig2img(fig):
    """Convert a Matplotlib figure to a PIL Image."""
    fig.canvas.draw()
    buf = fig.canvas.tostring_rgb()
    width, height = fig.canvas.get_width_height()
    return Image.frombytes("RGB", (width, height), buf, "raw", "RGB")



# Number of frames in the GIF
num_frames = 5

# Create a sequence of figures
images = []
for i_plot in range(nfiles_solution):
    print("Plotting solution %03d out of %d"%(i_plot, nfiles_solution),end="\r")
    fig = plot_one_step(i_plot)
    images.append(fig2img(fig))
    plt.close('all')
print("")

# Create a GIF from the figures
gif_filename = "output.gif"
imageio.mimsave(gif_filename, images, fps=num_frames)

# Display the GIF filename
print(f"GIF saved to: {gif_filename}")
