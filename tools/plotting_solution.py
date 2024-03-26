import os, sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from optparse import OptionParser


parser = OptionParser()
parser.add_option("-d", "--dir",    action="store", type="string", dest="d", default="tsunami")
options, args = parser.parse_args()
folder = options.d

aspect_ratio="quad"

if aspect_ratio=="quad":
    aspect_ratio_dimensions=[4, 4, 4]
elif aspect_ratio=="rect_X":
    aspect_ratio_dimensions=[14, 4, 4]
else: 
    print("Aspcet ratio undefined")
    quit()

Gmm=1.4

print(folder)
if not os.path.isdir(folder):
    raise ValueError("Folder not valid. Run as 'python plotting.py --dir folder_name'")

i_plot = 0
def plotting(folder_name, what_plot='all'):
    global i_plot
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

    fig = plt.figure(1, figsize=(20,8))

    if what_plot=='density':
        ax_ro = fig.add_subplot(111, projection='3d')
    elif what_plot=='primitive':
        ax_ro = fig.add_subplot(141, projection='3d')
        ax_u  = fig.add_subplot(142, projection='3d')
        ax_v  = fig.add_subplot(143, projection='3d')
        ax_p  = fig.add_subplot(144, projection='3d')
    elif what_plot=='conservative':
        ax_ro = fig.add_subplot(141, projection='3d')
        ax_qx = fig.add_subplot(142, projection='3d')
        ax_qy = fig.add_subplot(143, projection='3d')
        ax_E  = fig.add_subplot(144, projection='3d')
    elif what_plot=='all':
        ax_ro = fig.add_subplot(241, projection='3d')
        ax_u  = fig.add_subplot(242, projection='3d')
        ax_v  = fig.add_subplot(243, projection='3d')
        ax_p  = fig.add_subplot(244, projection='3d')
        ax_qx = fig.add_subplot(245, projection='3d')
        ax_qy = fig.add_subplot(246, projection='3d')
        ax_E  = fig.add_subplot(247, projection='3d')
    else:
        print("Plotting option not available")
        print("It was")
        print(what_plot)
        quit()


    def plot_one_step(i):
        filename = import_solution[i]
        fig.suptitle(filename + " right:next, left:prev, up:-10, down:+10")
        mydata_solution = np.loadtxt(folder_name + filename, skiprows=headerlines_in)
        ro   = mydata_solution[:, 0]
        u    = mydata_solution[:, 1]
        v    = mydata_solution[:, 2]
        p    = mydata_solution[:, 3]
        phi  = mydata_solution[:, 4]

        RO  = np.reshape(ro , (X.shape[0], Y.shape[1]))
        U   = np.reshape(u  , (X.shape[0], Y.shape[1]))
        V   = np.reshape(v  , (X.shape[0], Y.shape[1]))
        P   = np.reshape(p  , (X.shape[0], Y.shape[1]))
        PHI = np.reshape(phi, (X.shape[0], Y.shape[1]))



        # Surface plots
        ax_ro.clear()
        ax_ro.plot_surface(X, Y, RO, facecolor=plt.cm.get_cmap('Blues')(0.7), alpha=0.7)
        ax_ro.set_title(r"$\rho$")
        ax_ro.set_xlabel('X')
        ax_ro.set_ylabel('Y')
        # ax_eta.set_zlabel(r'$\rho$')
        ax_ro.set_box_aspect(aspect_ratio_dimensions)  # Set aspect ratio
        # ax_eta.view_init(elev=30, azim=175) 

        if what_plot=='primitive' or what_plot=='all':
            ax_u.clear()
            ax_u.plot_surface(X, Y, U, cmap='viridis')
            ax_u.set_title('u')
            ax_u.set_xlabel('X')
            ax_u.set_ylabel('Y')
            ax_u.set_box_aspect(aspect_ratio_dimensions)  # Set aspect ratio
            # ax_u.set_zlabel('u')

            ax_v.clear()
            ax_v.plot_surface(X, Y, V, cmap='viridis')
            ax_v.set_title('v')
            ax_v.set_xlabel('X')
            ax_v.set_ylabel('Y')
            ax_v.set_box_aspect(aspect_ratio_dimensions)  # Set aspect ratio
            # ax_v.set_zlabel('v')

            ax_p.clear()
            ax_p.plot_surface(X, Y, P, cmap='viridis')
            ax_p.set_title('p')
            ax_p.set_xlabel('X')
            ax_p.set_ylabel('Y')
            ax_p.set_box_aspect(aspect_ratio_dimensions)  # Set aspect ratio
            # ax_p.set_zlabel('p')

        if what_plot=='conservative' or what_plot=='all':
            ax_qx.clear()
            ax_qx.plot_surface(X, Y, RO * U, cmap='viridis')
            ax_qx.set_title(r'$q_x$')
            ax_qx.set_xlabel('X')
            ax_qx.set_ylabel('Y')
            ax_qx.set_box_aspect(aspect_ratio_dimensions)  # Set aspect ratio
            # ax_qx.set_zlabel(r'$q_x$')

            ax_qy.clear()
            ax_qy.plot_surface(X, Y, RO * V, cmap='viridis')
            ax_qy.set_title(r'$q_y$')
            ax_qy.set_xlabel('X')
            ax_qy.set_ylabel('Y')
            # ax_qy.set_zlabel(r'$q_y$')

            ax_E.clear()
            ax_E.plot_surface(X, Y, P/(Gmm-1.)+0.5*RO*(U**2+V**2), cmap='viridis')
            ax_E.set_title(r'$E$')
            ax_E.set_xlabel('X')
            ax_E.set_ylabel('Y')
            # ax_E.set_zlabel(r'$E$')

    
    plot_one_step(i_plot)


    def change_file(event):
        global i_plot
        sys.stdout.flush()
        if event.key=='right':
            i_plot+=1
        elif event.key=='left':
            if i_plot>0:
                i_plot-=1
            else:
                return
        elif event.key == 'down':
            i_plot=min(i_plot+10, nfiles_solution)
        elif event.key == 'up':
            i_plot=max(i_plot-10,0)

        plot_one_step(i_plot)
        fig.canvas.draw_idle()

    fig.canvas.mpl_connect('key_press_event', change_file)
    plt.show()

    # for i in range(nfiles_solution):
    #     plot_one_step(i)


    #     plt.pause(0.1)

    # plt.show()





if __name__=='__main__':
    
    # plotting(folder,'density')
    plotting(folder,'primitive')
    # plotting(folder,'conservative')
    # plotting(folder,'all')