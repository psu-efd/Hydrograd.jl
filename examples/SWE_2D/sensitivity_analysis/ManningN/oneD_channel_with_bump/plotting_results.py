import numpy as np
import json

import matplotlib.pyplot as plt
import cv2

#import pyHMT2D
#from pyHMT2D.Misc.tools import setNumpyArrayValueToNaN

plt.rc('text', usetex=True)  #allow the use of Latex for math expressions and equations
plt.rc('font', family='serif') #specify the default font family to be "serif"

def plot_one_step_forward_result(iStep, x_Coordinates, zb_cells, WSE_step):
    """

    Parameters
    ----------

    Returns
    -------

    """

    # create a figure and axis  
    fig, ax = plt.subplots(figsize=(6, 4))
    
    # plot bed
    plt.plot(x_Coordinates, zb_cells, color='black', label=r'Bed')
    
    # plot the wse
    plt.plot(x_Coordinates, WSE_step, color='aqua', label='WSE inversion')    
    
    # annotate the loss and iteration number
    plt.text(0.5, 0.55, 'Iteration: {}'.format(iStep), fontsize=14, color='black')
       

    # set the limit for the x and y axes
    plt.xlim([min(x_Coordinates), max(x_Coordinates)])
    plt.ylim([0.0, 0.6])

    # set x and y axes label and font size
    plt.xlabel(r'$x$ (m)', fontsize=16)
    plt.ylabel(r'$z$ (m)', fontsize=16)

    # show the ticks on both axes and set the font size
    plt.tick_params(axis='both', which='major', labelsize=14)

    # show legend, set its location, font size, and turn off the frame
    plt.legend(loc='upper right', fontsize=14, frameon=False)

    #save to file
    plt.savefig('forward_simulation_result_iter_{}.png'.format(iStep), dpi=300, bbox_inches='tight', pad_inches=0)

    #plt.show()
    plt.close()


#plot all steps of the forward simulation
def plot_all_forward_simulation_steps():
    #Load the solution in JSON format from the forward simulation
    with open('forward_simulation_results.json', 'r') as file:

        data = json.load(file)

        forward_simulation_results = data['forward_simulation_results']
        zb_cells = data['zb_cells']

    print("len(forward_simulation_results): {}".format(len(forward_simulation_results)))
    #exit()

    nCells = len(zb_cells)
    print("Number of cells: {}".format(nCells))

    #get the number of time steps
    nTimeSteps = int(len(forward_simulation_results)/3/nCells)
    print("Number of time steps: {}".format(nTimeSteps))

    #get the x coordinates 
    x_Coordinates = np.linspace(0, 25, nCells)

    #loop through all steps of the forward simulation
    for iStep in range(0, nTimeSteps):
        print("Processing step {}".format(iStep))

        h_step = forward_simulation_results[3*nCells*iStep:3*nCells*iStep+nCells]
        hu_step = forward_simulation_results[3*nCells*iStep+nCells:3*nCells*iStep+2*nCells]
        hv_step = forward_simulation_results[3*nCells*iStep+2*nCells:3*nCells*iStep+3*nCells]

        #convert to numpy array
        h_step = np.array(h_step)
        hu_step = np.array(hu_step)
        hv_step = np.array(hv_step)
        
        WSE_step = h_step + zb_cells

        print("size of h_step: {}".format(len(h_step)))       
        print("size of zb_cells: {}".format(len(zb_cells))) 
        print("size of WSE_step: {}".format(len(WSE_step)))

        #plot the inversion result against the truth
        plot_one_step_forward_result(iStep, x_Coordinates, zb_cells, WSE_step)


#make a video from the image sequences
def make_video_from_image_sequences(iStart, iEnd):

    # video file name
    output_video = 'oneD_channel_with_bump_sensitivity_forward.mp4'

    # Read the first image to determine the frame size
    first_image = cv2.imread('forward_simulation_result_iter_{}.png'.format(iStart))
    height, width, layers = first_image.shape

    # Define the video codec and create VideoWriter object
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')  # Codec for MP4
    fps = 5  # Frames per second
    video = cv2.VideoWriter(output_video, fourcc, fps, (width, height))

    # Add each image to the video
    for i in range(iStart, iEnd+1):
        #image file name
        image_path = 'forward_simulation_result_iter_{}.png'.format(i)
        
        frame = cv2.imread(image_path)
        video.write(frame)

    # Release the video writer
    video.release()

    print(f"Video saved as {output_video}")


# plot the sensitivity results
def plot_sensitivity():
   
    #Load the solution in JSON format from the forward simulation
    with open('forward_simulation_results.json', 'r') as file:

        data = json.load(file)

        forward_simulation_results = data['forward_simulation_results']
        wstill = data['wstill']
        zb_cells = data['zb_cells']

    nCells = len(zb_cells)
    print("Number of cells: {}".format(nCells))

    #get the number of time steps
    nTimeSteps = int(len(forward_simulation_results)/3/nCells)
    print("Number of time steps: {}".format(nTimeSteps))

    #only get the result from the last time step
    xi_step = forward_simulation_results[3*nCells*(nTimeSteps-1):3*nCells*(nTimeSteps-1)+nCells]
    hu_step = forward_simulation_results[3*nCells*(nTimeSteps-1)+nCells:3*nCells*(nTimeSteps-1)+2*nCells]
    hv_step = forward_simulation_results[3*nCells*(nTimeSteps-1)+2*nCells:3*nCells*(nTimeSteps-1)+3*nCells]

    #convert to numpy array
    xi_step = np.array(xi_step)
    hu_step = np.array(hu_step)
    hv_step = np.array(hv_step)

    WSE_step = xi_step + wstill

    #get the x coordinates 
    x_Coordinates = np.linspace(0, 25, nCells)        

    #number of parameters
    nParams = 3
    parameter_name = 'ManningN'

    #loop through each parameter
    for iParam in range(0, nParams):
        #load the sensitivity results from JSON file
        with open('sensitivity_results_ManningN_{}.json'.format(iParam+1), 'r') as file:
            data = json.load(file)

            dh_dparam = data['dh_dparam']
            dhu_dparam = data['dhu_dparam']
            dhv_dparam = data['dhv_dparam']

        #convert to numpy array
        dh_dparam = np.array(dh_dparam)
        dhu_dparam = np.array(dhu_dparam)
        dhv_dparam = np.array(dhv_dparam)

        #plot the sensitivity results
        fig, ax1 = plt.subplots(figsize=(6, 4))

        #plot the sensitivity on the left axis
        if iParam == 1:
            ax1.plot(x_Coordinates, dh_dparam, color='blue', linestyle='--', label=r'$\frac{\partial h}{\partial n_{1}}$')
        elif iParam == 2:
            ax1.plot(x_Coordinates, dh_dparam, color='blue', linestyle='--',
                     label=r'$\frac{\partial h}{\partial n_{2}}$')
        elif iParam == 3:
            ax1.plot(x_Coordinates, dh_dparam, color='blue', linestyle='--',
                     label=r'$\frac{\partial h}{\partial n_{3}}$')
        else:
            ax1.plot(x_Coordinates, dh_dparam, color='blue', linestyle='--',
                     label=r'$\frac{\partial h}{\partial n}$')



        #ax1.plot(x_Coordinates, dhu_dparam, color='green', label=r'$\frac{\partial hu}{\partial n}$')
        #ax1.plot(x_Coordinates, dhv_dparam, color='blue', label=r'$\frac{\partial hv}{\partial n}$')

        # Set the labels and limits for the primary axis
        ax1.set_xlabel(r'$x$ (m)', fontsize=16)
        ax1.set_ylabel('Sensitivity', fontsize=16, color='black')
        ax1.tick_params(axis='both', which='major', labelsize=14)
        ax1.set_xlim([min(x_Coordinates), max(x_Coordinates)])
        ax1.set_ylim([-1, 12])

        #plot the bed and the WSE on the right axis
        ax2 = ax1.twinx()
        ax2.plot(x_Coordinates, zb_cells, color='black', alpha=0.5, label='Bed')
        ax2.plot(x_Coordinates, WSE_step, color='aqua', label='WSE')

        ax2.set_ylim([-0.1, 0.6])

        # Set the labels for the secondary axis
        ax2.set_ylabel(r'$z$ (m)', fontsize=16, color='black')
        ax2.tick_params(axis='both', which='major', labelsize=14)

        #add a shaded region to indicate the bump
        ax2.fill_between([0, 12.5], 0.0, 0.2, color='gray', alpha=0.5)
        ax2.fill_between([12.5, 25], 0.0, 0.2, color='green', alpha=0.5)

        #add text in the shaded region (center aligned)
        ax2.text(4, 0.1, 'zone $n_1$ = 0.02', fontsize=14, color='black', ha='center')
        ax2.text(19, 0.1, 'zone $n_2$ = 0.03', fontsize=14, color='black', ha='center')

        # Add legends for both axes
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1, labels1, loc='upper left', fontsize=14, frameon=False)
        ax2.legend(lines2, labels2, loc='upper right', fontsize=14, frameon=False)

        #save to file
        plt.savefig('sensitivity_ManningN_{}.png'.format(iParam+1), dpi=300, bbox_inches='tight', pad_inches=0)

        #plt.show()
        


if __name__ == "__main__":

    #show current working directory
    import os
    print("Current working directory: {}".format(os.getcwd()))

    #change the working directory to the directory where this script is located
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    print("Changed working directory to: {}".format(os.getcwd()))

    iStart = 1
    iEnd = 100

    #plot all the inversion results
    #plot_all_forward_simulation_steps()

    #make a video
    #make_video_from_image_sequences(iStart, iEnd)

    #plot the sensitivity
    plot_sensitivity()
    

    print("Done!")