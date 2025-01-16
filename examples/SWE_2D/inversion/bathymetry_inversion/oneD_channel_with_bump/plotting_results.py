import numpy as np
import meshio
import vtk
from vtk.util import numpy_support as VN
import h5py
import json

import matplotlib.pyplot as plt
import cv2

#import pyHMT2D
#from pyHMT2D.Misc.tools import setNumpyArrayValueToNaN

plt.rc('text', usetex=True)  #allow the use of Latex for math expressions and equations
plt.rc('font', family='serif') #specify the default font family to be "serif"

def plot_inversion_result(iter_number, loss_total, x_truth, zb_inversion, wse_inversion, zb_truth, wse_truth):
    """

    Parameters
    ----------
    zb_inversion: numpy array
        The inversion result of the bathymetry
    zb_truth: numpy array
        The truth of the bathymetry

    Returns
    -------

    """

    # create a figure and axis  
    fig, ax = plt.subplots(figsize=(6, 4))
    
    # plot bed
    plt.plot(x_truth, zb_inversion, color='black', label=r'$z_b$ inversion')
    plt.plot(x_truth, zb_truth, color='black', linestyle='--', label=r'$z_b$ truth')

    # plot the wse
    plt.plot(x_truth, wse_inversion, color='aqua', label='WSE inversion')    
    plt.plot(x_truth, wse_truth, color='aqua', linestyle='--', label='WSE truth')    

    # annotate the loss and iteration number
    plt.text(0.5, 0.55, 'Iteration: {}'.format(iter_number), fontsize=14, color='black')
    plt.text(0.5, 0.5, 'Loss: {:.4f}'.format(loss_total), fontsize=14, color='black')
    

    # set the limit for the x and y axes
    plt.xlim([min(x_truth), max(x_truth)])
    plt.ylim([-0.1, 0.6])

    # set x and y axes label and font size
    plt.xlabel(r'$x$ (m)', fontsize=16)
    plt.ylabel(r'$z$ (m)', fontsize=16)

    # show the ticks on both axes and set the font size
    plt.tick_params(axis='both', which='major', labelsize=14)

    # show legend, set its location, font size, and turn off the frame
    plt.legend(loc='upper right', fontsize=14, frameon=False)

    #save to file
    plt.savefig('oneD_channel_with_bump_inversion_result_iter_{}.png'.format(iter_number), dpi=300, bbox_inches='tight', pad_inches=0)

    #plt.show()
    plt.close()


#plot all the inversion results
def plot_all(iStart, iEnd):
    #Load the solution in JSON format from the forward simulation
    with open('forward_simulation_solution_truth.json', 'r') as file:

        data = json.load(file)

        wse_truth = data['wse_truth']
        zb_truth = data['zb_cell_truth']

    #x coordinates of the truth
    x_truth = np.linspace(0, 25, len(zb_truth))

    #open inversion_results_losses_and_ODE_solution.json file
    with open('inversion_results_losses_and_ODE_solution.json', 'r') as file:
        data = json.load(file)

        #get the loss history
        nIter = data['nIter']
        loss_totals = data['loss_totals']
        loss_preds = data['loss_preds']
        loss_pred_WSEs = data['loss_pred_WSEs']
        loss_pred_uvs = data['loss_pred_uvs']
        loss_bounds = data['loss_bounds']
        loss_slopes = data['loss_slopes']

        inverted_WSEs = data['inverted_WSEs']
        
    iterations = np.arange(1, nIter+1)

    #make sure iEnd is within iterations
    if iEnd > nIter:
        print("iEnd is greater than the number of iterations.")
        exit()

    #loop through inversion_callback_save_iter_i files to get the solution
    for i in range(iStart, iEnd+1):
        print("Processing iteration {}".format(i))
        with open('inversion_callback_save_iter_{}.json'.format(i), 'r') as file:
            data = json.load(file)

            zb_inversion = data['theta']
            iter_number = data['iter_number']
            loss_total = data['loss_total']

        WSE_inversion = inverted_WSEs[i-1]

        #plot the inversion result against the truth
        plot_inversion_result(iter_number, loss_total, x_truth, zb_inversion, WSE_inversion, zb_truth, wse_truth)


#make a video from the image sequences
def make_video_from_image_sequences(iStart, iEnd):

    # video file name
    output_video = 'oneD_channel_with_bump_inversion.mp4'

    # Read the first image to determine the frame size
    first_image = cv2.imread('oneD_channel_with_bump_inversion_result_iter_{}.png'.format(iStart))
    height, width, layers = first_image.shape

    # Define the video codec and create VideoWriter object
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')  # Codec for MP4
    fps = 5  # Frames per second
    video = cv2.VideoWriter(output_video, fourcc, fps, (width, height))

    # Add each image to the video
    for i in range(iStart, iEnd+1):
        #image file name
        image_path = 'oneD_channel_with_bump_inversion_result_iter_{}.png'.format(i)
        
        frame = cv2.imread(image_path)
        video.write(frame)

    # Release the video writer
    video.release()

    print(f"Video saved as {output_video}")


#plot the loss history
def plot_loss_history(bTotalLossOnly=False):

    #open inversion_results_losses_and_ODE_solution.json file
    with open('inversion_results_losses_and_ODE_solution.json', 'r') as file:
        data = json.load(file)

        #get the loss history
        nIter = data['nIter']
        loss_totals = data['loss_totals']
        loss_preds = data['loss_preds']
        loss_pred_WSEs = data['loss_pred_WSEs']
        loss_pred_uvs = data['loss_pred_uvs']
        loss_bounds = data['loss_bounds']
        loss_slopes = data['loss_slopes']
        
    iterations = np.arange(1, nIter+1)

    #print("loss_totals: ", loss_totals)
    #print("loss preds: ", loss_preds)
    #print("loss bounds: ", loss_bounds)
    #print("loss slopes: ", loss_slopes)

    #make the plot
    fig, ax = plt.subplots(figsize=(6, 4))

    
    if bTotalLossOnly:
        plt.plot(iterations, loss_totals, color='black', label='Total loss')    

    else:
        plt.plot(iterations, loss_totals, color='black', alpha=0.5, label='Total loss')    
        plt.plot(iterations, loss_preds, color='black', linestyle='-.', marker='o', markersize=1, label='Loss (prediction)')
        plt.plot(iterations, loss_bounds, color='red', linestyle='--', label='Loss (bounds)')
        plt.plot(iterations, loss_slopes, color='blue', linestyle='--', label='Loss (slope)')
        

    # set x and y axes label and font size
    plt.xlabel('Iteration', fontsize=16)
    plt.ylabel('Loss', fontsize=16)

     # set the limit for the x and y axes
    plt.xlim([0, len(iterations)])
    if bTotalLossOnly:
        plt.ylim([1e-3, 0.1])
    else:
        plt.ylim([1e-6, 0.1])
    

    # make y axis log scale
    plt.yscale('log')

     # show the ticks on both axes and set the font size
    plt.tick_params(axis='both', which='major', labelsize=14)

    

    #save to file
    if bTotalLossOnly:
        # show legend, set its location, font size, and turn off the frame
        plt.legend(loc='lower left', fontsize=14, frameon=True, facecolor='white')

        plt.savefig('oneD_channel_with_bump_inversion_loss_history_total_loss_only.png', dpi=300, bbox_inches='tight', pad_inches=0)
    else:
        # show legend, set its location, font size, and turn off the frame
        plt.legend(loc='lower left', fontsize=14, frameon=True, facecolor='white')

        plt.savefig('oneD_channel_with_bump_inversion_loss_history_all_losses.png', dpi=300, bbox_inches='tight', pad_inches=0)

    plt.close()


if __name__ == "__main__":

    #show current working directory
    import os
    print("Current working directory: {}".format(os.getcwd()))

    #change the working directory to the directory where this script is located
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    print("Changed working directory to: {}".format(os.getcwd()))

    iStart = 1
    iEnd = 320

    #plot all the inversion results
    #plot_all(iStart, iEnd)

    #make a video
    make_video_from_image_sequences(iStart, iEnd)

    #plot the loss history
    #plot_loss_history(bTotalLossOnly=True)
    #plot_loss_history(bTotalLossOnly=False)

    

    print("Done!")