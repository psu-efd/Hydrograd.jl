"""
Test the SRH_2D_Model class

Run SRH-2D simulation with pyHMT2D
"""
import json

import numpy as np
import pyHMT2D

if __name__ == "__main__":

    #number of cells
    nCells = 20

    h = np.zeros(nCells)

    for iCell in range(nCells):
        if iCell < nCells/2:
            h[iCell] = 1.0
        else:
            h[iCell] = 0.0

    q_x = h * 0.0
    q_y = h * 0.0

    print('h', h, 'q_x', q_x, 'q_y', q_y)

    # Organize arrays into a dictionary
    data = {
        "h": [h.tolist()],
        "q_x": [q_x.tolist()],
        "q_y": [q_y.tolist()]
    }

    # Write to a JSON file with pretty-print
    with open("forward_simulation_initial_condition.json", "w") as json_file:
        json.dump(data, json_file, indent=4)

    print("All done!")
