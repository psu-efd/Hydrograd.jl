"""
Test the SRH_2D_Model class

Run SRH-2D simulation with pyHMT2D
"""
import json
import pyHMT2D

def run_SRH_2D():
    """ Run the SRH-2D case

    Returns
    -------

    """

    #the follow should be modified based on your installation of SRH-2D
    version = "3.6.5"
    srh_pre_path = r"C:\Program Files\SMS 13.3 64-bit\python\Lib\site-packages\srh2d_exe\SRH_Pre_Console.exe"
    srh_path = r"C:\Program Files\SMS 13.3 64-bit\python\Lib\site-packages\srh2d_exe\SRH-2D_Console_v365.exe"
    extra_dll_path = r"C:\Program Files\SMS 13.3 64-bit\python\Lib\site-packages\srh2d_exe"

    #create a SRH-2D model instance
    my_srh_2d_model = pyHMT2D.SRH_2D.SRH_2D_Model(version, srh_pre_path,
                       srh_path, extra_dll_path, faceless=False)

    #initialize the SRH-2D model
    my_srh_2d_model.init_model()

    print("Hydraulic model name: ", my_srh_2d_model.getName())
    print("Hydraulic model version: ", my_srh_2d_model.getVersion())

    #open a SRH-2D project
    my_srh_2d_model.open_project("Muncie.srhhydro")

    #run SRH-2D Pre to preprocess the case
    my_srh_2d_model.run_pre_model()

    #run the SRH-2D model's current project
    my_srh_2d_model.run_model()

    #close the SRH-2D project
    my_srh_2d_model.close_project()

    #quit SRH-2D
    my_srh_2d_model.exit_model()


def convert_SRH_2D_to_VTK():
    """ Convert SRH-2D results to VTK

    Returns
    -------

    """

    my_srh_2d_data = pyHMT2D.SRH_2D.SRH_2D_Data("Muncie.srhhydro")

    #read SRH-2D result in XMDF format (*.h5)
    #wether the XMDF result is nodal or cell center
    bNodal = False

    my_srh_2d_data.readSRHXMDFFile("Muncie_XMDFC.h5", bNodal)

    print(type(my_srh_2d_data.xmdfAllData_Cell))
    print(my_srh_2d_data.xmdfAllData_Cell.keys())

    h = my_srh_2d_data.xmdfAllData_Cell['Water_Depth_ft']
    u = my_srh_2d_data.xmdfAllData_Cell['Vel_X_ft_p_s']
    v = my_srh_2d_data.xmdfAllData_Cell['Vel_Y_ft_p_s']

    hu = h * u
    hv = h * v

    print(type(h))
    print(type(hu))
    print(type(hv))

    # Organize arrays into a dictionary
    data = {
        "h": h.tolist(),
        "hu": hu.tolist(),
        "hv": hv.tolist()
    }

    # Write to a JSON file with pretty-print
    with open("forward_simulation_initial_condition.json", "w") as json_file:
        json.dump(data, json_file, indent=4)


    #export to VTK
    #my_srh_2d_data.outputXMDFDataToVTK(bNodal, lastTimeStep=False, dir='')

if __name__ == "__main__":

    #run_SRH_2D()

    convert_SRH_2D_to_VTK()

    print("All done!")
