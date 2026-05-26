import shutil
import os
import pydicom

def structure_frames(
        src : str, 
        dst : str, 
        Nb_of_frames : int): 
    """
    Organizes frames from source directory into the destination directory based on specified parameters.

    Parameters:
    ----------
        src (str): Path to the source directory containing PET frames (dicom).
        dst (str): Path to the destination directory where frames will be organized.
        Nb_of_frames (int): The number of frames to be organized (from end frame)

    Returns:
        None
    """

    if os.path.isdir(src):

        dcm_all = os.listdir(src) 
        
        dcm = os.path.join(src, dcm_all[0])
        dcm_hdr = pydicom.dcmread(dcm, force=True, stop_before_pixels=True) 
        
        # Define nb of frames and create list of wanted dataframes
        NbFrames = dcm_hdr.NumberOfTimeSlices                           
        dfNb = [i for i in range((NbFrames+1)-Nb_of_frames, NbFrames+1)]

        # Copy .dcm files to designated frame directories
        for dcmFile in dcm_all:
            dcm = os.path.join(src,dcmFile)           
            header = pydicom.dcmread(dcm, force=True, stop_before_pixels=True)
            
            num_slices = header.get((0x0054, 0x0081), None)
            image_index = header.get((0x0054, 0x1330), None)
            FNb = (image_index.value - 1) // num_slices.value + 1

            if FNb in dfNb:
                FR_path = os.path.join(dst,"FR"+str(FNb))
                if not os.path.isdir(FR_path):
                    os.makedirs(FR_path)
                shutil.copy2(dcm, FR_path)
