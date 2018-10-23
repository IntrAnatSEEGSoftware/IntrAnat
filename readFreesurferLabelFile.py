import numpy
from random import randint

def readFreesurferLabelFile(freesurfer_label_file, nLabels):
    
    # Initialize empty dictionnary
    dict_label = {}
    
    # Initialize palette with random colors or zeros (when a real input file exists)
    colors = [0,0,0]
    for x in range(nLabels):
        if freesurfer_label_file:
            colors.extend([0, 0, 0])
        else:
            colors.extend([randint(0,255), randint(0,255), randint(0,255)])
        
        
    # No file to read: return only the random colors
    if freesurfer_label_file:
        # Read file line by line
        with open(freesurfer_label_file,'r') as inf:
            for line in inf:
                # Line is empty or is a comment: skip
                if (not line) or (line[0] == "#"):
                    continue
                # Split line: must contain at least 5 elements (label index, label name, R, G, B)
                info_line = line.split()
                if (len(info_line) < 5):
                    continue
                # Get label
                dict_label.update({info_line[0]:info_line[1:]})
                # Get colors (R,G,B)
                iColor = 3 * int(info_line[0])
                colors[iColor] = int(info_line[2])
                colors[iColor+1] = int(info_line[3])
                colors[iColor+2] = int(info_line[4])
            
    return (dict_label, colors)
