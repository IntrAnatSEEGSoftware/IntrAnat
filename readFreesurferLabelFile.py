def readFreesurferLabelFile(freesurfer_label_file):
    dict_label = {}
    with open(freesurfer_label_file,'r') as inf:
       for line in inf:
           info_line = line.split()
           dict_label.update({info_line[0]:info_line[1:]})
    return dict_label
