def readSulcusLabelTranslationFile(sulcus_label_file):
    sulc_labels = []
    with open(sulcus_label_file,'r') as inf:
        for line in inf:
            sulc_labels.append(line.split())
    sulc_labels_dict = dict((int(value), key) for (value, key) in sulc_labels)
    return sulc_labels_dict
