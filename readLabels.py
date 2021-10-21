#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# (c) 2018-2021 Inserm
#
# License GNU GPL v3

def readLabels(sulcus_label_file):
    """
        Read sulcus label files as dictionary
    """
    sulc_labels = []
    with open(sulcus_label_file,'r') as inf:
        for line in inf:
            sulc_labels.append(line.split())
    sulc_labels_dict = dict((int(value), key) for (value, key) in sulc_labels)
    return sulc_labels_dict
