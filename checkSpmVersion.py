def checkSpmVersion(spm_path):
        import os
        op_cont = open(spm_path+os.sep+'Contents.m')
        line = op_cont.readline()
        line = op_cont.readline()
        spm_version = line.split(' ')[3]
        return spm_version
