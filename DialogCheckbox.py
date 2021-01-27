import sys
from soma.qt_gui.qt_backend import QtGui, QtCore
from soma.qt_gui.qt_backend.QtCore import *
from soma.qt_gui.qt_backend.QtGui import *


class DialogCheckbox(QMessageBox):

    def __init__(self, listOpt, titleWin="", titleList="", defaultSel=None):
        super(DialogCheckbox, self).__init__()

        # Get the Layout of the MessageBox
        layout = self.layout()
        # Create a layout to contain all the checkboxes
        layoutList = QGridLayout();
        layout.addLayout(layoutList, 1, 1)
        self.checkbox = [None] * len(listOpt);
        # Create all the the checkboxes
        Nrows = 30
        for i in range(len(listOpt)):
            self.checkbox[i] = QCheckBox()
            self.checkbox[i].setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Expanding)
            self.checkbox[i].setText(listOpt[i])
            layoutList.addWidget(self.checkbox[i], i%Nrows, i/Nrows)
            if (defaultSel is not None) and (i <= len(defaultSel)) and (defaultSel[i]):
                self.checkbox[i].setCheckState(Qt.Checked)
        # Configure window
        self.setWindowTitle(titleWin)
        self.setText(titleList)
        self.setStandardButtons(QMessageBox.Cancel | QMessageBox.Ok)
        self.setDefaultButton(QMessageBox.Cancel)
        
    def exec_(self, *args, **kwargs):
        """
        Override the exec_ method so you can return the value of the checkbox
        """
        # Wait for user input
        res = QMessageBox.exec_(self, *args, **kwargs)
        # Check selection of checkboxes
        if res == QMessageBox.Ok:
            resList = [None] * len(self.checkbox);
            for i in range(0, len(self.checkbox)):
                resList[i] = self.checkbox[i].isChecked()
        else:
            resList = None
        return resList


# TEST FUNCTIONS
if __name__ == '__main__':
    app = QApplication(sys.argv)

    dialog = DialogCheckbox(["XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX", "XXX_XXXX_XXXX"], "Export", "Select options to run:", None)
    answer = dialog.exec_()
    print answer

    app.quit() 
    sys.exit(app.exec_())
