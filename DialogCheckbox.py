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
        self.layoutList = QVBoxLayout();
        layout.addLayout(self.layoutList, 1, 1)
        self.checkbox = [None] * len(listOpt);
        # Create all the the checkboxes
        for i in range(0, len(listOpt)):
            self.checkbox[i] = QCheckBox()
            self.checkbox[i].setText(listOpt[i])
            self.layoutList.addWidget(self.checkbox[i])
            if (defaultSel is not None) and (i <= len(defaultSel)) and (defaultSel[i]):
                self.checkbox[i].setCheckState(Qt.Checked)

        self.setWindowTitle(titleWin)
        self.setText(titleList)

        self.setStandardButtons(QMessageBox.Cancel |QMessageBox.Ok)
        self.setDefaultButton(QMessageBox.Cancel)
        # self.setIcon(QMessageBox.Warning)


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

    dialog = DialogCheckbox(["Test 1", "Test 2", "Test 3", "Test 4"], "Export", "Select options to run:", [True, True, False, False])
    answer = dialog.exec_()

    # Example how to handle the response
    print answer

    app.quit() 
    sys.exit(app.exec_())
