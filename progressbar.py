# Persisent wait bar
#
# (c) Inserm U1216 2018 - Francois Tadel
# License GNU GPL v3

import sys
from soma.qt_gui.qt_backend import QtCore
from soma.qt_gui.qt_backend import QtGui


# ===== PROGRESS WINDOW =====
class ProgressDialog(QtGui.QDialog):
    
    def __init__(self, parent=None, progressText='...', progressTitle='Process', cancelEnabled=True):
        super(ProgressDialog, self).__init__(parent)
        # Dialog window
        self.setWindowTitle(progressTitle)
        self.resize(300, 100)
        # Layout: Vertical list
        vbox = QtGui.QVBoxLayout()
        # Text label
        self.labelMessage = QtGui.QLabel(progressText)
        self.labelMessage.setAlignment(QtCore.Qt.AlignCenter)
        vbox.addWidget(self.labelMessage)
        # Progress bar
        self.progressBar = QtGui.QProgressBar()
        self.progressBar.setAlignment(QtCore.Qt.AlignCenter)
        self.progressBar.setObjectName("progressBar")
        self.progressBar.setMinimum(0)
        self.progressBar.setMaximum(0)
        self.progressBar.setMaximum(0)
        vbox.addWidget(self.progressBar)
        # Cancel button
        self.cancelEnabled = cancelEnabled
        self.buttonCancel = QtGui.QPushButton("Cancel")
        self.buttonCancel.setObjectName("buttonCancel")
        self.buttonCancel.setEnabled(self.cancelEnabled)
        vbox.addWidget(self.buttonCancel)
        # Set layout
        self.setLayout(vbox)
        # Associated worker
        self.worker = None

    # Start a thread
    def start(self, func, isModal=False):
        self.worker = ProgressThread(func)
        self.connect(self.worker, QtCore.SIGNAL("PROGRESS"), self.setProgress)
        self.connect(self.worker, QtCore.SIGNAL("PROGRESS_TEXT"), self.setText)
        self.connect(self.worker, QtCore.SIGNAL("finished()"), self.terminated)
        self.worker.start()
        if self.cancelEnabled:
            self.buttonCancel.clicked.connect(self.cancel)
        if isModal:
            res = self.exec_()
            return self.worker.output()
        else:
            self.show()
            return None
        
    # Callback: Cancel button clicked
    def cancel(self):
        if self.worker is not None:
            self.worker.terminate()
        self.hide()
    
    # Event: Progress dialog closed
    def closeEvent(self, event):
        if self.cancelEnabled:
            QtGui.QFrame.closeEvent(self, event)
            self.cancel()
        
    # Event: Thread terminated
    def terminated(self):
        self.hide()
         
    # Set progress bar value
    def setProgress(self, progress):
        # Valid progress bar values
        if (progress >= 0) and (progress <= 100):
            self.progressBar.setMaximum(100)
            self.progressBar.setValue(progress)
        # Infinite progress bar
        else:
            self.progressBar.setMaximum(0)
            self.progressBar.setValue(0)           
    
    # Set progress window text
    def setText(self, text):
        self.labelMessage.setText(text)
    
    @staticmethod
    def call(func, isModal=True, parent=None, progressText='...', progressTitle='Process', cancelEnabled=True):
        progress = ProgressDialog(parent, progressText, progressTitle, cancelEnabled)
        res = progress.start(func, isModal)
        return res
        
        
# ===== PROCESSING THREAD =====
class ProgressThread(QtCore.QThread):
    """ Executes a python function in a separate thread. """
    def __init__(self, func, parent=None):
        """:param func : The function that will be run in a separate thread
           :param parent : a parent object for the QThread (default is None)
        """
        QtCore.QThread.__init__(self,parent)
        self.func = func
        self.out = None

    def output(self):
        """ Returns the output value of the function when execution is terminated"""
        return self.out

    def run(self):
        """
           Reimplementation of the run function of the QThread.
           This SHOULD NOT BE CALLED DIRECTLY as it would run in the current thread.
           Run self.start() to start the execution in a separate thread
        """
        self.out = self.func(self)


# ===== EXAMPLE =====
class Example(QtGui.QWidget):
    def __init__(self):
        super(Example, self).__init__()
        btn = QtGui.QPushButton('Start', self)
        btn.resize(btn.sizeHint())
        btn.move(50, 50)
        self.setWindowTitle('Example')
        self.connect(btn, QtCore.SIGNAL("clicked()"), self.test)

    def test(self):
        res = ProgressDialog.call(self.activeWait, False, self, "Processing...", "Example process")
        print res
        
    def activeWait(self, thread):
        import time
        for n in range(0,200):
            print n
            if (n > 20):
                thread.emit(QtCore.SIGNAL("PROGRESS"), n)
            thread.emit(QtCore.SIGNAL("PROGRESS_TEXT"), "Processing event #{}".format(n))
            n += 1
            time.sleep(0.1)


# ===== MAIN =====
def main():
    app = QtGui.QApplication(sys.argv)
    ex = Example()
    ex.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()


