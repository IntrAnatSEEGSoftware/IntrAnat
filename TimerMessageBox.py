#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from PyQt4 import QtCore
from PyQt4 import QtGui


class TimerMessageBox(QtGui.QMessageBox):
    def __init__(self, timeout=3, parent=None):
        super(TimerMessageBox, self).__init__(parent)
        self.setWindowTitle("NOT FOR MEDICAL USAGE")
        self.time_to_wait = timeout
        self.changeContent()
        #self.setText("wait (closing automatically in {0} secondes.)".format(timeout))
        self.setStandardButtons(QtGui.QMessageBox.NoButton)
        self.timer = QtCore.QTimer(self)
        self.timer.setInterval(1000)
        self.connect(self.timer, QtCore.SIGNAL("timeout()"), self.changeContent)
        self.timer.start()

    def changeContent(self):
        self.setText("NOT FOR MEDICAL USAGE (closing automatically in {0} secondes.)".format(self.time_to_wait))

        if self.time_to_wait <= 0:
            self.close()
        
        self.time_to_wait -= 1

    def closeEvent(self, event):
        self.timer.stop()
        event.accept()


class Example(QtGui.QWidget):
    def __init__(self):
        super(Example, self).__init__()
        btn = QtGui.QPushButton('Button', self)
        btn.resize(btn.sizeHint())
        btn.move(50, 50)
        self.setWindowTitle('Example')
        self.connect(btn, QtCore.SIGNAL("clicked()"), self.warning)

    def warning(self):
        messagebox = TimerMessageBox(5, self)
        messagebox.exec_()


def main():
    app = QtGui.QApplication(sys.argv)
    ex = Example()
    ex.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
