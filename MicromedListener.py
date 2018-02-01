#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Reveiving data from a micromed SEEG system
#
# (c) Inserm U1216 2012-2016 - Manik Bhattacharjee
#
# License GNU GPL v3
#
# Some notes about this :
# - data arrives with a header
#typedef struct
  #{
  #char fixCode[4];
  #short int infoType;
  #int lenData;
  #} hdr_struct_t;
# followed by lenData bytes.
# We should check the header was received and is valid (see server_micromed.c)
#

from PyQt4.QtNetwork import QTcpServer, QTcpSocket
import struct

class MicromedListener(QTcpServer):
    def __init__(self, address = None, port = 5000):
        QTcpServer.__init__(self)
        self.tcpSocket = None
        self.trcHeader = None
        self.expectedBytes = 0
        self.expectedDatatype = None
        self.newConnection.connect(self.gotConnection)
        self.listen(port=port)

    def gotConnection(self):
        if self.tcpSocket is not None:
            self.tcpSocket.abort() #disconnectFromHost() woulk keep current data for processing
        self.tcpSocket = self.nextPendingConnection () #   setReadBufferSize
        print ("Receiving from : "+str(self.tcpSocket.peerAddress().toString()))
        self.tcpSocket.readyRead.connect(self.gotData)
        self.tcpSocket.disconnected.connect(self.closeSocket)
        self.tcpSocket.error.connect(lambda x: self.closeSocket())

    def gotData(self):
        bytesAvailable = self.tcpSocket.bytesAvailable()
        if self.expectedBytes == 0:
            if bytesAvailable < 10: # We don't even have the full mini header, wait for more
                return
            else:
                # Read mini header
                miniHeader = self.tcpSocket.read(10)
                if not miniHeader[:4] == 'MICM':
                    self.closeSocket("Incorrect mini-header code")
                    return
                else:
                    self.expectedDatatype = struct.unpack('h', miniHeader[4:6])
                    self.expectedBytes = struct.unpack('i', miniHeader[6:])

        elif bytesAvailable >= self.expectedBytes:#   -> is it what we expect
            newData = self.tcpSocket.read(self.expectedBytes) # QByteArray
            self.expectedBytes = 0
            if self.expectedDatatype == 0: # trcHeader
                self.decodeTrcHeader(newData)
            elif self.expectedDatatype == 1:
                self.analyzeRawData(newData)
            else:
                print('Received unknow datatype' + repr(self.expectedDatatype))
        else:
            print("Waiting for more data : expecting %s, got %s"%(str(self.expectedBytes), str(bytesAvailable)))

    def decodeTrcHeader(self, data):
      self.trcHeader = {}
      print('Got Trc Header !')

    def analyzeRawData(self, data):
      #Separate data from each electrode ?
      # Signal filtering
      # We have displayable data !
      # Should we send a signal to refresh the display (max 30fps) ?
      print('Got raw data')

    def closeSocket(self, message="no reason given"):
      print ('Closing socket : '+message)
      self.tcpSocket.abort()
      self.expectedBytes = 0
      self.expectedDatatype = -1
