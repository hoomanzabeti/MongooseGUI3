# GUI SUCCESSFULLY DISPLAYS IN THIS VERSION
# This application requires Qt environment: http://pyqt.sourceforge.net/Docs/PyQt4/installation.html
# NEED TO DO:
# -DOWNLOAD QSOPTEX

import os
import sys
import logging
import shelve
from MetaMerge import *
from PyQt4 import QtCore, QtGui, uic

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class QtHandler(logging.Handler):
    def __init__(self):
        logging.Handler.__init__(self)
    def emit(self, record):
        record = self.format(record)
        if record: XStream.stdout().write('%s\n'%record)
        # originally: XStream.stdout().write("{}\n".format(record))

logger = logging.getLogger(__name__)
handler = QtHandler()
handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)

class XStream(QtCore.QObject):
    _stdout = None
    _stderr = None
    messageWritten = QtCore.pyqtSignal(str)
    def flush( self ):
        pass
    def fileno( self ):
        return -1
    def write( self, msg ):
        if ( not self.signalsBlocked() ):
            self.messageWritten.emit(str(msg))
    @staticmethod
    def stdout():
        if ( not XStream._stdout ):
            XStream._stdout = XStream()
            sys.stdout = XStream._stdout
        return XStream._stdout
    @staticmethod
    def stderr():
        if ( not XStream._stderr ):
            XStream._stderr = XStream()
            sys.stderr = XStream._stderr
        return XStream._stderr



class Ui_MainWindow(object):
    def setupUi(self, MainWindow):

        # size constraints
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(600, 470)

        # central widget ( main window for buttons )
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))

        # output
        self.resultOutput = QtGui.QTextBrowser(self.centralwidget)
        self.resultOutput.setGeometry(QtCore.QRect(230, 50, 361, 370))
        self.resultOutput.setObjectName(_fromUtf8("resultOutput"))
        XStream.stdout().messageWritten.connect( self.resultOutput.insertPlainText )
        XStream.stderr().messageWritten.connect( self.resultOutput.insertPlainText )

        # intializes executeAction button
        self.executeAction = QtGui.QPushButton(self.centralwidget)
        self.executeAction.setGeometry(QtCore.QRect(20, 330, 200, 32))
        self.executeAction.setObjectName(_fromUtf8("executeAction"))
        self.executeAction.setEnabled(False)

        # intializes saveContent button
        self.saveContent = QtGui.QPushButton(self.centralwidget)
        self.saveContent.setGeometry(QtCore.QRect(20, 360, 200, 32))
        self.saveContent.setObjectName(_fromUtf8("saveContent"))
        self.saveContent.clicked.connect(self.save_content)
        self.saveContent.setEnabled(False)

        # intializes saveModel button
        self.saveModel = QtGui.QPushButton(self.centralwidget)
        self.saveModel.setGeometry(QtCore.QRect(20, 390, 200, 32))
        self.saveModel.setObjectName(_fromUtf8("saveModel"))
        self.saveModel.clicked.connect(self.save_model)
        self.saveModel.setEnabled(False)


        # initializes file select button
        self.selectFile = QtGui.QPushButton(self.centralwidget)
        self.selectFile.setGeometry(QtCore.QRect(20, 50, 200, 32))
        self.selectFile.setObjectName(_fromUtf8("selectFile"))
        self.selectFile.clicked.connect(self.open)


        # intializes first dropdown menu, which chooses the model
        #self.chooseModel = QtGui.QTextEdit(self.centralwidget)
        self.chooseModel = QtGui.QLineEdit(self.centralwidget)
        self.chooseModel.setGeometry(QtCore.QRect(22, 95, 145, 33))
        self.chooseModel.setObjectName(_fromUtf8("chooseModel"))
        self.chooseModel.setPlaceholderText('Enter Model Name')

        # initializes select model button
        self.selectModel = QtGui.QPushButton(self.centralwidget)
        self.selectModel.setGeometry(QtCore.QRect(162, 91, 59, 45))
        self.selectModel.setObjectName(_fromUtf8("selectModel"))
        self.selectModel.clicked.connect(self.select_model)
        self._running = False

        # intializes title of app
        self.title = QtGui.QLabel(self.centralwidget)
        self.title.setGeometry(QtCore.QRect(120, 5, 420, 31))
        font = QtGui.QFont()
        font.setPointSize(22)
        self.title.setFont(font)
        self.title.setAlignment(QtCore.Qt.AlignCenter)
        self.title.setObjectName(_fromUtf8("title"))

        # intializes second dropdown menu, chooses a function
        self.chooseFunction2 = QtGui.QComboBox(self.centralwidget)
        self.chooseFunction2.setGeometry(QtCore.QRect(20, 174, 200, 50))
        self.chooseFunction2.setObjectName(_fromUtf8("chooseFunction2"))
        self.chooseFunction2.addItem(_fromUtf8(""))
        self.chooseFunction2.addItem(_fromUtf8(""))
        self.chooseFunction2.addItem(_fromUtf8(""))
        self.chooseFunction2.addItem(_fromUtf8(""))
        self.chooseFunction2.addItem(_fromUtf8(""))
        self.chooseFunction2.addItem(_fromUtf8(""))
        self.chooseFunction2.addItem(_fromUtf8(""))
        self.chooseFunction2.addItem(_fromUtf8(""))
        self.chooseFunction2.setVisible(False)
        self.chooseFunction2.currentIndexChanged.connect(self.hide2)

        # intializes first dropdown menu, chooses a function
        self.chooseFunction1 = QtGui.QComboBox(self.centralwidget)
        self.chooseFunction1.setGeometry(QtCore.QRect(20, 137, 200, 50))
        self.chooseFunction1.setEditable(True)
        self.chooseFunction1.setObjectName(_fromUtf8("chooseFunction1"))
        self.chooseFunction1.addItem(_fromUtf8(""))
        self.chooseFunction1.addItem(_fromUtf8(""))
        self.chooseFunction1.addItem(_fromUtf8(""))
        self.chooseFunction1.addItem(_fromUtf8(""))
        self.chooseFunction1.addItem(_fromUtf8(""))
        self.chooseFunction1.addItem(_fromUtf8(""))
        self.chooseFunction1.addItem(_fromUtf8(""))
        self.chooseFunction1.addItem(_fromUtf8(""))
        self.chooseFunction1.addItem(_fromUtf8(""))
        self.chooseFunction1.currentIndexChanged.connect(self.hide1)

        # intializes third dropdown menu, chooses a function
        self.chooseFunction3 = QtGui.QComboBox(self.centralwidget)
        self.chooseFunction3.setGeometry(QtCore.QRect(20, 198, 200, 50))
        self.chooseFunction3.setObjectName(_fromUtf8("chooseFunction3"))
        self.chooseFunction3.addItem(_fromUtf8(""))
        self.chooseFunction3.addItem(_fromUtf8(""))
        self.chooseFunction3.setVisible(False)

        # calls MONGOOSE functions
        self.executeAction.clicked.connect(self.chooseFunction)
        '''
        # intializes third dropdown menu, chooses a function
        self.varInputA = QtGui.QTextEdit(self.centralwidget)
        self.varInputA.setGeometry(QtCore.QRect(20, 220, 200, 50))
        self.varInputA.setObjectName(_fromUtf8("chooseFunction3"))
        '''
        '''
        global dropdown1_open
        dropdown1_open = False
        global dropdown2_open
        dropdown2_open = False
        global dropdown3_open
        dropdown3_open = False
        '''

        self.textEdit1 = QtGui.QTextEdit(self.centralwidget)
        self.textEdit1.setGeometry(QtCore.QRect(100, 240, 115, 25))
        self.textEdit1.setObjectName(_fromUtf8("textEdit"))
        self.label1 = QtGui.QLabel(self.centralwidget)
        self.label1.setGeometry(QtCore.QRect(30, 245, 56, 13))
        self.label1.setObjectName(_fromUtf8("label"))
        self.textEdit1.setVisible(False)
        self.label1.setVisible(False)

        '''
        self.textEdit2 = QtGui.QTextEdit(self.centralwidget)
        self.textEdit2.setGeometry(QtCore.QRect(100, 270, 115, 25))
        self.textEdit2.setObjectName(_fromUtf8("textEdit"))
        self.label2 = QtGui.QLabel(self.centralwidget)
        self.label2.setGeometry(QtCore.QRect(30, 275, 56, 13))
        self.label2.setObjectName(_fromUtf8("label"))
        self.textEdit2.setVisible(False)
        self.label2.setVisible(False)

        self.textEdit3 = QtGui.QTextEdit(self.centralwidget)
        self.textEdit3.setGeometry(QtCore.QRect(100, 300, 115, 25))
        self.textEdit3.setObjectName(_fromUtf8("textEdit"))
        self.label3 = QtGui.QLabel(self.centralwidget)
        self.label3.setGeometry(QtCore.QRect(30, 305, 56, 13))
        self.label3.setObjectName(_fromUtf8("label"))
        self.textEdit3.setVisible(False)
        self.label3.setVisible(False)
        '''

        # intializes central widget
        MainWindow.setCentralWidget(self.centralwidget)
        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)


    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow", None))
        self.executeAction.setText(_translate("MainWindow", "Execute Action", None))
        self.saveContent.setText(_translate("MainWindow", "Save Contents", None))
        self.selectFile.setText(_translate("MainWindow", "Upload File", None))
        self.saveModel.setText(_translate("MainWindow", "Save Model", None))
        self.selectModel.setText(_translate("MainWindow", "Enter", None))
        #self.chooseModel.setItemText(0, _translate("MainWindow", "<Upload Model>", None))
        #self.chooseModel.setItemText(1, _translate("MainWindow", "model", None))
        self.title.setText(_translate("MainWindow", "MONGOOSE Exact Arithmetic Toolbox", None))
        self.title.setToolTip(_translate("MainWindow", "<html><head/><body><p> <b>MONGOOSE Author:</b><br> Leonid Chindelevitch <br><br> <b>GUI Author:</b  > <br> Christopher Le </p></body></html>", None))
        self.chooseFunction2.setItemText(0, _translate("MainWindow", "<Choose>", None))
        self.chooseFunction2.setItemText(1, _translate("MainWindow", "name", None))
        self.chooseFunction2.setItemText(2, _translate("MainWindow", "pairs", None))
        self.chooseFunction2.setItemText(3, _translate("MainWindow", "reductionStatus", None))
        self.chooseFunction2.setItemText(4, _translate("MainWindow", "reversible", None))
        self.chooseFunction2.setItemText(5, _translate("MainWindow", "length", None))
        self.chooseFunction2.setItemText(6, _translate("MainWindow", "species", None))
        self.chooseFunction2.setItemText(7, _translate("MainWindow", "external", None))
        self.chooseFunction1.setToolTip(_translate("MainWindow", "<html><head/><body><p>If the desired function requires an [index] parameter, please specify what index you would like to analyze.</p></body></html>", None))
        self.chooseFunction1.setItemText(0, _translate("MainWindow", "<Choose>", None))
        self.chooseFunction1.setItemText(1, _translate("MainWindow", "reduceNetwork", None))
        self.chooseFunction1.setItemText(2, _translate("MainWindow", "addReaction", None))
        self.chooseFunction1.setItemText(3, _translate("MainWindow", "findEssentialReaction", None))
        self.chooseFunction1.setItemText(4, _translate("MainWindow", "findBiomassReaction", None))
        self.chooseFunction1.setItemText(5, _translate("MainWindow", "reactions[index]", None))
        self.chooseFunction1.setItemText(6, _translate("MainWindow", "reactionSubsets[index]", None))
        self.chooseFunction1.setItemText(7, _translate("MainWindow", "metabolites[index]", None))
        self.chooseFunction1.setItemText(8, _translate("MainWindow", "printReactionFormula[index]", None))
        self.chooseFunction3.setItemText(0, _translate("MainWindow", "<Choose>", None))
        self.chooseFunction3.setItemText(1, _translate("MainWindow", "name", None))
        self.label1.setText(_translate("MainWindow", "Var Name", None))
        #self.label2.setText(_translate("MainWindow", "Var Name", None))
        #self.label3.setText(_translate("MainWindow", "Var Name", None))

    def select_model(self):
        global modelInput
        # checks valid input
        keyList = list(s.keys())

        global modelInput
        modelInput = self.chooseModel.text()

        if modelInput in keyList:
            print('Valid option, please select functions from the dropdown menu(s)')
            self.saveModel.setEnabled(True)
            self.executeAction.setEnabled(True)
            self.saveContent.setEnabled(True)
            self.selectModel.setEnabled(False)
            self.selectFile.setEnabled(False)
            self._running = False
        else:
            print('Not a valid option, please choose again')
            self._running = True

    def open (self):
        #filename = QtGui.QFileDialog().setFilter("Text files (*.txt)").getOpenFileName(self.centralwidget, 'Open File', '.')
        filename = QtGui.QFileDialog.getOpenFileName(self.centralwidget, 'Open File', '.')
        global s
        s = shelve.open(str(filename))
        print('Please type which of the following keys you would like to analyze: ')
        #print 'PRINTING QLINEEDIT CURRENT TEXT'
        #print str(self.chooseModel.text())
        print((list(s.keys())))
        #print 'Path file :', filename
        #print 'Please select a following model from the first dropdown menu'
        #raw_input("Press Enter to continue...")
        self._running = True
        #loops until select_model button has been clicked
        while self._running:
            QtGui.qApp.processEvents()

        #apostrophe = "''"
        #desiredModel = str(apostrophe + self.chooseModel.text() + apostrophe)
        desiredModel = str(modelInput)
        global model
        model = s[desiredModel]

    # saves contents of console to desired text file
    def save_content (self):
        choice = QtGui.QMessageBox.question(QtGui.QMainWindow(), 'Save Contents', 'Are you sure you want to save the contents?', QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        if choice == QtGui.QMessageBox.Yes:
            name = QtGui.QFileDialog.getSaveFileName(QtGui.QMainWindow(), 'Save File')
            file = open(name,'w')
            text = self.resultOutput.toPlainText()
            file.write(str(self.resultOutput.toPlainText()))
            file.close()
        else:
            pass


    # saves model
    def save_model (self):
        #!!!! prompt for key then assign
        choice = QtGui.QMessageBox.question(QtGui.QMainWindow(), 'Save Model', 'Are you sure you want to save the model?', QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        if choice == QtGui.QMessageBox.Yes:
            savedModel = "'" + modelInput + "'"
            s[modelInput] = model
            print(("Saving model: " + modelInput))
        else:
            pass
        # and close the shelve file
        s.close()

    # controls visibility of dropdown menus based on user input
    def hide1(self):
        global dropdown1_open
        dropdown1_open = True
        global dropdown2_open
        dropdown2_open = False
        global dropdown3_open
        dropdown3_open = False
        index = self.chooseFunction1.findText(self.chooseFunction1.currentText())
        index2 = self.chooseFunction2.findText(self.chooseFunction2.currentText())
        #print "You have selected %s " %self.chooseFunction1.currentText()
        #print "You have selected %s " %self.chooseFunction2.currentText()
        if( index != 7 or index2 != 6):
            self.chooseFunction3.setVisible(False)
            dropdown3_open = False
        if ( index != -1 and index != 0 and index != 1 and index != 2 and index != 3 and index != 4 and index!= 8):
            #print(index)
            self.chooseFunction2.setVisible(True)
            dropdown2_open = True
            if(index == 5):
                self.chooseFunction2.setCurrentIndex(0)
                self.chooseFunction2.model().item(1).setEnabled(True)
                self.chooseFunction2.model().item(2).setEnabled(True)
                self.chooseFunction2.model().item(3).setEnabled(True)
                self.chooseFunction2.model().item(4).setEnabled(False)
                self.chooseFunction2.model().item(5).setEnabled(False)
                self.chooseFunction2.model().item(6).setEnabled(False)
                self.chooseFunction2.model().item(7).setEnabled(False)
            if(index == 6):
                self.chooseFunction2.setCurrentIndex(0)
                self.chooseFunction2.model().item(1).setEnabled(False)
                self.chooseFunction2.model().item(2).setEnabled(True)
                self.chooseFunction2.model().item(3).setEnabled(False)
                self.chooseFunction2.model().item(4).setEnabled(True)
                self.chooseFunction2.model().item(5).setEnabled(True)
                self.chooseFunction2.model().item(6).setEnabled(False)
                self.chooseFunction2.model().item(7).setEnabled(False)
            if(index == 7):
                self.chooseFunction2.setCurrentIndex(0)
                self.chooseFunction2.model().item(1).setEnabled(False)
                self.chooseFunction2.model().item(2).setEnabled(False)
                self.chooseFunction2.model().item(3).setEnabled(False)
                self.chooseFunction2.model().item(4).setEnabled(False)
                self.chooseFunction2.model().item(5).setEnabled(False)
                self.chooseFunction2.model().item(6).setEnabled(True)
                self.chooseFunction2.model().item(7).setEnabled(True)
        else:
            self.chooseFunction2.setVisible(False)
            self.chooseFunction3.setVisible(False)
            dropdown2_open = False
            dropdown3_open = False

    # changes visibility of the third function selection dropdown menu
    def hide2(self):
        global dropdown1_open
        dropdown1_open = True
        global dropdown2_open
        dropdown2_open = True
        global dropdown3_open
        dropdown3_open = False
        index1 = self.chooseFunction1.findText(self.chooseFunction1.currentText())
        index2 = self.chooseFunction2.findText(self.chooseFunction2.currentText())
        if ( index1 == 7 and index2 == 6):
            #print "You have selected %s " %self.chooseFunction1.currentText()
            #print "You have selected %s " %self.chooseFunction2.currentText()
            self.chooseFunction3.setVisible(True)
            dropdown3_open = True
        else:
            #print "You have selected %s " %self.chooseFunction1.currentText()
            #print "You have selected %s " %self.chooseFunction2.currentText()
            self.chooseFunction3.setVisible(False)
            dropdown3_open = False

    # want to display dropdown menus depending on the choice
    def chooseFunction(self):
        print()
        # if the function has been chosen from the dropdown menu,
        # assign a var to the text
        if dropdown1_open == True:
            if dropdown2_open == True and self.chooseFunction2.findText(self.chooseFunction2.currentText()) != 0:
                if dropdown3_open == True and self.chooseFunction3.findText(self.chooseFunction3.currentText()) != 0:
                    #all three dropdown menus are open

                    #checks if the first dropdown menu contains an index parameter
                    #splits the string two string vars,
                    # one containing the function name, other
                    # containing the parameter
                    if '[' in self.chooseFunction1.currentText() and ']' in self.chooseFunction1.currentText():
                        function1, dummyVar1 = self.chooseFunction1.currentText().split('[')
                        param1, dummyVar2 = dummyVar1.split(']')
                        function1 = str(function1)
                        if(param1 == "index"):
                            print("Insert proper index")
                            return
                        #function2 = function_mappings[ self.chooseFunction2.currentText() ]
                        #function3 = function_mappings[ self.chooseFunction3.currentText() ]
                        function2 = str(self.chooseFunction2.currentText())
                        function3 = str(self.chooseFunction3.currentText())
                        #model.function1[int(param1)].function2.function3
                        #print '3 Dropdown menu with index parameters'
                        #print function1
                        #print param1
                        #print function2
                        #print function3
                        a = getattr(model, function1)[int(param1)]
                        print(getattr( getattr( a , function2) , function3))
                    '''
                    else:
                        # three functions chained to model
                        #function1 = function_mappings[ self.chooseFunction1.currentText() ]
                        #function2 = function_mappings[ self.chooseFunction2.currentText() ]
                        #function3 = function_mappings[ self.chooseFunction3.currentText() ]
                        function1 = str(self.chooseFunction1.currentText())
                        function2 = str(self.chooseFunction2.currentText())
                        function3 = str(self.chooseFunction3.currentText())
                        #model.function1.function2.function3
                        print '3 Dropdown menu without index parameters'
                        print function1
                        print function2
                        print function3
                        a = getattr(model, function1)[]
                        print getattr( getattr( getattr(model, function1)() , function2)() , function3)()
                    '''

                    #pseudo code ends

                #otherwise, dropdownmenu3 is closed, 1 and 2 are open
                else:
                    # two functions chained to model
                    if '[' in self.chooseFunction1.currentText() and ']' in self.chooseFunction1.currentText():
                        function1, dummyVar1 = self.chooseFunction1.currentText().split('[')
                        param1, dummyVar2 = dummyVar1.split(']')
                        function1 = str(function1)

                        #function2 = function_mappings[ self.chooseFunction2.currentText() ]
                        function2 =  str(self.chooseFunction2.currentText())
                        #print function2
                        #print str(function2)
                        #model.function1.function2
                        if(param1 == "index"):
                            print("Insert proper index")
                            return
                        param1 = int(param1)
                        #print '2 Dropdown menu with index parameters'
                        #print function1
                        #print param1
                        #print function2
                        a = getattr(model, function1)[param1]
                        #print a
                        print(getattr( a , function2))
                    '''
                    else:
                        # three functions chained to model
                        #function1 = function_mappings[ self.chooseFunction1.currentText() ]
                        #function2 = function_mappings[ self.chooseFunction2.currentText() ]
                        function1 = str(self.chooseFunction1.currentText())
                        function2 = str(self.chooseFunction2.currentText())

                        #model.function1.function2
                        print '2 Dropdown menu without index parameters'
                        print function1
                        print param1
                        print function2
                        print getattr( getattr(model, function1)() , function2)()
                    '''


            #otherwise, dropdownmenu2 is closed but dropdownmenu1 is open
            else:
                # contains an index parameter
                if '[' in self.chooseFunction1.currentText() and ']' in self.chooseFunction1.currentText():
                    function1, dummyVar1 = self.chooseFunction1.currentText().split('[')
                    param1, dummyVar2 = dummyVar1.split(']')
                    function1 = str(function1)
                    param1 = int(param1)
                    #model.function1
                    if(param1 == "index"):
                        print("Insert proper index")
                        return
                    #print '1 Dropdown menu with index parameters'
                    #print function1
                    #print param1
                    print(getattr(model, function1)(param1))

                #does not contain an index parameter
                else:
                    # three functions chained to model
                    #function1 = function_mappings[ self.chooseFunction1.currentText() ]
                    function1 = str(self.chooseFunction1.currentText())
                    #model.function1
                    #print '1 Dropdown menu without index parameters'
                    #print function1
                    print(getattr(model, function1)())

        #otherwise, dropdownmenu1 is closed
        else:
            print('You have not selected a function')


# need to add all below to display gui
def main(): # defines main function
    app = QtGui.QApplication(sys.argv)
    MainWindow = QtGui.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

main()
