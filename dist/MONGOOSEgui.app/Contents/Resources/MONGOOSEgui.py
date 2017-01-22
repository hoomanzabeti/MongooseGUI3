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
        self.chooseFunction1.setGeometry(QtCore.QRect(20, 137, 145, 50))
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

        # intializes index input
        self.chooseIndex = QtGui.QLineEdit(self.centralwidget)
        self.chooseIndex.setGeometry(QtCore.QRect(165, 140, 53, 45))
        self.chooseIndex.setObjectName(_fromUtf8("chooseIndex"))
        self.chooseIndex.setPlaceholderText('[index]')
        self.chooseIndex.setVisible(False)

        # intializes third dropdown menu, chooses a function
        self.chooseFunction3 = QtGui.QComboBox(self.centralwidget)
        self.chooseFunction3.setGeometry(QtCore.QRect(20, 198, 200, 50))
        self.chooseFunction3.setObjectName(_fromUtf8("chooseFunction3"))
        self.chooseFunction3.addItem(_fromUtf8(""))
        self.chooseFunction3.addItem(_fromUtf8(""))
        self.chooseFunction3.setVisible(False)

        # calls MONGOOSE functions
        self.executeAction.clicked.connect(self.chooseFunction)

        self.textEdit1 = QtGui.QTextEdit(self.centralwidget)
        self.textEdit1.setGeometry(QtCore.QRect(100, 240, 115, 25))
        self.textEdit1.setObjectName(_fromUtf8("textEdit"))
        self.label1 = QtGui.QLabel(self.centralwidget)
        self.label1.setGeometry(QtCore.QRect(30, 245, 56, 13))
        self.label1.setObjectName(_fromUtf8("label"))
        self.textEdit1.setVisible(False)
        self.label1.setVisible(False)



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
        self.chooseFunction1.setItemText(5, _translate("MainWindow", "reactions", None)) # []
        self.chooseFunction1.setItemText(6, _translate("MainWindow", "reactionSubsets", None)) # []
        self.chooseFunction1.setItemText(7, _translate("MainWindow", "metabolites", None))# []
        self.chooseFunction1.setItemText(8, _translate("MainWindow", "printReactionFormula", None))# []
        self.chooseFunction3.setItemText(0, _translate("MainWindow", "<Choose>", None))
        self.chooseFunction3.setItemText(1, _translate("MainWindow", "name", None))
        self.label1.setText(_translate("MainWindow", "Var Name", None))

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
          #boolean variables to flag if dropdown menus are visible
          global dropdown1_open
          dropdown1_open = True
          global dropdown2_open
          dropdown2_open = False
          global dropdown3_open
          dropdown3_open = False

          #fetches current index of dropdown menu
          index = self.chooseFunction1.findText(self.chooseFunction1.currentText())
          index2 = self.chooseFunction2.findText(self.chooseFunction2.currentText())

          #initializes constant values corresponding to index value of dropdown menu
          ERROR = -1
          CHOOSE = 0

          #first dropdown menu
          REDUCE_NETWORK = 1
          ADD_RXN = 2
          FIND_ESSENTIAL_RXN = 3
          FIND_BIOMASS_RXN = 4
          REACTIONS = 5
          REACTION_SUBSETS = 6
          METABOLITES = 7
          PRINT_RXN_FORMULA = 8

          #second dropdown menu
          NAME = 1
          PAIRS = 2
          REDUCTION_STATUS = 3
          REVERSIBLE = 4
          LENGTH = 5
          SPECIES = 6
          EXTERNAL = 7

          if( index != METABOLITES or index2 != SPECIES):
              self.chooseFunction3.setVisible(False)
              dropdown3_open = False
          if ( index != ERROR and index != CHOOSE and index != REDUCE_NETWORK and
                index != ADD_RXN and index != FIND_ESSENTIAL_RXN and
                index != FIND_BIOMASS_RXN and index!= PRINT_RXN_FORMULA):
              #print(index)
              self.chooseFunction2.setVisible(True)
              dropdown2_open = True
              #enables which options from dropdown menu are selectable
              if(index == REACTIONS):
                  self.chooseIndex.setVisible(True)
                  self.chooseFunction2.setCurrentIndex(CHOOSE)
                  self.chooseFunction2.model().item(NAME).setEnabled(True)
                  self.chooseFunction2.model().item(PAIRS).setEnabled(True)
                  self.chooseFunction2.model().item(REDUCTION_STATUS).setEnabled(True)
                  self.chooseFunction2.model().item(REVERSIBLE).setEnabled(False)
                  self.chooseFunction2.model().item(LENGTH).setEnabled(False)
                  self.chooseFunction2.model().item(SPECIES).setEnabled(False)
                  self.chooseFunction2.model().item(EXTERNAL).setEnabled(False)
              if(index == REACTION_SUBSETS):
                  self.chooseIndex.setVisible(True)
                  self.chooseFunction2.setCurrentIndex(CHOOSE)
                  self.chooseFunction2.model().item(NAME).setEnabled(False)
                  self.chooseFunction2.model().item(PAIRS).setEnabled(True)
                  self.chooseFunction2.model().item(REDUCTION_STATUS).setEnabled(False)
                  self.chooseFunction2.model().item(REVERSIBLE).setEnabled(True)
                  self.chooseFunction2.model().item(LENGTH).setEnabled(True)
                  self.chooseFunction2.model().item(SPECIES).setEnabled(False)
                  self.chooseFunction2.model().item(EXTERNAL).setEnabled(False)
              if(index == METABOLITES):
                  self.chooseIndex.setVisible(True)
                  self.chooseFunction2.setCurrentIndex(CHOOSE)
                  self.chooseFunction2.model().item(NAME).setEnabled(False)
                  self.chooseFunction2.model().item(PAIRS).setEnabled(False)
                  self.chooseFunction2.model().item(REDUCTION_STATUS).setEnabled(False)
                  self.chooseFunction2.model().item(REVERSIBLE).setEnabled(False)
                  self.chooseFunction2.model().item(LENGTH).setEnabled(False)
                  self.chooseFunction2.model().item(SPECIES).setEnabled(True)
                  self.chooseFunction2.model().item(EXTERNAL).setEnabled(True)
          else:
              if(index == PRINT_RXN_FORMULA):
                  self.chooseIndex.setVisible(True)
              else:
                  self.chooseIndex.setVisible(False)
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

        ERROR = -1
        CHOOSE = 0

        REDUCE_NETWORK = 1
        ADD_RXN = 2
        FIND_ESSENTIAL_RXN = 3
        FIND_BIOMASS_RXN = 4
        REACTIONS = 5
        REACTION_SUBSETS = 6
        METABOLITES = 7
        PRINT_RXN_FORMULA = 8

        NAME = 1
        PAIRS = 2
        REDUCTION_STATUS = 3
        REVERSIBLE = 4
        LENGTH = 5
        SPECIES = 6
        EXTERNAL = 7

        if ( index1 == METABOLITES and index2 == SPECIES):
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
        #declares the options for corresponding index
        CHOOSE = 0
        REACTIONS = 5
        REACTION_SUBSETS = 6
        METABOLITES = 7
        PRINT_RXN_FORMULA = 8
        #grabs the index of the current dropdown option
        index1 = self.chooseFunction1.findText(self.chooseFunction1.currentText())
        # if the function has been chosen from the dropdown menu,
        # assign a var to the text
        if dropdown1_open == True:
            if dropdown2_open == True and self.chooseFunction2.findText(self.chooseFunction2.currentText()) != CHOOSE:
                if dropdown3_open == True and self.chooseFunction3.findText(self.chooseFunction3.currentText()) != CHOOSE:
                    #all three functions
                    function1 = str(self.chooseFunction1.currentText())
                    function2 = str(self.chooseFunction2.currentText())
                    function3 = str(self.chooseFunction3.currentText())
                    param1 = self.chooseIndex.text()
                    if param1:
                        param1 = int ( str( param1 ) )
                    else:
                        print("Choose an index")

                    #call function and display output
                    chain1 = getattr(model, function1)[param1]
                    chain2 = getattr(chain1, function2)
                    print(getattr(chain2, function3))

                else: # only two fuctions
                    function1 = str(self.chooseFunction1.currentText())
                    function2 = str(self.chooseFunction2.currentText())
                    param1 = self.chooseIndex.text()
                    if param1:
                        param1 = int ( str( param1 ) )
                    else:
                        print("Choose an index")
                    #call function and display output
                    chain1 = getattr(model, function1)[param1]
                    print(getattr(chain1, function2))

            else: # only one function

                if( index1 == REACTIONS or index1 == REACTION_SUBSETS or index1 == METABOLITES or index1 == PRINT_RXN_FORMULA):
                    function1 = str(self.chooseFunction1.currentText())
                    param1 = self.chooseIndex.text()
                    if param1:
                        param1 = int ( str( param1 ) )
                    else:
                        print("Choose an index")

                    if( index1 != PRINT_RXN_FORMULA):
                        #call function and display output
                        print(getattr(model, function1)[param1])
                    else:
                        #call function and display output
                        print(getattr(model, function1)(param1))

                else:
                    function1 = str(self.chooseFunction1.currentText())
                    #call function and display output
                    print(getattr(model, function1)())


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
