'''
Coding Standard:
----------------
-> lower/upper case for multi word variables
-> all upper case for constants


Most Recent Changes:
--------------------

1. Changes addMetabolite: (line 661)
addMetabolite now takes the following parameters: (string, character)
2. Fixed error that would display array index when not necessary

'''
import os
import sys
import logging
import dbm
import shelve
import warnings
from functools import partial
from MetaMerge import *
from PyQt4 import QtCore, QtGui, uic
from PyQt4.QtCore import QThread


with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    warnings.simplefilter("ignore")


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

# redirects console output to GUI
class QtHandler(logging.Handler):
    def __init__(self):
        logging.Handler.__init__(self)
    def emit(self, record):
        record = self.format(record)
        if record: XStream.stdout().write('%s\n'%record)
        # originally: XStream.stdout().write("{}\n".format(record))
    #ignore warnings from terminal
    def handler(msg_type, msg_log_context, msg_string):
        pass

# redirects console output to GUI
logger = logging.getLogger(__name__)
handler = QtHandler()
handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
logger.addHandler(handler)
logger.setLevel(logging.CRITICAL)

# redirects console output to GUI
class XStream(QtCore.QObject):
    _stdout = None
    _stderr = None
    process = QtCore.QProcess()
    messageWritten = QtCore.pyqtSignal(str)
    #messageWritten = str(process.readAllStandardOutput())
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

# provided multi threading capabilities
class WorkerThread(QThread):

    def __init__(self,model_name, function_name, index_name):
        QThread.__init__(self)
        self.model_name = model_name
        self.function_name = function_name
        self.index_name = index_name
        self.reduceNetwork = 1
        self.findSyntheticLethalPairs = 19

    def __del__(self):
        self.wait()


    def run(self):
        if(self.index_name == self.reduceNetwork):
            print(getattr(self.model_name, self.function_name)())
            print("Finished reducing network")
        elif(self.index_name == self.findSyntheticLethalPairs):
            print(getattr(self.model_name, self.function_name)())
            print("Finished finding synthetic lethal pairs")
        else:
            print("Thread could not find work")
            print("indexes name: %s" %(self.index_name))

# defines UI
class Ui_MainWindow(object):
    def setupUi(self, MainWindow):

        # size constraints
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        #MainWindow.resize(600, 470)
        MainWindow.resize(650, 470)

        # central widget ( main window for buttons )
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))

        # output
        self.resultOutput = QtGui.QTextEdit(self.centralwidget)
        #self.resultOutput.setGeometry(QtCore.QRect(230, 50, 361, 370))
        self.resultOutput.setGeometry(QtCore.QRect(280, 50, 361, 390))
        self.resultOutput.setObjectName(_fromUtf8("resultOutput"))

        #disables the output box to be changed if user clicks, new lines printed out at bottom
        self.resultOutput.setReadOnly(True)
        self.resultOutput.ensureCursorVisible()
        self.resultOutput.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        cursor = self.resultOutput.textCursor()
        cursor.movePosition(QtGui.QTextCursor.EndOfBlock, QtGui.QTextCursor.MoveAnchor)
        self.resultOutput.setTextCursor(cursor)
        self.resultOutput.moveCursor(QtGui.QTextCursor.End)
        self.resultOutput.verticalScrollBar().setValue(self.resultOutput.verticalScrollBar().maximum())

        #XStream.stdout().messageWritten.connect( self.resultOutput.insertPlainText )
        #XStream.stderr().messageWritten.connect( self.resultOutput.insertPlainText )
        XStream.stdout().messageWritten.connect( self.resultOutput.append )
        #XStream.stderr().messageWritten.connect( self.resultOutput.append )
        self.resultOutput.setTextInteractionFlags(QtCore.Qt.NoTextInteraction) #crucial line

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

        # intializes writeSBML button
        self.writeSBML = QtGui.QPushButton(self.centralwidget)
        self.writeSBML.setGeometry(QtCore.QRect(20, 420, 200, 32))
        self.writeSBML.setObjectName(_fromUtf8("writeSBML"))
        self.writeSBML.clicked.connect(self.write_sbml)
        self.writeSBML.setEnabled(False)


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
        #self.chooseFunction1.setGeometry(QtCore.QRect(20, 137, 145, 50))
        self.chooseFunction1.setGeometry(QtCore.QRect(20, 137, 200, 50))
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
        self.chooseFunction1.addItem(_fromUtf8(""))
        self.chooseFunction1.addItem(_fromUtf8(""))
        self.chooseFunction1.addItem(_fromUtf8(""))
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
        self.chooseIndex.setGeometry(QtCore.QRect(220, 140, 53, 45))
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

        self.rxnParam1 = QtGui.QLineEdit(self.centralwidget)
        self.rxnParam1.setGeometry(QtCore.QRect(20, 240, 190, 25))
        self.rxnParam1.setObjectName(_fromUtf8("rxnParam1"))
        self.label1 = QtGui.QLabel(self.centralwidget)
        self.label1.setGeometry(QtCore.QRect(30, 245, 56, 13))
        self.label1.setObjectName(_fromUtf8("label1"))
        self.rxnParam1.setPlaceholderText('Reaction Name')
        self.rxnParam1.setVisible(False)
        self.label1.setVisible(False)

        self.rxnParam2 = QtGui.QLineEdit(self.centralwidget)
        self.rxnParam2.setGeometry(QtCore.QRect(20, 270, 190, 25))
        self.rxnParam2.setObjectName(_fromUtf8("rxnParam2"))
        self.label2 = QtGui.QLabel(self.centralwidget)
        self.label2.setGeometry(QtCore.QRect(30, 275, 56, 13))
        self.label2.setObjectName(_fromUtf8("label2"))
        self.rxnParam2.setPlaceholderText('List Of Pairs')
        self.rxnParam2.setVisible(False)
        self.label2.setVisible(False)

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
        self.writeSBML.setText(_translate("MainWindow", "Write SBML", None))
        self.selectModel.setText(_translate("MainWindow", "Enter", None))
        self.title.setText(_translate("MainWindow", "MONGOOSE Exact Arithmetic Toolbox", None))
        self.title.setToolTip(_translate("MainWindow", "<html><head/><body><p> <b>MONGOOSE Author:</b><br> Leonid Chindelevitch <br><br> <b>GUI Author:</b  > <br> Christopher Le </p></body></html>", None))
        self.chooseFunction2.setItemText(0, _translate("MainWindow", "Choose", None))
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
        self.chooseFunction1.setItemText(2, _translate("MainWindow", "addReaction", None))#name and list of pairs
        self.chooseFunction1.setItemText(3, _translate("MainWindow", "deleteReactions", None))
        self.chooseFunction1.setItemText(4, _translate("MainWindow", "findBiomassReaction", None))
        self.chooseFunction1.setItemText(5, _translate("MainWindow", "reactions", None)) # [index]
        self.chooseFunction1.setItemText(6, _translate("MainWindow", "reactionSubsets", None)) # [index]
        self.chooseFunction1.setItemText(7, _translate("MainWindow", "metabolites", None))# [index]
        self.chooseFunction1.setItemText(8, _translate("MainWindow", "printReactionFormula", None))# [index]
        self.chooseFunction1.setItemText(9, _translate("MainWindow", "deleteMetabolites", None))# [index]
        self.chooseFunction1.setItemText(10, _translate("MainWindow", "addMetabolite", None))# name and a single pair
        self.chooseFunction1.setItemText(11, _translate("MainWindow", "checkElementalBalance", None))# new
        self.chooseFunction1.setItemText(12, _translate("MainWindow", "findTopoBlockedMetabolites", None))# new
        self.chooseFunction1.setItemText(13, _translate("MainWindow", "findTopoBlockedReactions", None))# new
        self.chooseFunction1.setItemText(14, _translate("MainWindow", "findStoichBlockedReactions", None))# new
        self.chooseFunction1.setItemText(15, _translate("MainWindow", "findIrrevBlockedReactions", None))# new
        self.chooseFunction1.setItemText(16, _translate("MainWindow", "findSemiBlockedReactions", None))# new
        self.chooseFunction1.setItemText(17, _translate("MainWindow", "unblockBiomassReaction", None))# new
        self.chooseFunction1.setItemText(18, _translate("MainWindow", "findEssentialReactions", None))# new
        self.chooseFunction1.setItemText(19, _translate("MainWindow", "findSyntheticLethalPairs", None))# new
        self.chooseFunction1.setItemText(20, _translate("MainWindow", "findMinimalMedia", None))# new
        self.chooseFunction3.setItemText(0, _translate("MainWindow", "<Choose>", None))
        self.chooseFunction3.setItemText(1, _translate("MainWindow", "name", None))
        self.label1.setText(_translate("MainWindow", "Name", None))
        self.label2.setText(_translate("MainWindow", "Pairs", None))



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
            self.writeSBML.setEnabled(True)
            self.selectModel.setEnabled(False)
            self.selectFile.setEnabled(False)
            self._running = False
        else:
            print('Not a valid option, please choose again')
            self._running = True



    def open (self):
        #extracts the desired file, and parses it into proper 3.5 db format
        filename_temp = QtGui.QFileDialog.getOpenFileName(self.centralwidget, 'Open File', '.')
        global s
        filename_temp = str(filename_temp)
        filename = filename_temp.split('.')
        s = shelve.open(str(filename[0])) #3.5

        #prompts user for the desired model
        print('Please type which of the following models you would like to analyze: ')
        print((list(s.keys())))
        self._running = True
        #loops until select_model button has been clicked
        while self._running:
            QtGui.qApp.processEvents()
        global desiredModel
        desiredModel = str(modelInput)
        global model
        model = s[desiredModel]

    # saves contents of console to desired text file
    def save_content (self):
        choice = QtGui.QMessageBox.question(QtGui.QMainWindow(), 'Save Contents', 'Are you sure you want to save the contents?', QtGui.QMessageBox.No | QtGui.QMessageBox.Yes)
        #choice = QtGui.QMessageBox.question(QtGui.QMainWindow(), 'Save Contents', 'Are you sure you want to save the contents?', QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
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
        #choice = QtGui.QMessageBox.question(QtGui.QMainWindow(), 'Save Model', 'Are you sure you want to save the model?', QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        choice = QtGui.QMessageBox.question(QtGui.QMainWindow(), 'Save Model', 'Are you sure you want to save the model?', QtGui.QMessageBox.No | QtGui.QMessageBox.Yes)
        if choice == QtGui.QMessageBox.Yes:
            savedModel = "'" + modelInput + "'"
            s[desiredModel] = model
            modelSaved = True
            print(("Saving model: " + modelInput))
        else:
            pass

        choice = QtGui.QMessageBox.question(QtGui.QMainWindow(), 'Save & Quit', 'Would you like to exit the application?', QtGui.QMessageBox.No | QtGui.QMessageBox.Yes)
        #choice = QtGui.QMessageBox.question(QtGui.QMainWindow(), 'Save & Quit', 'Would you like to exit the application?', QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        if choice == QtGui.QMessageBox.Yes:
            print("Closing shelve file")
            # close the shelve file
            s.close()
            sys.exit()
        else:
            pass

    def write_sbml(self):
        choice = QtGui.QMessageBox.question(QtGui.QMainWindow(), 'Write SBML', 'Are you sure you want to write to SBML?', QtGui.QMessageBox.No | QtGui.QMessageBox.Yes)
        #choice = QtGui.QMessageBox.question(QtGui.QMainWindow(), 'Save Contents', 'Are you sure you want to save the contents?', QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        if choice == QtGui.QMessageBox.Yes:
            name = QtGui.QFileDialog.getSaveFileName(QtGui.QMainWindow(), 'Save File')
            model.writeSBML(name)
        else:
            pass



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
          DELETE_RXN = 3
          FIND_BIOMASS_RXN = 4
          REACTIONS = 5
          REACTION_SUBSETS = 6
          METABOLITES = 7
          PRINT_RXN_FORMULA = 8
          DELETE_METAB = 9
          ADD_METAB = 10
          CHECK_ELEM_BAL = 11 #checkElementalBalance
          FIND_TOP_BLK_METAB = 12 #findTopoBlockedMetabolites
          FIND_TOP_BLK_RXN = 13 #findTopoBlockedReactions
          FIND_STOICH_BLK_RXN = 14 #findStoichBlockedReactions
          FIND_IRREV_BLK_RXN = 15 #findIrrevBlockedReactions
          FIND_SEMI_BLK_RXN = 16 #findSemiBlockedReactions
          UNBLK_BIOM_RXN = 17 #unblockBiomassReaction
          FIND_ESS_RXN = 18 #findEssentialReactions
          FIND_SYNTH_LETH_PAIRS = 19 #findSyntheticLethalPairs
          FIND_MIN_MED = 20 #findMinimalMedia

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
              self.rxnParam1.setVisible(False)
              self.rxnParam2.setVisible(False)
              dropdown3_open = False
          if ( index != ERROR and index != CHOOSE and index != REDUCE_NETWORK and
                index != ADD_RXN and index != DELETE_RXN and
                index != FIND_BIOMASS_RXN and index!= PRINT_RXN_FORMULA and index!=DELETE_RXN
                and index!=DELETE_METAB and index!=ADD_METAB and index!=CHECK_ELEM_BAL and index!=FIND_TOP_BLK_METAB
                 and index!=FIND_TOP_BLK_RXN and index!=FIND_STOICH_BLK_RXN and index!=FIND_IRREV_BLK_RXN and index!=FIND_SEMI_BLK_RXN
                  and index!=UNBLK_BIOM_RXN and index!=FIND_ESS_RXN and index!=FIND_SYNTH_LETH_PAIRS and index!=FIND_MIN_MED):
              #print(index)
              self.chooseFunction2.setVisible(True)
              self.rxnParam1.setVisible(False)
              self.rxnParam2.setVisible(False)
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
            if(index == PRINT_RXN_FORMULA or index==DELETE_RXN or index==DELETE_METAB):
                self.chooseIndex.setVisible(True)

            elif(index == ADD_RXN or index == ADD_METAB):
                self.rxnParam1.setVisible(True)
                self.rxnParam2.setVisible(True)
                self.chooseIndex.setVisible(False)

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
        DELETE_RXN = 3
        FIND_BIOMASS_RXN = 4
        REACTIONS = 5
        REACTION_SUBSETS = 6
        METABOLITES = 7
        PRINT_RXN_FORMULA = 8
        DELETE_METAB = 9
        ADD_METAB = 10

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
        #print

        #declares the options for corresponding index
        ERROR = -1
        CHOOSE = 0

        REDUCE_NETWORK = 1
        ADD_RXN = 2
        DELETE_RXN = 3
        FIND_BIOMASS_RXN = 4
        REACTIONS = 5
        REACTION_SUBSETS = 6
        METABOLITES = 7
        PRINT_RXN_FORMULA = 8
        DELETE_METAB = 9
        ADD_METAB = 10
        CHECK_ELEM_BAL = 11 #checkElementalBalance
        FIND_TOP_BLK_METAB = 12 #findTopoBlockedMetabolites
        FIND_TOP_BLK_RXN = 13 #findTopoBlockedReactions
        FIND_STOICH_BLK_RXN = 14 #findStoichBlockedReactions
        FIND_IRREV_BLK_RXN = 15 #findIrrevBlockedReactions
        FIND_SEMI_BLK_RXN = 16 #findSemiBlockedReactions
        UNBLK_BIOM_RXN = 17 #unblockBiomassReaction
        FIND_ESS_RXN = 18 #findEssentialReactions
        FIND_SYNTH_LETH_PAIRS = 19 #findSyntheticLethalPairs
        FIND_MIN_MED = 20 #findMinimalMedia

        NAME = 1
        PAIRS = 2
        REDUCTION_STATUS = 3
        REVERSIBLE = 4
        LENGTH = 5
        SPECIES = 6
        EXTERNAL = 7


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
                    print(">>> model.%s[%s].%s.%s" % (function1,param1,function2,function3))
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
                    print(">>> model.%s[%s].%s" % (function1,param1,function2))
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
                        print(">>> model.%s(%d)" % (function1, param1))
                        print(getattr(model, function1)[param1])
                    else:
                        #call function and display output
                        print(">>> model.%s(%d)" % (function1, param1))
                        print(getattr(model, function1)(param1))

                elif(index1 == ADD_RXN):
                    function1 = str(self.chooseFunction1.currentText())
                    addReactionName = str(self.rxnParam1.text())
                    addReactionPairs = str(self.rxnParam2.text())
                    addReactionPairs = addReactionPairs.replace('[', "")
                    addReactionPairs = addReactionPairs.replace(']', "")
                    count = 0
                    if(count==0):
                        list = [[int(el) for el in addReactionPairs.split(',')]]
                        count = count + 1
                    else:
                        list = [[Fraction(el) for el in addReactionPairs.split(',')]]
                    print(">>> model.%s('%s',%s)" % (function1,addReactionName,list))
                    print(getattr(model, function1)(addReactionName, list))

                elif(index1 == ADD_METAB):
                    function1 = str(self.chooseFunction1.currentText())
                    addMetaboliteName = str(self.rxnParam1.text())
                    addMetaboliteCompartment = str(self.rxnParam2.text())
                    print(">>> model.%s('%s',%s)" % (function1,addMetaboliteName,addMetaboliteCompartment))
                    print(getattr(model, function1)(addMetaboliteName, addMetaboliteCompartment))

                elif(index1 == DELETE_METAB or index1 == DELETE_RXN):
                    function1 = str(self.chooseFunction1.currentText())
                    param1 = self.chooseIndex.text()
                    if param1:
                        param1 = str( param1 )
                    else:
                        print("Choose an index")
                    count = 0
                    if( count == 0):
                        list = [int(el) for el in param1.split(',')]
                        count = count + 1
                    else:
                        list = [Fraction(el) for el in param1.split(',')]
                    print(">>> model.%s(%s)" % (function1,list))
                    print(getattr(model, function1)(list))
                else:

                    function1 = str(self.chooseFunction1.currentText())
                    #call function and display output
                    print(">>> model.%s()" % (function1))
                    if(index1 == REDUCE_NETWORK):
                        print("Reducing network - this may take some time!")
                        self.myThread = WorkerThread(model,function1, REDUCE_NETWORK)
                        self.myThread.start()
                    elif(index1 == FIND_SYNTH_LETH_PAIRS):
                        print("Finding synthetic lethal pairs - this may take some time!")
                        self.myThread = WorkerThread(model,function1, FIND_SYNTH_LETH_PAIRS)
                        self.myThread.start()
                    else:
                        print(getattr(model, function1)())


        else:
            print('You have not selected a function')


#handles errors (specifically unwanted QNSView mouse dragged warning)
def handler(msg_type, msg_string):
    pass
QtCore.qInstallMsgHandler(handler)



# need to add all below to display gui
def main(): # defines main function
    app = QtGui.QApplication(sys.argv)
    app.processEvents()
    MainWindow = QtGui.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())


main()
