"""Qt compatibility helpers for chromatoPy.

The desktop wrapper prefers PySide6, but several existing modules were written
against PyQt5. This module exposes a small common surface so the rest of the
package can stay focused on workflow logic.
"""

from __future__ import annotations

QT_API = ""

# try:
#     from PySide6 import QtCore, QtGui, QtWidgets
#
#     QT_API = "PySide6"
#     Signal = QtCore.Signal
# except ImportError:  # pragma: no cover - fallback for older environments
#     from PyQt5 import QtCore, QtGui, QtWidgets
#
#     QT_API = "PyQt5"
#     Signal = QtCore.pyqtSignal
from PySide6 import QtCore, QtGui, QtWidgets


Qt = QtCore.Qt
ApplicationModal = getattr(Qt, "ApplicationModal", Qt.WindowModality.ApplicationModal)
StrongFocus = getattr(Qt, "StrongFocus", Qt.FocusPolicy.StrongFocus)
ShiftModifier = getattr(Qt, "ShiftModifier", Qt.KeyboardModifier.ShiftModifier)
Key_Backspace = getattr(Qt, "Key_Backspace", Qt.Key.Key_Backspace)
Key_Delete = getattr(Qt, "Key_Delete", Qt.Key.Key_Delete)
WaitCursor = getattr(Qt, "WaitCursor", Qt.CursorShape.WaitCursor)

QApplication = QtWidgets.QApplication
QDialog = QtWidgets.QDialog
QGroupBox = QtWidgets.QGroupBox
QWidget = QtWidgets.QWidget
QMainWindow = QtWidgets.QMainWindow

QLabel = QtWidgets.QLabel
QLineEdit = QtWidgets.QLineEdit
QPushButton = QtWidgets.QPushButton
QCheckBox = QtWidgets.QCheckBox
QSpinBox = QtWidgets.QSpinBox
QDoubleSpinBox = QtWidgets.QDoubleSpinBox
QTextEdit = QtWidgets.QTextEdit
QPlainTextEdit = QtWidgets.QPlainTextEdit
QListWidget = QtWidgets.QListWidget
QListWidgetItem = QtWidgets.QListWidgetItem
QStackedWidget = QtWidgets.QStackedWidget
QFrame = QtWidgets.QFrame
QScrollArea = QtWidgets.QScrollArea
QFileDialog = QtWidgets.QFileDialog
QComboBox = QtWidgets.QComboBox
QFormLayout = QtWidgets.QFormLayout
QVBoxLayout = QtWidgets.QVBoxLayout
QHBoxLayout = QtWidgets.QHBoxLayout
QGridLayout = QtWidgets.QGridLayout
QMessageBox = QtWidgets.QMessageBox
QDialogButtonBox = QtWidgets.QDialogButtonBox
QSplitter = QtWidgets.QSplitter
QSizePolicy = QtWidgets.QSizePolicy
QEventLoop = QtCore.QEventLoop
QTimer = QtCore.QTimer


def exec_dialog(dialog: QDialog) -> int:
    """Run a modal dialog across Qt bindings."""

    if hasattr(dialog, "exec"):
        return dialog.exec()
    return dialog.exec_()


def run_application(app: QApplication) -> int:
    """Start the Qt event loop across Qt bindings."""

    if hasattr(app, "exec"):
        return app.exec()
    return app.exec_()
