"""Desktop application entry point."""

from __future__ import annotations

import os
import sys

os.environ.setdefault("QT_API", "pyside6")

import matplotlib

matplotlib.use("QtAgg")

from ..qt_compat import QApplication, run_application
from .views import ChromatoPyMainWindow


APP_STYLE = """
QWidget {
    background: #f5efe2;
    color: #2f3b43;
    font-size: 13px;
}
QLabel#appTitle {
    font-size: 30px;
    font-weight: 700;
    color: #24313a;
}
QLabel#appSubtitle {
    font-size: 14px;
    color: #51616d;
}
QLabel#pageTitle {
    font-size: 22px;
    font-weight: 700;
    color: #24313a;
}
QLabel#pageDescription {
    color: #51616d;
}
QLabel#hintLabel {
    color: #6b4f2f;
    background: #efe0c8;
    border: 1px solid #d2bf9b;
    border-radius: 10px;
    padding: 10px;
}
QFrame#moduleCard {
    background: #fcfaf6;
    border: 1px solid #d8ceb8;
    border-radius: 18px;
    padding: 16px;
}
QLabel#cardTitle {
    font-size: 18px;
    font-weight: 700;
    color: #24313a;
}
QLabel#cardBody {
    color: #51616d;
}
QLineEdit,
QPlainTextEdit,
QListWidget,
QSpinBox,
QDoubleSpinBox {
    background: #ffffff;
    border: 1px solid #c8bda5;
    border-radius: 8px;
    padding: 6px;
}
QPushButton {
    background: #3b4954;
    color: #f7ece1;
    border: 0;
    border-radius: 10px;
    padding: 9px 14px;
    font-weight: 600;
}
QPushButton:hover {
    background: #30404c;
}
QPushButton:pressed {
    background: #24313a;
}
"""


def main() -> int:
    app = QApplication.instance()
    owns_app = app is None
    if app is None:
        app = QApplication(sys.argv)

    app.setApplicationName("chromatoPy")
    app.setStyleSheet(APP_STYLE)

    window = ChromatoPyMainWindow()
    window.show()

    if owns_app:
        return run_application(app)
    return 0
