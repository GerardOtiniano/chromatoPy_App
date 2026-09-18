#!/bin/zsh

set -u

PROJECT_ROOT="$(cd "$(dirname "$0")" && pwd)"
LOCAL_PYTHON="$PROJECT_ROOT/.venv/bin/python"
PYTHON_LAUNCHER="$PROJECT_ROOT/start_chromatopy.py"
QT_PLUGINS="$PROJECT_ROOT/.venv/lib/python3.11/site-packages/PySide6/Qt/plugins"

cd "$PROJECT_ROOT" || exit 1

if [[ ! -f "$PYTHON_LAUNCHER" ]]; then
    print -u2 "Could not find start_chromatopy.py in:"
    print -u2 "  $PROJECT_ROOT"
    exit 1
fi

export QT_API="pyside6"

if [[ -x "$LOCAL_PYTHON" ]]; then
    # macOS may mark Qt plugins as hidden, which prevents Qt from
    # discovering the Cocoa platform plugin.
    if [[ -d "$QT_PLUGINS" ]]; then
        chflags -R nohidden "$QT_PLUGINS"
    fi

    exec "$LOCAL_PYTHON" "$PYTHON_LAUNCHER"
fi

if command -v poetry >/dev/null 2>&1; then
    exec poetry run python "$PYTHON_LAUNCHER"
fi

print -u2 "chromatoPy's local Python environment was not found."
print -u2 "From this directory, create it once with:"
print -u2 "  poetry install"
print -u2 "Then run this launcher again."
exit 1
