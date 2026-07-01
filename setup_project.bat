@echo off
setlocal enabledelayedexpansion

echo.
echo ============================================================
echo   BrainSafe AI v6  --  Project Setup Script
echo   Working directory: D:\BRAINSAFE_AI
echo ============================================================
echo.

:: ── Step 1: Move to project root ────────────────────────────────────────────
cd /d D:\BRAINSAFE_AI
if %errorlevel% neq 0 (
    echo ERROR: D:\BRAINSAFE_AI does not exist.
    echo Creating it now...
    mkdir D:\BRAINSAFE_AI
    cd /d D:\BRAINSAFE_AI
)

:: ── Step 2: Create all directories ──────────────────────────────────────────
echo [1/6] Creating folder structure...
for %%d in (
    data
    models_v5
    scripts
    tests
    manuscript_final\figures
    manuscript_final\tables
    logs
) do (
    if not exist "%%d" (
        mkdir "%%d"
        echo       Created: %%d
    )
)
echo       Done.

:: ── Step 3: Check Python version ────────────────────────────────────────────
echo.
echo [2/6] Checking Python version...
python --version 2>nul
if %errorlevel% neq 0 (
    echo ERROR: Python not found. Please install Python 3.10+ from https://python.org
    pause
    exit /b 1
)

:: ── Step 4: Create virtual environment ──────────────────────────────────────
echo.
echo [3/6] Creating virtual environment (brainsafe_env)...
if not exist "brainsafe_env" (
    python -m venv brainsafe_env
    echo       Virtual environment created.
) else (
    echo       Virtual environment already exists.
)

:: ── Step 5: Activate and upgrade pip ────────────────────────────────────────
echo.
echo [4/6] Activating environment and upgrading pip...
call brainsafe_env\Scripts\activate.bat
python -m pip install --upgrade pip --quiet
echo       pip upgraded.

:: ── Step 6: Install requirements ────────────────────────────────────────────
echo.
echo [5/6] Installing dependencies (this may take 5-10 minutes)...
pip install -r requirements.txt
if %errorlevel% neq 0 (
    echo WARNING: Some packages may have failed. Check output above.
    echo Trying RDKit separately...
    pip install rdkit
)

:: ── Step 7: Clone/update repo ───────────────────────────────────────────────
echo.
echo [6/6] Checking GitHub repo...
if not exist ".git" (
    echo       No git repo found. If you want to clone fresh, run:
    echo       git clone https://github.com/krishna-g-999/brainsafe-ai.git .
) else (
    echo       Existing git repo found. To pull latest, run: git pull
)

echo.
echo ============================================================
echo   Setup complete!
echo.
echo   To activate the environment in future sessions:
echo     brainsafe_env\Scripts\activate.bat
echo.
echo   Next steps:
echo     1. python scripts\validate_all_fixes.py
echo     2. python scripts\generate_training_data.py
echo     3. python ml_v5_training.py --data data/brainsafe_training_set.csv --out models_v5/
echo     4. streamlit run app.py
echo ============================================================
pause
