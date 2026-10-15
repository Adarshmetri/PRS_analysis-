@echo off
REM ==========================================================
REM Master Runner for Windows Command Prompt
REM This launches the pipeline using Git Bash / Bash on Windows
REM ==========================================================

where bash >nul 2>nul
if %ERRORLEVEL% EQU 0 (
    echo Launching PRS Pipeline via bash...
    bash run_pipeline.sh
) else (
    echo Bash was not found in your PATH.
    echo If you have Git installed, please run this using Git Bash:
    echo     bash run_pipeline.sh
    pause
)
