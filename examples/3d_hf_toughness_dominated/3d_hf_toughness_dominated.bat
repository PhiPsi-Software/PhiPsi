@echo off

echo ************************************************************************
echo *                             VERSION 1.1                              *
echo ************************************************************************
echo *   This batch file is used to generate PhiPsi input files using ANSYS.*
echo *   This batch file is written by Fang Shi.                            *
echo *   The user should modify files and folders as the case may be.       *
echo *   Website: http://phipsi.top/                                        *
echo *   Email:   shifang@hyit.edu.cn / shifang@ustc.edu.cn                 *
echo *   Initial Date: 2021-08-01; Modification Date: 2022-05-18.           *
echo ************************************************************************

pause

@echo on

:: Set files and folders.
::---------------------     USER DEFINED PART   ----------------------------
set PhiPsi_File_Name=3d_hf_toughness_dominated
set PhiPsi_Folder=X:\PhiPsi_Work\3d_hf_toughness_dominated
if exist "C:\Program Files\ANSYS Inc\v221\ansys\bin\winx64\MAPDL.exe" (set ANSYS_Program=C:\Program Files\ANSYS Inc\v221\ansys\bin\winx64\MAPDL.exe)
if exist "D:\Program Files\ANSYS Inc\v221\ansys\bin\winx64\MAPDL.exe" (set ANSYS_Program=D:\Program Files\ANSYS Inc\v221\ansys\bin\winx64\MAPDL.exe)
if exist "E:\Program Files\ANSYS Inc\v221\ansys\bin\winx64\MAPDL.exe" (set ANSYS_Program=E:\Program Files\ANSYS Inc\v221\ansys\bin\winx64\MAPDL.exe)
set ANSYS_Work=X:\ansys_work
set OUT_FILE=X:\ansys_work\file.out
::--------------------- END OF USER DEFINED PART ---------------------------

::RUN ANSYS.
::This code line is generated according to the ANSYS Product luncher ANSYS batch, after setting, select Tool < Display Command Line.
"%ANSYS_Program%" -p ane3flds -np 4 -lch -dir "%ANSYS_Work%" -j "file" -s read -l en-us -b -i "%PhiPsi_Folder%\%PhiPsi_File_Name%.apdl" -o "%OUT_FILE%"

::COPY AND REPLACE GENERATED INPUT FILES.
if exist "%ANSYS_Work%\%PhiPsi_File_Name%.elem" ( xcopy /y "%ANSYS_Work%\%PhiPsi_File_Name%.elem" %PhiPsi_Folder%)
if exist "%ANSYS_Work%\%PhiPsi_File_Name%.node" ( xcopy /y "%ANSYS_Work%\%PhiPsi_File_Name%.node" %PhiPsi_Folder%)
if exist "%ANSYS_Work%\%PhiPsi_File_Name%.focx" ( xcopy /y "%ANSYS_Work%\%PhiPsi_File_Name%.focx" %PhiPsi_Folder%)
if exist "%ANSYS_Work%\%PhiPsi_File_Name%.focy" ( xcopy /y "%ANSYS_Work%\%PhiPsi_File_Name%.focy" %PhiPsi_Folder%)
if exist "%ANSYS_Work%\%PhiPsi_File_Name%.focz" ( xcopy /y "%ANSYS_Work%\%PhiPsi_File_Name%.focz" %PhiPsi_Folder%)
if exist "%ANSYS_Work%\%PhiPsi_File_Name%.boux" ( xcopy /y "%ANSYS_Work%\%PhiPsi_File_Name%.boux" %PhiPsi_Folder%)
if exist "%ANSYS_Work%\%PhiPsi_File_Name%.bouy" ( xcopy /y "%ANSYS_Work%\%PhiPsi_File_Name%.bouy" %PhiPsi_Folder%)
if exist "%ANSYS_Work%\%PhiPsi_File_Name%.bouz" ( xcopy /y "%ANSYS_Work%\%PhiPsi_File_Name%.bouz" %PhiPsi_Folder%)

::DELETE FILES.
if exist "%ANSYS_Work%\%PhiPsi_File_Name%.elem" ( del /f /s /q "%ANSYS_Work%\%PhiPsi_File_Name%.elem")
if exist "%ANSYS_Work%\%PhiPsi_File_Name%.node" ( del /f /s /q "%ANSYS_Work%\%PhiPsi_File_Name%.node")
if exist "%ANSYS_Work%\%PhiPsi_File_Name%.focx" ( del /f /s /q "%ANSYS_Work%\%PhiPsi_File_Name%.focx")
if exist "%ANSYS_Work%\%PhiPsi_File_Name%.focy" ( del /f /s /q "%ANSYS_Work%\%PhiPsi_File_Name%.focy")
if exist "%ANSYS_Work%\%PhiPsi_File_Name%.focz" ( del /f /s /q "%ANSYS_Work%\%PhiPsi_File_Name%.focz")
if exist "%ANSYS_Work%\%PhiPsi_File_Name%.boux" ( del /f /s /q "%ANSYS_Work%\%PhiPsi_File_Name%.boux")
if exist "%ANSYS_Work%\%PhiPsi_File_Name%.bouy" ( del /f /s /q "%ANSYS_Work%\%PhiPsi_File_Name%.bouy")
if exist "%ANSYS_Work%\%PhiPsi_File_Name%.bouz" ( del /f /s /q "%ANSYS_Work%\%PhiPsi_File_Name%.bouz")

echo ALL DONE.

