@echo off
chcp 437
Set ProcessName=send_agent.exe
echo    Is process %ProcessName% running?
TaskList /FI "ImageName EQ %ProcessName%" 2>nul|Find /I "%ProcessName%">nul||(
echo "Process %ProcessName% is NOT running!"
)

Set ProcessName=a2009_reg.exe
echo    Is process %ProcessName% running?
TaskList /FI "ImageName EQ %ProcessName%" 2>nul|Find /I "%ProcessName%">nul||(
echo "Process %ProcessName% is NOT running!"
)

Set ProcessName=datamonitor.exe
echo    Is process %ProcessName% running?
TaskList /FI "ImageName EQ %ProcessName%" 2>nul|Find /I "%ProcessName%">nul||(
echo "Process %ProcessName% is NOT running!"
)

echo ""
echo Directory c:/angara/files contains:
dir d:\Work\seis\data\irkut\kel\2016\16_02\
