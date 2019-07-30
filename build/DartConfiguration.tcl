# This file is configured by CMake automatically as DartConfiguration.tcl
# If you choose not to use CMake, this file may be hand configured, by
# filling in the required variables.


# Configuration directories and files
SourceDirectory: /afs/desy.de/user/s/sdharani/flc/eeZHqqbb/Hbb
BuildDirectory: /afs/desy.de/user/s/sdharani/flc/eeZHqqbb/Hbb/build

# Where to place the cost data store
CostDataFile: 

# Site is something like machine.domain, i.e. pragmatic.crd
Site: naf-ilc11.desy.de

# Build name is osname-revision-compiler, i.e. Linux-2.4.2-2smp-c++
BuildName: Linux-g++

# Submission information
IsCDash: 
CDashVersion: 
QueryCDashVersion: 
DropSite: 
DropLocation: 
DropSiteUser: 
DropSitePassword: 
DropSiteMode: 
DropMethod: http
TriggerSite: 
ScpCommand: /usr/bin/scp

# Dashboard start time
NightlyStartTime: 00:00:00 EDT

# Commands for the build/test/submit cycle
ConfigureCommand: "/cvmfs/ilc.desy.de/sw/x86_64_gcc49_sl6/CMake/3.6.3/bin/cmake" "/afs/desy.de/user/s/sdharani/flc/eeZHqqbb/Hbb"
MakeCommand: /cvmfs/ilc.desy.de/sw/x86_64_gcc49_sl6/CMake/3.6.3/bin/cmake --build . --config "${CTEST_CONFIGURATION_TYPE}" -- -i
DefaultCTestConfigurationType: Release

# version control
UpdateVersionOnly: 

# CVS options
# Default is "-d -P -A"
CVSCommand: /usr/bin/cvs
CVSUpdateOptions: -d -A -P

# Subversion options
SVNCommand: /usr/bin/svn
SVNOptions: 
SVNUpdateOptions: 

# Git options
GITCommand: /usr/bin/git
GITInitSubmodules: 
GITUpdateOptions: 
GITUpdateCustom: 

# Perforce options
P4Command: P4COMMAND-NOTFOUND
P4Client: 
P4Options: 
P4UpdateOptions: 
P4UpdateCustom: 

# Generic update command
UpdateCommand: 
UpdateOptions: 
UpdateType: 

# Compiler info
Compiler: /cvmfs/sft.cern.ch/lcg/contrib/gcc/4.9.3/x86_64-slc6/bin/g++
CompilerVersion: 4.9.3

# Dynamic analysis (MemCheck)
PurifyCommand: 
ValgrindCommand: 
ValgrindCommandOptions: 
MemoryCheckType: 
MemoryCheckSanitizerOptions: 
MemoryCheckCommand: /usr/bin/valgrind
MemoryCheckCommandOptions: 
MemoryCheckSuppressionFile: 

# Coverage
CoverageCommand: /cvmfs/sft.cern.ch/lcg/contrib/gcc/4.9.3/x86_64-slc6/bin/gcov
CoverageExtraFlags: -l

# Cluster commands
SlurmBatchCommand: SLURM_SBATCH_COMMAND-NOTFOUND
SlurmRunCommand: SLURM_SRUN_COMMAND-NOTFOUND

# Testing options
# TimeOut is the amount of time in seconds to wait for processes
# to complete during testing.  After TimeOut seconds, the
# process will be summarily terminated.
# Currently set to 25 minutes
TimeOut: 1500

# During parallel testing CTest will not start a new test if doing
# so would cause the system load to exceed this value.
TestLoad: 

UseLaunchers: 
CurlOptions: 
# warning, if you add new options here that have to do with submit,
# you have to update cmCTestSubmitCommand.cxx

# For CTest submissions that timeout, these options
# specify behavior for retrying the submission
CTestSubmitRetryDelay: 5
CTestSubmitRetryCount: 3
