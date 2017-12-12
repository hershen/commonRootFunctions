#!/usr/bin/env python3

import os, sys, glob
import subprocess

OKGREEN = '\033[92m'
RED = '\u001b[31;1m'

testFiles = glob.glob('tests/executables/**/*.exe')
failedTests = []

for testFile in testFiles:
    completedProcess = subprocess.run([testFile])
    if(completedProcess.returncode):
        failedTests.append(testFile)

if not failedTests:
    print(OKGREEN + 'All tests completed successfully!')
else:
    for failedTest in failedTests:
        print(RED + 'Tests ' + failedTest + ' failed!')

print('\u001b[0m')
