#!/usr/bin/env python3

import os, sys, glob
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--files', help='Comma seperated list of files to process', )
args = parser.parse_args()

if args.files:
    files = [x.strip() for x in args.files.split(',')]
else:
    files = []

OKGREEN = '\033[92m'
RED = '\u001b[31;1m'

testFiles = glob.glob('tests/executables/**/*.exe', recursive=True)
failedTests = []

for testFile in testFiles:
    testName = os.path.basename(testFile).split('test_')[1].split('.exe')[0]
    if(not files or testName in files):
        print("Testing " + testFile + ":")
        completedProcess = subprocess.run([testFile])
        if(completedProcess.returncode):
            failedTests.append(testFile)

if not failedTests:
    print(OKGREEN + 'All tests completed successfully!')
else:
    for failedTest in failedTests:
        print(RED + 'Tests ' + failedTest + ' failed!')

print('\u001b[0m')
