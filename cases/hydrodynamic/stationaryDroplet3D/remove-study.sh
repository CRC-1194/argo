#!/usr/bin/env bash
STUDY_PATTERN=$1

# Remove database files, we do not use them
rm *.database
rm -rf "$STUDY_PATTERN"*
