#!/usr/bin/env sh
cd $(dirname "$0")
pydoc -w ./dice.py
pydoc ./dice.py > dice.txt

