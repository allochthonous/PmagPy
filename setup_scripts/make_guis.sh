#!/bin/bash

# usage: ./make_guis.sh "new_commit_name"
# removes all pre-existing dist/ stuff
# makes new frozen binaries
# commits them to PmagPy-Standalone-OSX
# pushes the new commit to Github

pyv="$(which python)"
echo $pyv
if [ $pyv == "/usr/local/bin/python" ]; then
    echo 'Using correct Python distribution'
elif [[ $pyv == *"Canopy"* ]]; then
    echo "-W- Using incorrect Python distribution (py2app can't run with Canopy)"
    echo "-W- Please try again with a non-Canopy distribution"
    exit
else
    echo "-W- Using unknown Python distribution.  This may or may not work!"
fi

new="$1"

if [ -z "$new" ]; then
    echo "-W- You need to provide a commit name"
    echo 'Correct usage is: make_guis.sh "new commit name"'
    exit
fi

my_commit="$new"
echo $my_commit

# change to PmagPy main directory

echo "starting make_guis script"
cp setup_scripts/setup_pmag_gui.py .
#cp setup_scripts/setup_magic_gui2.py .
cp setup_scripts/setup_magic_gui.py .
echo "copied in setup scripts"
rm -rf build dist
echo "removed build & dist"

python setup_pmag_gui.py py2app
echo "ran setup script for Pmag GUI"
rm -rf ../PmagPy-Standalone-OSX/pmag_gui.app
echo "deleted old Pmag GUI"
mv dist/pmag_gui.app ../PmagPy-Standalone-OSX
echo "moved pmag_gui.app to PmagPy-Standalone-OSX"

rm -rf build dist
echo "removed build & dist"

#python setup_magic_gui2.py py2app
#echo "ran setup script for Magic GUI 2"
#rm -rf ../PmagPy-Standalone-OSX/magic_gui2.app
#echo "deleted old Magic GUI"
#mv dist/magic_gui2.app ../PmagPy-Standalone-OSX
#echo "moved magic_gui2.app " #to PmagPy-Standalone-OSX"

#rm -rf build dist
#echo "removed build & dist"

python setup_magic_gui.py py2app
echo "ran setup script for MagIC GUI"
rm -rf ../PmagPy-Standalone-OSX/magic_gui.app
echo "deleted old Pmag GUI"
mv dist/magic_gui.app ../PmagPy-Standalone-OSX
echo "moved magic_gui.app to PmagPy-Standalone-OSX"

echo "clean up"
rm setup_pmag_gui.py
#rm setup_magic_gui2.py
rm setup_magic_gui.py

cd ../PmagPy-Standalone-OSX

echo 'adding git stuff'
git add .
echo 'adding .'
git commit -m "$my_commit"
echo 'committed'


echo 'Done!  Changes have been committed to PmagPy-Standalone-OSX'
