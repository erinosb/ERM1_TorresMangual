#!/bin/bash

#execute with:
# bash makeContactSheet.sh <contactsheetname.jpg>

contactsheetname=$1

# User imagemagick to make a contact sheet
montage -verbose -label '%f' -pointsize 14 -font /Library/Fonts/Arial\ Unicode.ttf -background '#000000' -fill 'gray' -geometry 300x300+2+2 -auto-orient *.jpg $contactsheetname

#-define jpeg:size=200x200