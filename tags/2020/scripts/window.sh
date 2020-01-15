#!/bin/bash
#Script to pipe the results of a given script into a zenity window
##check that zenity is in the linux distro
#nl=$(whereis zenity | awk '{print NF}')
#if [ $nl -eq 1 ]; then
#   echo "Sorry, but it seems that zenity is not installed in this linux distribution"
#   echo "window.sh cannot be employed"
#   exit 1
#else
   if [[ $1 == *.sh ]]; then
      $* |  zenity --text-info --font='DejaVu Sans Mono' --title="$1" --width=800 --height=500   2>/dev/null
   else
      zenity --text-info --font='DejaVu Sans Mono' --title="$1" --filename=$1 --width=800 --height=500   2>/dev/null
   fi
#fi

#Add a logo into {SHARE and put is instead of foto.png}
#      zenity --text-info --font='DejaVu Sans Mono' --title="$1" --filename=$1 --width=800 --height=500   --window-icon=$HOME/foto.png  2>/dev/null

