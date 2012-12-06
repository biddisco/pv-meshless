#!/bin/bash

# $1 from tag
# $2 to tag

# 05729f831f   SPHERIC-2011
# SPHERIC-2011 release-3.10
# release-3.10 SPHERIC-2012
# SPHERIC-2012 release-3.14

# echo "SPHERIC-2011" > tmp.txt && cat tmp.txt ChangeLog.txt > tmp2.txt && mv tmp2.txt ChangeLog.txt

if [ -n "$2" ]
 then 
  echo -e "\n"$2"\n" 
  git log $1..$2 --no-merges --pretty=format:"%ad %x09 %an%x09 %h %s" --date=short | gawk '{printf("%s %-24s %s ", $1, $2" "$3, $4); $1=$2=$3=$4=""; print substr($0,0,80)}'
else
  echo -e "\n"$1"\n"
  git log $1 --no-merges --pretty=format:"%ad %x09 %an%x09 %h %s" --date=short | gawk '{printf("%s %-24s %s ", $1, $2" "$3, $4); $1=$2=$3=$4=""; print substr($0,0,80)}'
fi
