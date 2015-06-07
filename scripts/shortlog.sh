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
  echo "================"
  echo -e $2
  echo "================"
  git log $1..$2 --reverse --no-merges --pretty=format:"commit %h%x0AAuthor: %an <%ae> %x0A%x0A %x09 %ad %s %x0A" --date=short | git shortlog -n -w130,2,13

#  git shortlog --no-merges --pretty=format:"%ad %h %s" --date=short $1..$2
else
  echo "================"
  echo -e $1
  echo "================"
  git log $1 --reverse --no-merges --pretty=format:"commit %h%x0AAuthor: %an <%ae> %x0A%x0A %x09 %ad %s %x0A" --date=short | git shortlog -n -w130,2,13

#  git shortlog --no-merges --pretty=format:"%ad %h %s" --date=short $1
fi
