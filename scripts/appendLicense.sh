#!/bin/sh
# script to copy the headers to all the source files and header files

for f in ../README; do
  if grep -q 'Software License Agreement (New BSD License)' "$f";then 
    echo "No need to copy the License Header to $f"
  else
    cat ../License $f > $f.new
    mv $f.new $f
    echo "License Header copied to $f"
  fi 
done  

for f in ../include/*.hpp; do
  if grep -q 'Software License Agreement (New BSD License)' "$f";then 
    echo "No need to copy the License Header to $f"
  else
    cat ../License $f > $f.new
    mv $f.new $f
    echo "License Header copied to $f"
  fi 
done  

for f in ../src/*.cpp; do
  if grep -q 'Software License Agreement (New BSD License)' "$f";then 
    echo "No need to copy the License Header to $f"
  else
    cat ../License $f > $f.new
    mv $f.new $f
    echo "License Header copied to $f"
  fi 
done  