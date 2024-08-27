#!/bin/bash
#include<stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <unistd.h>
[ -e /tmp/fd1 ] || mkfifo /tmp/fd1
exec 3<>/tmp/fd1
rm -rf /tmp/fd1
for ((i=1;i<=10;i++));
do
        echo >&3
done


for i in gene_TES/*;
do
read -u3
{
  name=${i/gene_TES/gene_TES_cluster};
  if [[  -f "$name" ]]; then
    echo "$name exists"; echo >&3; continue
  fi
  echo "$name start"
  Rscript DBscan-cluster.R $i $name
  echo >&3
}&
done

wait
exec 3<&-
exec 3>&-

