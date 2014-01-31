#!/bin/bash

outputfile=`echo $1 | sed 's/xml/txt/g'`
lwtprint $1 -t sim_inspiral -c h_end_time > ${outputfile}
