#!/bin/bash

awk -v M='1.98892e30' -v R='6.955e8' \
	'BEGIN {
		Inorm=1e7*M*R*R; 
		printf"#%24s %25s %25s %25s %25s %25s %25s\n", "star_age", "conv_mx1_bot", "conv_mx1_bot_r", "L", "R", "Iconv", "Irad";
		printf "#%24s %25s %25s %25s %25s %25s %25s\n", "[1]", "[2]", "[3]", "[4]", "[5]", "[6]", "[7]"
	}

	{
		if(NR>6 && $2>100.0 && ($67>0 || $2<3e7)) 
			printf "%25.16e %25.16e %25.16e %25.16e %25.16e %25.16e %25.16e\n", $2/1e9, $9, $17, 10.0**$42, 10.0**$43, $66/Inorm, $67/Inorm;
	}' $1
