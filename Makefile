#* **************************************************************************
# * This file is part of Timothy
# *
# * Copyright (c) 2014/15 Maciej Cytowski
# * Copyright (c) 2014/15 ICM, University of Warsaw, Poland
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# *
# * *************************************************************************/

.PHONY: all clean doc clean_doc clean_all

IS_GMAKE = $(shell command -v gmake > /dev/null; echo $$?)
IS_DOXYGEN = $(shell command -v doxygen > /dev/null; echo $$?)


ifeq "$(IS_GMAKE)" "0"
	MAKE = gmake
else
	MAKE = make
endif 

all:
	make -C src

validator:
	make -C src ../validator

timothy:
	make -C src ../timothy

clean:
	make -C src clean 

clean_doc:
	rm -rf doc/html doc/latex 

clean_all: clean clean_doc 

doc:
ifeq "$(IS_DOXYGEN)" "0"
		make -C src doc 
else 
		@echo "Doxygen not found"
endif 
