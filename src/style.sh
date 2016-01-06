#!/bin/sh
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

# Default style for Timothy source files can be set by Artistic Style (astyle) tool.
# Default style is Kernighan&Ritchie style with 2 spaces indent.

if command -v astyle >/dev/null 2>&1; then
  astyle --style=kr --indent=spaces=2 *.c
  astyle --style=kr --indent=spaces=2 *.h
else 
  echo "Error: astyle not found." 
  echo "More info about astyle: http://astyle.sourceforge.net/."
  exit 1
fi

