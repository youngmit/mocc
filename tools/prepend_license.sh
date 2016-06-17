#!/bin/bash
# script to copy the headers to all the source files and header files
# http://unix.stackexchange.com/questions/20641/how-to-prepend-a-license-header-recursively-for-all-h-and-cpp-files-in-a-direc

printf "/*
   Copyright 2016 Mitchell Young

   Licensed under the Apache License, Version 2.0 (the \"License\");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an \"AS IS\" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

" > /tmp/$USER-license

shopt -s globstar nullglob extglob

for f in **/*.@(cpp|hpp) ;do
  if (grep Copyright $f);then 
    echo "No need to copy the License Header to $f"
  else
	[[ -f $f ]] && cat "/tmp/$USER-license" "$f" | sponge "$f"
    echo "License Header copied to $f"
  fi 
done 
