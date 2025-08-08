#!/bin/bash
cd "$(dirname "$0")"
export RENV_AUTOLOADER_ENABLED=TRUE
export RENV_CONFIG_PAK_ENABLED=TRUE
rig rstudio ./renv.lock 
