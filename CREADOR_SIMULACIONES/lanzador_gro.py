# -*- coding: utf-8 -*-

from subprocess import call
import os

# Script sencillo que de forma automática invoca a gro y ejecuta los archivos .gro
# que se hayan creado y almacenado en la ruta que se le indique.

archivos = os.listdir("C:/experimentos/")
os.chdir('C:/')
for i in range(0,len(archivos)):
     call("C:/EXE-elegro-1-2-3-2rel/gro/gro.exe C:/experimentos/"+archivos[i])