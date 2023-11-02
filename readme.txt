1.- Se ejecuta el programa conf.py 

python3 conf.py

 crear FCC              [1]
 crear foto             [2]
 STOP                   [0]
1
 nc  N=4 nc^3
3
 número de moleculas          108
 density
0.5
 Temperatura reducida
2
 reduced temperature   2.0000000000000000     
 crear FCC              [1]
 crear foto             [2]
 STOP                   [0]
0
 nombre del archivo de salida
fort.dat
 boxx,boxy,boxz       108  5.99  5.99  5.99 
 
 2.- Se corre el programa DM_python.py
 
 python3 DM_python.py
 
 3.- (Opcional)
 
 Para cambiar los valores de entrada se deben editar los documentos md.dat y param.mix
 
 El archivo md.dat tiene la siguiente estructura:
 
 Nombre del documento que contiene la configuración inicial
 Número de ciclos
 Frecuencia de impresión de los datos
 Frecuencia de la ngr
 nfilm
 Tiempo de pasos 
 Temperatura
 Radio de corte
 
 Y el archivo param.mix:
 
 Número de especies
 Naturaleza de los átomos H,O,Ar...
 Número de moléculas por especie Masa de las moléculas Sigma Eps
 
 donde:
 sigma=diámetro de las moléculas
 eps=profundidad del pozo de potencial
  
