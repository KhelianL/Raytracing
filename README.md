# Projet Raytracing par Khélian LARVET

# RENDER
![alt text](https://github.com/KhelianL/Raytracing/blob/main/CornellBox_[1000x1000]_[100samples]_[5shadows].png?raw=true)

# EXECUTION
- Makefile produit un executable ./main ;
- Commande rapide pour compiler et executer : make clean ; make && ./main ;
- Le render s'exécute à l'appuie de la touche "r" (attention majuscule) et le résultat s'affiche dans le fichier "rendu.ppm" ;
- "raytracingParam.h" contient tous les paramètres modifiables ;

# PHASES & AJOUTS
- Ajout        : Threads sur chaque ligne de l'image (taille statique! voir raytracingParam.h)
- Phase 1      : Toutes les intersections sont calculées (Line - Plane - Triangle - Square - Mesh)
- Phase 2.1    : Implémentation Phong
- Phase 2.2    : Ombres hard
- Phase 2.3    : Ombres soft
- Phase 3.1    : Importation Mesh (Très lent en fonction des paramètres!)
- Phase 3.2    : Reflexion / Refraction

# COMMANDES
- 'r' : render
- 'q' : quit
- 'w' : polygon mode
- '+' : next scene
