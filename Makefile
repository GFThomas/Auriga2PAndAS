all: density mask foreground down_res down_cat
		mkdir results
		chmod +x exec.py

down_res:
		cd results
		wget "https://ws-cadc.canfar.net/vault/synctransTARGET=vos%3A%2F%2Fcadc.nrc.ca%21vault%2Fgthomas%2FAuriga2PAndAS%2Fresults.tar.gz&DIRECTION=pullFromVoSpace&PROTOCOL=ivo%3A%2F%2Fivoa.net%2Fvospace%2Fcore%23httpget"
		cd ..

down_cat:
		cd catalogues
		wget "https://ws-cadc.canfar.net/vault/synctrans?TARGET=vos%3A%2F%2Fcadc.nrc.ca%21vault%2Fgthomas%2FAuriga2PAndAS%2Fcatalogues.tar.gz&DIRECTION=pullFromVoSpace&PROTOCOL=ivo%3A%2F%2Fivoa.net%2Fvospace%2Fcore%23httpget"
		cd ..

density:
		gcc -fPIC -fopenmp -shared ./pipeline/density.c -o ./pipeline/density.so -lm
		

mask:
		gcc -fPIC -fopenmp -shared ./pipeline/mask_footprint.c -o ./pipeline/mask_footprint.so -lm

foreground:
		gcc ./pipeline/foreground/foreground.c -o ./pipeline/foreground/foreground.x -lm

