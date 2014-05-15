all: graphs.png

graphs.png: LineasEmision.py energy_counts.dat
	python LineasEmision.py

energy_counts.dat: 
	wget https:/raw.githubusercontent.com/forero/ComputationalMethodsData/master/homework/hw_8/lotka_volterra_obs.dat