This is my development package of EOF and related routines for understanding variance in climate variables. Contributions, bug reports, suggestions, reqpests, etc are extremely welcome.

EOF requires my waveprac package as well. Download both to your machine. Then in R

In R:  

	install.package("devtools")
 	require(devtools)
	load_all('your.path.t0/waveprac')
	load_data('your.path.to/waveprac')
	load_all('your.path.to/EOF')
	load_data('your.path.to/EOF')  


The current version of EOF is uber-beta. Much needs to be done.   



