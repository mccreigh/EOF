This is my development package of EOF and related routines for understanding variance in climate variables. Contributions, bug reports, suggestions, reqpests, etc are extremely welcome.

EOF requires my waveprac package as well (https://github.com/mccreigh/waveprac). Download both to your machine. The simplest way is to click on the ZIP download above (if you dont plan to contribute, else you'll need to fork it).

In R:  

	install.package("devtools")
 	require(devtools)
	load_all('your.path.to/waveprac')
	load_data('your.path.to/waveprac')
	load_all('your.path.to/EOF')
	load_data('your.path.to/EOF')  


You will likely get some messages about required packges. Please install those from CRAN. 

The current version of EOF is uber-beta. Much needs to be done.   



