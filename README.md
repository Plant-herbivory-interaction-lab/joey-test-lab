# CammaranoMonardaPollinator2024
Repo for pollination visit surveys at four sites in Wisconsin conducted in summer 2022. The manuscript is under revisions at Natural Areas Journal:
Cammarano, J.H., S. Kroening, E.S. Calixto, and P.G. Hahn. Co-blooming neighbor plant diversity and floral display color similarity associated with higher flower visitation to focal species. Natural Areas Journal (in revision).

All the data are contained within the excel file 'Pollination_neighborhood_data21.xlsx'.

## Sheet 'plots': 
This sheet contains information on the quadrats (n = 20) surveyed each date.

variables	  units	            description
----------  --------          --------------
plot	      #	                1 m diameter circular quadrat w/ central focal sp plant
focal_stam	count	            no. of focal staminate flowers in plot
focal_pist	count	            no. of focal pistillate flowers in plot
grass_ht	  cm	              apparent average height of grass in plot
grass_cover	%	                apparent cover by grass in plot
temp	      °F	              temperature
sun	        categorical	      relative sunlight
wind	      Beaufort scale    value	wind speed
start 	    hr:min	          time starting site flower visitor survey
end	        hr:min	          time ending site flower visitor survey


## Sheet 'plants': 
This sheet contains information on flowering plant species within each quadrat, with number of flowers, height of tallest flower (cm), and height of shortest flower (cm).

variables	          units	                      description
----------          --------                    --------------
plot                # 	                        1 m diameter circular quadrat w/ central focal sp plant
flower_sp	          sp name                     blooming sp within plot
flower_no	          count	                      no. of inflorescences of sp in plot
flower_tall	        cm  	                      height of tallest peduncle of sp in plot
flowe_short	        cm	                        height of shortest peduncle of sp in plot
[neighbor sp name]  sequence no. (flower no.)   order of visit to neighbor plant inflorescence and number of flowers visited

## Sheet 'visitors': 
This sheet contains information on pollinator vistors to each quadrat, the species they visited, and the number of visits.

variables	      units	      description
----------      --------    --------------
plot	          # 	        1 m diameter circular quadrat w/ central focal sp plant
visitor_sp      sp name	    insect visiting flowers within plot during 10-minute timed period
visitor_no	    count	      no. of individuals of insect sp visiting flowers within plot during timed period
flower_sp	      sp name	    blooming sp visited by insect
flower_visits	  count	      no. of inflorescences visited by insect
switch	        count	      change in sp being visited by insect
flower_sp_2	sp  name	      blooming sp switched to
flower_visits_2	count	      no. of inflorescences visited by insect of sp switched to


## Sheet 'plant_spp_traits': 
This sheet contains information on floral traits of each species, along with the references from where they were obtained.

variables	      units	      description
----------      --------    --------------
petal_color	    categorical	color of petal or colored sepal
petal_wvl	      scale	      violet=1, blue=2, green=3, yellow=4, orange=5, red=6
corolla_depth	  mm	        length of corolla/floret from base to tip of petal
flor_per_inflor	count	      average flowers per inflorescence
inflor_diam	    mm	        diameter or length of inflorescence
orientation	    categorical	orientation of floral surface
references		              source of information


## Sheet 'visitor_spp_traits': 
This sheet contains incomplete information on pollinator species traits.



