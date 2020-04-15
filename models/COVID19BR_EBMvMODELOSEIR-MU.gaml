/***
* Name: Municipios
* Author: hqngh
* Adaptado por: Augusto Arcela
* Description: 
* Tags: Tag1, Tag2, TagN
***/
model Municipios

global {
	shape_file municipios_shp_file <- shape_file("../includes/gadm36_RIDE_shp/GADM36_BRA2_RIDE.shp");
	file pop_csv_file <- csv_file("../includes/pop_ride.csv");
	geometry shape <- envelope(municipios_shp_file);
	float time_simulation <- 365.0;
	init {
		create Municipio from: municipios_shp_file with: [h::0.01, N::3015268, I::3.0];
		ask Municipio {
			neighbours <- Municipio where (each touches self);
		}

		matrix data <- matrix(pop_csv_file.contents);
		loop i from: 1 to: data.rows - 1 {
			Municipio p <- first(Municipio where (each.GID_2 = data[0, i]));
			ask p {
				N <- int(data[1, i]);
				I <- float(data[2, i]);
				S <- N - I;
				rgb null <- mycolor;
			}
		}
	}
	
	reflex pause {
		if (time > time_simulation) {
			do pause;	
		}
	}
	// imprimindo cabecalho
    reflex save_data { save string(cycle) +   " , S= ," + 
    	string(Municipio collect each.GID_2) + " , E= ," + 
        string(Municipio collect each.GID_2) + " , I= ," +    	
        string(Municipio collect each.GID_2) + " , R= ," +    	
        string(Municipio collect each.GID_2)   
    	to: "cabecalho.txt" type: "text"; } 	
    	
}

species Municipio {
	float t <- 1 #d; 
	int N;
	float S <- N - I;
	float E <- 0.0;
	float I;
	float R <- 0.0;
	float h <- 0.01;
	float rzero <-  2.0;
	float gamma <- 1/14;
	float sigma <- 1/3;	
	float beta <- rzero * gamma;
//	float mu <- 0.001;
	string GID_2;
	rgb mycolor -> {hsb(0, (I > 25?0.1:0)+(I > 25 ? 25 : I) / 29, 1)};
//	rgb mycolor -> {hsb(0, I/N, 1)};

	// must be followed with exact order S, E, I, R, t  and N,beta,gamma,sigma,mu
equation classicSEIR {
	diff(S,t) = (- beta * S * I / N);
	diff(E,t) = (beta * S * I / N) - (sigma * E) ;
	diff(I,t) = (sigma * E) - (gamma * I);
	diff(R,t) = (gamma * I);
    }
	list<Municipio> neighbours <- [];

	reflex solving {
		solve classicSEIR method: "rk4" step_size: h;
	}

	reflex transmission {
		Municipio candi <- any(neighbours);
		if (candi != nil) {
			int r <- rnd(2);
			switch r {
				match 0 {
					if (N > 1) and (S > 1) {
						N <- N - 1;
						candi.N <- candi.N + 1;
						S <- S - 1;
						candi.S <- candi.S + 1;
					}
				}
				match 1 {
					if (N > 1) and (E > 1) {
						N <- N - 1;
						candi.N <- candi.N + 1;
						E <- E - 1;
						candi.E <- candi.E + 1;
					}
				}
			}
		}
	}

	aspect default {
		draw shape color: mycolor border: #black;
	}

}


experiment Pandemic2020 type: gui {
	output {
			layout horizontal([0::5000, 1::5000]) tabs: true editors: false;
//		layout horizontal([0::5000, 1::5000]) tabs: true editors: false;
		display "Municipios" synchronized:true autosave: true {
			species Municipio;
		}

		display "Statistic" synchronized:true autosave: true {
			chart 'SEIR' type: series size: {1,1} {
				data "S" value: sum(Municipio collect each.S) color: #green;
				data "E" value: sum(Municipio collect each.E) color: #yellow;
				data "I" value: sum(Municipio collect each.I) color: #red;
				data "R" value: sum(Municipio collect each.R) color: #blue;
			}
		}
		
		monitor  total_S value: sum(Municipio collect each.S) refresh: every(10#cycle);
		monitor  total_E value: sum(Municipio collect each.E) refresh: every(10#cycle);
		monitor  total_I value: sum(Municipio collect each.I) refresh: every(10#cycle);
		monitor  total_R value: sum(Municipio collect each.R) refresh: every(10#cycle);
		monitor  total_N value: sum(Municipio collect each.N) refresh: every(10#cycle);

		monitor  Nome_mun value: (Municipio collect each.GID_2) refresh: every(10#cycle);
		monitor  Infected_mun value: (Municipio collect each.I) refresh: every(10#cycle);

	}
	
	// Imprimindo todos os dados 

    reflex save_data { save string(cycle) + " , S= ," +  
    	string(Municipio collect each.S)  + " , E= ," +
    	string(Municipio collect each.E)  + " , I= ," +
    	string(Municipio collect each.I)  + " , R= ," +
    	string(Municipio collect each.R)      	    	
    	to: "Dados_SEIR_ro2.txt" type: "text" rewrite: false ;
    }
}




