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

	init {
		create Municipio from: municipios_shp_file with: [h::0.1, N::500, I::1.0];
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

}

species Municipio {
	float t;
	int N;
	float S <- N - I;
	float E <- 0.0;
	float I;
	float R <- 0.0;
	float h;
	float beta <- 0.4;
	float gamma <- 0.01;
	float sigma <- 0.05;
	float mu <- 0.01;
	string GID_2;
	rgb mycolor -> {hsb(0, (I > 25?0.1:0)+(I > 25 ? 25 : I) / 29, 1)};
//	rgb mycolor -> {hsb(0, I/N, 1)};

	// must be followed with exact order S, E, I, R, t  and N,beta,gamma,sigma,mu
	equation eqSEIR type: SEIR vars: [S, E, I, R, t] params: [N, beta, gamma, sigma, mu];
	list<Municipio> neighbours <- [];

	reflex solving {
		solve eqSEIR method: "rk4" step_size: h;
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
		display "Municipios" synchronized:true{
			species Municipio;
		}

		display "Statistic" {
			chart 'SEIR' type: series {
				data "S" value: sum(Municipio collect each.S) color: #green;
				data "E" value: sum(Municipio collect each.E) color: #yellow;
				data "I" value: sum(Municipio collect each.I) color: #red;
				data "R" value: sum(Municipio collect each.R) color: #blue;
			}

		}

	}

}