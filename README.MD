# Projeto PrEpidemia
Simulacao do COVID-19 no Distrito Federal e RIDE no Ambiente Gama


Gama é uma plataforma de desenvolvimento de modelagem e simulação disponível no <https://gama-platform.github.io> 

Essa plataforma foi utilizado para implementar o modelo para simulação espaciotemporal da epidemia do COVID-19 no Distrito Federal e região de entorno. 

O roteiro de instalação de instalação do Gama e do modelo está no arquivo: 

`Tutorial_instalar_GAMA_e_modelo.pdf`

Os arquivos de entrada estão na pasta:  

`includes`

Ao rodar o modelo são gerados, além da saída gráfica no computador, dois arquivos texto na pasta:

` 
models
`

O arquivo de saido está no formato CSV com os indivíduos SEIR (S=Susceptivel, E=exposto, I=Infeccioso,R =Recuperado) organizado por  municipio, por dia, num període de 365 dias. 
Cada linha do arquivo representa um dia, e nas colunas temos os 4 grupos de pessoas para cada municipio e o DF. O modelo foi desenvolvido usando  equações diferenciais e a progressão espacial foi feito por modelagem baseado em agentes. 

## Parametros
Os principais parametros da simulação devem ser ajustados no programa

	float rzero <-  4.00;
	float gamma <- 1/14;
	float sigma <- 1/3;	 

Onde rzero, gamma e sigma representam respecitvamente Número básico de reprodução, Período de latência médio e Período infeccioso médio.

Para rodar uma nova simulação é necessário mudar o valor de rzero e rodar o programa do novo.

## Resultados
Foram rodadas 4 cenários com valor de rzero de 4,3, 2.5 e 2. Os resultados foram colocados num arquivo excel na pasta:

`
resultados
`
 
Cada arquivo excell tem a quantidade de indivíduos organizados em 4 subpasta S,E,I e R de todos os municipios e o DF para um período de 365 dias.

