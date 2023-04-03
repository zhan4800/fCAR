%function fCAR(input_file)

load(input_file);

res = Gfmm_Cor_Car2(Y,model,wavespecs,MCMCspecs,sampleU,get_sigma,blocksize,paramroute,auto); 

save(strcat(output_path,'output.mat'),'res');
