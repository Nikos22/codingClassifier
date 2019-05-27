#!/usr/local/bin/Rscript
args = commandArgs(TRUE)

df_ss_cons = read.table(args[1],
	row.names=1,
	header=TRUE,
	stringsAsFactors = FALSE)

df = data.frame(
row.names = rownames(df_ss_cons),  
df_ss_cons[,c(1:64)],
Length = df_ss_cons$GENELENGTH, 
CAI = df_ss_cons$CAICLASS, 
GC= df_ss_cons$GC_CONT, 
GC3 = df_ss_cons$GC_CONT3,  
Cost = df_ss_cons$Cost,
perTM = (df_ss_cons$Exp_AA_in_transm)/(df_ss_cons$GENELENGTH/3),
disord = (df_ss_cons$Disordered)/(df_ss_cons$GENELENGTH/3),
Low_comp = (df_ss_cons$Complexity)/(df_ss_cons$GENELENGTH/3),
Gravy = df_ss_cons$Gravy,
Aromo = df_ss_cons$Aromo,
Helix = df_ss_cons$Helix/(df_ss_cons$GENELENGTH/3),
Sheet = df_ss_cons$Sheet/(df_ss_cons$GENELENGTH/3),
Aggregation = df_ss_cons$Aggregation/(df_ss_cons$GENELENGTH/3))

write.table(df, file=args[2], row.names = TRUE, col.names = TRUE, quote = FALSE)
