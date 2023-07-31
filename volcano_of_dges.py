# import packages
import pandas as pd
from bioinfokit import visuz
# import the DGE table (condition_infected_vs_control_dge.csv)
df = pd.read_csv("/home/scbb/abhijit/btp_4vsbtp_5_condition_infected_vs_control_dge.csv")
# drop NA values
df=df.dropna()

# create volcano plot
visuz.GeneExpression.volcano(df=df, lfc='logFC', pv='FDR', sign_line=True, plotlegend=True)
